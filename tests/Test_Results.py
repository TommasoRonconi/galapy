#------------------------------------------------------------------------------#
#                 Test of galapy.sampling.Results module.                       #
#------------------------------------------------------------------------------#
"""
Tests for the derived-property machinery of galapy.sampling.Results.

Covers:
- default derived quantities are all computed, tracked in ``_derived`` and
  exposed as attributes with the expected shapes (scalar -> (N,), SED -> 2D)
- add_property with a single callable + explicit name
- add_property with a single callable and automatic ``custom{N}`` naming
- add_property with a dict of callables
- add_property with an array-valued, finite quantity (linear avg attenuation)
  reachable through get_mean / get_quantile
- add_property input validation (non-callable, bad name type)
- dump / load round-trip preserves the custom derived quantities
"""

import numpy as np
import pytest

from galapy.sampling import Run
from galapy.sampling.Results import Results

# ---------------------------------------------------------------------------
# Shared test data (mirrors tests/Test_Run.py)
# ---------------------------------------------------------------------------

BANDS  = ['GOODS.b', 'GOODS.v', 'GOODS.i', 'GOODS.z']
FLUXES = np.array([1.0e-4, 3.0e-4, 9.0e-4, 2.0e-3])
ERRORS = FLUXES * 0.1
UPLIMS = np.zeros(len(BANDS), dtype=bool)

BASE_PARAMS = {
    'age'            : ([6., 11.], True),
    'redshift'       : 1.0,
    'sfh.psi_max'    : ([0.,  4.], True),
    'sfh.tau_star'   : ([6., 11.], True),
    'sfh.tau_quench' : 2e20,
    'ism.f_MC'       : 0.5,
    'ism.norm_MC'    : 100.,
    'ism.N_MC'       : ([0.,  5.], True),
    'ism.R_MC'       : ([0.,  5.], True),
    'ism.tau_esc'    : ([4.,  8.], True),
    'ism.dMClow'     : 1.3,
    'ism.dMCupp'     : 1.6,
    'ism.norm_DD'    : 1.0,
    'ism.Rdust'      : ([0.,  5.], True),
    'ism.f_PAH'      : 0.2,
    'ism.dDDlow'     : 0.7,
    'ism.dDDupp'     : 2.0,
}

DEFAULT_DERIVED = [ 'SED', 'Mstar', 'Mdust', 'Mgas',
                    'Zstar', 'Zgas', 'SFR', 'TMC', 'TDD' ]

NSAMP = 6


@pytest.fixture( scope = 'module' )
def results () :
    """Build a small but valid Results instance from a known-good model point.

    The samples are taken at the centre of the prior (a configuration that is
    guaranteed to evaluate without raising ``RuntimeError``), tiled NSAMP times
    so that every sample is finite and the machinery can be exercised
    deterministically.
    """
    state = Run.initialize( BANDS, FLUXES, ERRORS, UPLIMS,
                            BANDS, dict( BASE_PARAMS ) )
    par0       = state.handler.par_prior.mean( axis = 1 )
    sample_res = np.tile( par0, ( NSAMP, 1 ) )
    return Results( state.model, state.handler,
                    sample_res, np.zeros( NSAMP ),
                    sample_weights = np.ones( NSAMP ),
                    data = state.data, noise = state.noise,
                    sampler_name = 'dynesty' )


# ---------------------------------------------------------------------------
# Default derived quantities
# ---------------------------------------------------------------------------

class TestDefaultProperties :

    def test_all_defaults_tracked ( self, results ) :
        assert results._derived == DEFAULT_DERIVED

    def test_defaults_are_attributes ( self, results ) :
        for k in DEFAULT_DERIVED :
            assert hasattr( results, k )

    def test_scalar_shapes ( self, results ) :
        for k in [ 'Mstar', 'Mdust', 'Mgas', 'Zstar', 'Zgas', 'SFR', 'TMC', 'TDD' ] :
            assert getattr( results, k ).shape == ( NSAMP, )

    def test_sed_is_2d ( self, results ) :
        assert results.SED.ndim == 2
        assert results.SED.shape[0] == NSAMP

    def test_defaults_are_finite ( self, results ) :
        # the prior-centre configuration is valid, so nothing is sentinel'd
        for k in DEFAULT_DERIVED :
            assert np.isfinite( getattr( results, k ) ).all()

    def test_temperatures_are_positive ( self, results ) :
        # TMC/TDD are only meaningful once get_emission (via SED) has run;
        # a positive temperature confirms the ordering dependency holds
        assert ( results.TMC > 0. ).all()
        assert ( results.TDD > 0. ).all()


# ---------------------------------------------------------------------------
# add_property
# ---------------------------------------------------------------------------

class TestAddProperty :

    def test_single_callable_named ( self, results ) :
        res = _fresh( results )
        added = res.add_property( lambda model : model.age, name = 'theage' )
        assert added == [ 'theage' ]
        assert 'theage' in res._derived
        assert res.theage.shape == ( NSAMP, )
        assert np.isfinite( res.theage ).all()

    def test_single_callable_autoname ( self, results ) :
        res = _fresh( results )
        res.add_property( lambda model : model.age )
        res.add_property( lambda model : model.sfh.Mstar( model.age ) )
        assert hasattr( res, 'custom0' ) and 'custom0' in res._derived
        assert hasattr( res, 'custom1' ) and 'custom1' in res._derived

    def test_dict_of_callables ( self, results ) :
        res = _fresh( results )
        added = res.add_property( {
            'a' : lambda model : model.age,
            'b' : lambda model : model.sfh.Mstar( model.age ),
        } )
        assert set( added ) == { 'a', 'b' }
        assert res.a.shape == ( NSAMP, )
        assert res.b.shape == ( NSAMP, )

    def test_array_valued_finite_attenuation ( self, results ) :
        res = _fresh( results )

        def avg_attenuation ( model ) :
            # linear average attenuation in [0,1]: finite by construction,
            # unlike the magnitude form returned by get_avgAtt()
            _ = model.get_emission( store_attenuation = True )
            return model.Aavg

        res.add_property( avg_attenuation, name = 'Aavg' )
        assert res.Aavg.ndim == 2
        assert res.Aavg.shape[0] == NSAMP
        assert np.isfinite( res.Aavg ).all()
        # reachable through the statistics helpers
        assert np.isfinite( res.get_mean( 'Aavg' ) ).all()
        assert np.isfinite( res.get_quantile( 'Aavg' ) ).all()

    def test_mean_of_scalar_custom ( self, results ) :
        res = _fresh( results )
        res.add_property( lambda model : model.age, name = 'theage' )
        # identical samples => weighted mean equals the common value
        assert np.isfinite( res.get_mean( 'theage' ) )
        np.testing.assert_allclose( res.get_mean( 'theage' ), res.theage[0] )

    def test_non_callable_raises ( self, results ) :
        res = _fresh( results )
        with pytest.raises( AttributeError ) :
            res.add_property( { 'bad' : 42 } )

    def test_bad_name_type_raises ( self, results ) :
        res = _fresh( results )
        with pytest.raises( TypeError ) :
            res.add_property( lambda model : model.age, name = 123 )


# ---------------------------------------------------------------------------
# Serialisation round-trip
# ---------------------------------------------------------------------------

class TestDumpLoad :

    def test_defaults_roundtrip ( self, results ) :
        ret = Results.load( results.dump() )
        assert ret._derived == results._derived
        for k in DEFAULT_DERIVED :
            np.testing.assert_array_equal( getattr( ret, k ), getattr( results, k ) )

    def test_custom_property_survives_roundtrip ( self, results ) :
        res = _fresh( results )
        res.add_property( lambda model : model.age, name = 'theage' )
        ret = Results.load( res.dump() )
        assert 'theage' in ret._derived
        np.testing.assert_array_equal( ret.theage, res.theage )

    def test_loaded_object_stats_work ( self, results ) :
        res = _fresh( results )
        res.add_property( lambda model : model.age, name = 'theage' )
        ret = Results.load( res.dump() )
        np.testing.assert_allclose( ret.get_mean( 'theage' ),
                                    res.get_mean( 'theage' ) )


# ---------------------------------------------------------------------------
# helper: a fresh, independent Results so per-test add_property calls do not
# leak into the module-scoped fixture
# ---------------------------------------------------------------------------

def _fresh ( results ) :
    return Results.load( results.dump() )
