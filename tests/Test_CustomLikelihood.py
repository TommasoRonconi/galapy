#------------------------------------------------------------------------------#
#            Test of the custom log-likelihood machinery (galapy-fit).          #
#------------------------------------------------------------------------------#
"""
Tests for the pluggable ``loglikelihood`` hyper-parameter and its supporting
machinery.

Covers:
- Run._resolve_loglikelihood: None -> built-in, TypeError on non-callable or
  wrong signature, warnings on missing **kwargs and on unpicklable callables
  (lambda, closure, parameter-file namespace); functools.partial with bound
  per-object data accepted cleanly, unbound required data rejected, and the
  picklability checks see through partial wrappers
- Run.loglikelihood_name: provenance marker for functions and callable
  instances (stable, no memory address), partial wrappers unwrapped,
  consistency with the literal stored in Results._default_loglikelihood_name
- PipelineState: carries the resolved likelihood; validation at construction;
  Run.logprob dispatches to the custom likelihood while keeping the prior gate
- Run._expand_hyperpar: run-wide propagation to jobs, absent attribute ->
  None, validation in the main process, per-variant override -> ValueError
- Run._nautilus_loglikelihood: adapter restores the positional
  ``(par, state, **kwargs)`` contract under nautilus' functools.partial
  wrapping, is picklable, and matches the built-in likelihood exactly
- Run._sample_parallel: the sampler dispatch of the parallel path accepts
  nautilus (regression: it used to abort), forwards the pool and the state
  untouched, and still rejects unknown sampler names
- Results: loglikelihood_name attribute defaults to the built-in marker,
  survives dump/load, files without the key load as built-in
- model_comparison.bayes_factor: accepts raw floats and Results-like objects,
  warns on mismatched or custom likelihoods, raises on missing evidence
"""

import functools
import pickle
import warnings as _warnings

import numpy as np
import pytest

from types import SimpleNamespace

from galapy.sampling import Run
from galapy.sampling.Run import (
    PipelineState,
    _expand_hyperpar,
    _resolve_loglikelihood,
    _nautilus_loglikelihood,
    loglikelihood_name,
)
from galapy.sampling.Results import Results, _default_loglikelihood_name
from galapy.analysis.model_comparison import bayes_factor

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

NSAMP = 4


# ---------------------------------------------------------------------------
# Module-level custom likelihoods (importable, hence picklable by reference)
# ---------------------------------------------------------------------------

def constant_loglikelihood ( par, state, **kwargs ) :
    """A valid custom likelihood returning a recognisable constant."""
    return 42.


def loglikelihood_no_kwargs ( par, state ) :
    """Valid signature but missing **kwargs: should trigger a warning."""
    return 0.


def loglikelihood_wrong_signature ( par ) :
    """Cannot be called as f(par, state): should raise TypeError."""
    return 0.


def loglikelihood_needing_data ( par, state, logmstar_obs, logmstar_err,
                                 **kwargs ) :
    """Valid likelihood requiring per-object data, to be bound with
    functools.partial in the parameter file (the pattern documented in the
    custom-likelihood guide)."""
    return 0.


def _make_closure () :
    def inner ( par, state, **kwargs ) :
        return 0.
    return inner


class RecordingLoglikelihood :
    """Callable instance recording how it was invoked."""

    def __init__ ( self ) :
        self.calls = []

    def __call__ ( self, par, state, **kwargs ) :
        self.calls.append( ( par, state, kwargs ) )
        return 7.


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture( scope = 'module' )
def state () :
    """A default PipelineState (built-in likelihood)."""
    return Run.initialize( BANDS, FLUXES, ERRORS, UPLIMS,
                           BANDS, dict( BASE_PARAMS ) )


@pytest.fixture( scope = 'module' )
def custom_state () :
    """A PipelineState carrying a custom likelihood."""
    return Run.initialize( BANDS, FLUXES, ERRORS, UPLIMS,
                           BANDS, dict( BASE_PARAMS ),
                           loglikelihood = constant_loglikelihood )


@pytest.fixture( scope = 'module' )
def results ( state ) :
    """A small Results instance built without an explicit likelihood marker."""
    par0       = state.handler.par_prior.mean( axis = 1 )
    sample_res = np.tile( par0, ( NSAMP, 1 ) )
    return Results( state.model, state.handler,
                    sample_res, np.zeros( NSAMP ),
                    sample_weights = np.ones( NSAMP ),
                    data = state.data, noise = state.noise,
                    sampler_name = 'dynesty' )


# ===========================================================================
# _resolve_loglikelihood
# ===========================================================================

class TestResolveLoglikelihood :

    def test_none_selects_builtin ( self ) :
        assert _resolve_loglikelihood( None ) is Run.loglikelihood

    def test_valid_function_returned_unchanged_without_warnings ( self ) :
        with _warnings.catch_warnings( record = True ) as record :
            _warnings.simplefilter( 'always' )
            resolved = _resolve_loglikelihood( constant_loglikelihood )
        assert resolved is constant_loglikelihood
        assert len( record ) == 0

    def test_non_callable_raises ( self ) :
        with pytest.raises( TypeError, match = 'must be either None' ) :
            _resolve_loglikelihood( 42 )

    def test_wrong_signature_raises ( self ) :
        with pytest.raises( TypeError, match = 'cannot be called as' ) :
            _resolve_loglikelihood( loglikelihood_wrong_signature )

    def test_missing_kwargs_warns ( self ) :
        with pytest.warns( UserWarning, match = r'does not accept \*\*kwargs' ) :
            resolved = _resolve_loglikelihood( loglikelihood_no_kwargs )
        # still usable, only warned about
        assert resolved is loglikelihood_no_kwargs

    def test_lambda_warns_not_picklable ( self ) :
        with pytest.warns( UserWarning, match = 'lambda' ) :
            _resolve_loglikelihood( lambda par, state, **kwargs : 0. )

    def test_closure_warns_not_picklable ( self ) :
        with pytest.warns( UserWarning, match = 'closure' ) :
            _resolve_loglikelihood( _make_closure() )

    def test_parameter_file_namespace_warns_not_picklable ( self ) :
        # functions defined inside the parameter file are imported under the
        # throw-away module name "hyper_parameters" (see Run._run)
        def fake ( par, state, **kwargs ) :
            return 0.
        fake.__module__   = 'hyper_parameters'
        fake.__qualname__ = 'fake'
        with pytest.warns( UserWarning, match = 'cannot import' ) :
            _resolve_loglikelihood( fake )

    def test_callable_instance_accepted ( self ) :
        with _warnings.catch_warnings( record = True ) as record :
            _warnings.simplefilter( 'always' )
            resolved = _resolve_loglikelihood( RecordingLoglikelihood() )
        assert callable( resolved )
        assert len( record ) == 0

    def test_partial_bound_data_accepted_without_warnings ( self ) :
        # the documented per-object pattern: data bound in the parameter file
        bound = functools.partial( loglikelihood_needing_data,
                                   logmstar_obs = 10.65,
                                   logmstar_err = 0.15 )
        with _warnings.catch_warnings( record = True ) as record :
            _warnings.simplefilter( 'always' )
            resolved = _resolve_loglikelihood( bound )
        assert resolved is bound
        assert len( record ) == 0

    def test_unbound_required_data_raises ( self ) :
        # forgetting the functools.partial binding must fail at load time
        with pytest.raises( TypeError, match = 'cannot be called as' ) :
            _resolve_loglikelihood( loglikelihood_needing_data )

    def test_partial_of_lambda_still_warns_not_picklable ( self ) :
        # the picklability checks must see through the partial wrapper
        bound = functools.partial( lambda par, state, x, **kwargs : 0., x = 1. )
        with pytest.warns( UserWarning, match = 'lambda' ) :
            _resolve_loglikelihood( bound )


# ===========================================================================
# loglikelihood_name
# ===========================================================================

class TestLoglikelihoodName :

    def test_builtin_matches_results_literal ( self ) :
        # guards the literal kept in Results.py (to avoid a circular import)
        # against drifting away from the real qualified name
        assert loglikelihood_name( Run.loglikelihood ) \
            == _default_loglikelihood_name

    def test_module_level_function ( self ) :
        name = loglikelihood_name( constant_loglikelihood )
        assert name.endswith( '.constant_loglikelihood' )

    def test_callable_instance_is_stable ( self ) :
        # two instances of the same class must produce the same marker
        # (repr would embed the memory address)
        name1 = loglikelihood_name( RecordingLoglikelihood() )
        name2 = loglikelihood_name( RecordingLoglikelihood() )
        assert name1 == name2
        assert '0x' not in name1
        assert name1.endswith( '.RecordingLoglikelihood' )

    def test_partial_unwrapped_to_wrapped_function ( self ) :
        # the marker must identify the function, not the partial wrapper,
        # so per-object bindings share the provenance of their module
        bound = functools.partial( loglikelihood_needing_data,
                                   logmstar_obs = 10.65,
                                   logmstar_err = 0.15 )
        name = loglikelihood_name( bound )
        assert name == loglikelihood_name( loglikelihood_needing_data )
        assert name.endswith( '.loglikelihood_needing_data' )
        # rebinding different data keeps the same marker
        rebound = functools.partial( bound, logmstar_obs = 9.8 )
        assert loglikelihood_name( rebound ) == name


# ===========================================================================
# PipelineState carries the likelihood
# ===========================================================================

class TestPipelineStateLikelihood :

    def test_default_state_holds_builtin ( self, state ) :
        assert state.loglikelihood is Run.loglikelihood

    def test_custom_state_holds_custom ( self, custom_state ) :
        assert custom_state.loglikelihood is constant_loglikelihood

    def test_invalid_likelihood_fails_at_construction ( self ) :
        # validation must fire before any sampling machinery is touched
        with pytest.raises( TypeError, match = 'must be either None' ) :
            PipelineState( None, None, None, None, loglikelihood = 42 )

    def test_logprob_dispatches_to_custom ( self, custom_state ) :
        par = custom_state.handler.par_prior.mean( axis = 1 )
        assert Run.logprob( par, custom_state ) == pytest.approx( 42. )

    def test_logprob_prior_gate_still_applies ( self, custom_state ) :
        par_out = custom_state.handler.par_prior[:, 1] + 1.0
        assert Run.logprob( par_out, custom_state ) == -np.inf

    def test_state_without_attribute_falls_back_to_builtin ( self, state ) :
        # backward compatibility: states pickled before this feature have no
        # ``loglikelihood`` attribute and must keep working
        legacy = SimpleNamespace( data    = state.data,
                                  model   = state.model,
                                  noise   = state.noise,
                                  handler = state.handler )
        par = state.handler.par_prior.mean( axis = 1 )
        assert Run.logprob( par, legacy ) \
            == pytest.approx( Run.loglikelihood( par, state ) )


# ===========================================================================
# _expand_hyperpar
# ===========================================================================

def _make_hyperpar ( **overrides ) :
    """Minimal single-source SimpleNamespace mimicking a parameter file."""
    ns = SimpleNamespace(
        bands             = BANDS,
        fluxes            = FLUXES,
        errors            = ERRORS,
        uplims            = UPLIMS,
        filters           = BANDS,
        filters_custom    = None,
        galaxy_parameters = dict( BASE_PARAMS ),
        sfh_model         = 'insitu',
        ssp_lib           = 'parsec22.NT',
        do_AGN            = False,
        do_Radio          = False,
        do_Xray           = False,
        noise_model       = None,
        noise_parameters  = {},
        noise_kwargs      = {},
        lstep             = None,
        method_uplims     = 'chi2',
        sampler           = 'dynesty',
        nwalkers          = None,
        nsamples          = None,
        sampler_kw        = {},
        sampling_kw       = {},
        output_directory  = '/tmp',
        run_id            = '',
        store_method      = 'hdf5',
        store_lightweight = False,
        pickle_sampler    = False,
        pickle_raw        = False,
    )
    for k, v in overrides.items() :
        setattr( ns, k, v )
    return ns


class TestExpandHyperparLikelihood :

    def test_absent_attribute_gives_none ( self ) :
        # old parameter files without the ``loglikelihood`` variable
        jobs = _expand_hyperpar( _make_hyperpar() )
        assert jobs[0].loglikelihood is None

    def test_explicit_none_gives_none ( self ) :
        jobs = _expand_hyperpar( _make_hyperpar( loglikelihood = None ) )
        assert jobs[0].loglikelihood is None

    def test_custom_propagated_to_all_jobs ( self ) :
        N = 2
        fluxes = np.stack( [ FLUXES ] * N )
        errors = np.stack( [ ERRORS ] * N )
        uplims = np.zeros( ( N, len( BANDS ) ), dtype = bool )
        jobs = _expand_hyperpar( _make_hyperpar(
            fluxes = fluxes, errors = errors, uplims = uplims,
            loglikelihood = constant_loglikelihood,
        ) )
        assert len( jobs ) == N
        assert all( j.loglikelihood is constant_loglikelihood for j in jobs )

    def test_invalid_likelihood_fails_in_main_process ( self ) :
        with pytest.raises( TypeError, match = 'must be either None' ) :
            _expand_hyperpar( _make_hyperpar( loglikelihood = 42 ) )

    def test_per_variant_override_raises ( self ) :
        models = [ dict( do_AGN = False ),
                   dict( do_AGN = True,
                         loglikelihood = constant_loglikelihood ) ]
        with pytest.raises( ValueError, match = 'per-variant' ) :
            _expand_hyperpar( _make_hyperpar( models = models ) )


# ===========================================================================
# nautilus adapter
# ===========================================================================

class TestNautilusAdapter :

    def test_argument_order_restored ( self ) :
        # emulate the double partial-wrapping: galapy pre-binds (logl, state),
        # then nautilus wraps again with the likelihood keyword arguments
        logl    = RecordingLoglikelihood()
        wrapped = functools.partial( _nautilus_loglikelihood, logl, 'STATE' )
        nautilus_side = functools.partial( wrapped, method_uplims = 'chi2' )
        ret = nautilus_side( 'PAR' )
        assert ret == pytest.approx( 7. )
        ( par, st, kwargs ), = logl.calls
        assert par == 'PAR'
        assert st  == 'STATE'
        assert kwargs == { 'method_uplims' : 'chi2' }

    def test_matches_builtin_likelihood ( self, state ) :
        # the adapter must be a transparent reordering: same value as the
        # direct call used by the dynesty and emcee paths
        par     = state.handler.par_prior.mean( axis = 1 )
        wrapped = functools.partial( _nautilus_loglikelihood,
                                     Run.loglikelihood, state )
        assert wrapped( par ) == pytest.approx( Run.loglikelihood( par, state ) )

    def test_partial_is_picklable ( self ) :
        # nautilus' internal pool pickles the wrapped likelihood: the adapter
        # is module-level on purpose, so the partial serialises by reference
        wrapped  = functools.partial( _nautilus_loglikelihood,
                                      Run.loglikelihood, None )
        restored = pickle.loads( pickle.dumps( wrapped ) )
        assert restored.func is _nautilus_loglikelihood
        assert restored.args[0] is Run.loglikelihood


# ===========================================================================
# Parallel-path dispatch
# ===========================================================================

class _FakePool :
    """Stand-in for a multiprocessing pool: no worker process is started."""
    def __enter__ ( self ) : return self
    def __exit__ ( self, *args ) : return False

class _FakeContext :
    def __init__ ( self ) : self.requested = []
    def Pool ( self, Ncpu ) :
        self.requested.append( Ncpu )
        return _FakePool()


@pytest.fixture
def no_subprocess ( monkeypatch ) :
    """Intercept pool creation and result storage in _sample_parallel.

    The dispatch is what is under test, so neither a real pool nor a real
    sampling run is needed — spawning workers here would make the test slow
    and platform-dependent for no added coverage.
    """
    import multiprocessing

    ctx = _FakeContext()
    monkeypatch.setattr( multiprocessing, 'get_context', lambda _ : ctx )
    monkeypatch.setattr( Run, 'store_results', lambda *a, **kw : None )

    calls = []
    def _fake_sample ( state, **kwargs ) :
        calls.append( ( state, kwargs ) )
        return 'SAMPLER'
    monkeypatch.setattr( Run, 'sample', _fake_sample )

    return SimpleNamespace( ctx = ctx, calls = calls )


class TestParallelDispatch :

    def test_nautilus_is_accepted ( self, state, no_subprocess ) :
        # regression: nautilus was missing from the _sample_parallel dispatch
        # and every non-serial run aborted on the else branch
        Run._sample_parallel( state, which_sampler = 'nautilus', Ncpu = 2 )

        ( ( _, kwargs ), ) = no_subprocess.calls
        assert kwargs['sampler'] == 'nautilus'
        assert kwargs['Ncpu']    == 2
        assert isinstance( kwargs['pool'], _FakePool )
        assert no_subprocess.ctx.requested == [ 2 ]

    def test_custom_likelihood_travels_with_the_state ( self, custom_state,
                                                        no_subprocess ) :
        # sample() pre-binds the likelihood carried by the state, so the state
        # object itself must reach it unchanged on the parallel nautilus path
        Run._sample_parallel( custom_state, which_sampler = 'nautilus',
                              Ncpu = 1 )

        ( ( forwarded, _ ), ) = no_subprocess.calls
        assert forwarded is custom_state
        assert forwarded.loglikelihood is constant_loglikelihood

    def test_unknown_sampler_still_rejected ( self, state, no_subprocess ) :
        with pytest.raises( ValueError, match = 'is not valid' ) as excinfo :
            Run._sample_parallel( state, which_sampler = 'metropolis',
                                  Ncpu = 1 )
        # the message must advertise every supported sampler
        for name in ( 'dynesty', 'emcee', 'nautilus' ) :
            assert name in str( excinfo.value )


# ===========================================================================
# Results provenance marker
# ===========================================================================

class TestResultsProvenance :

    def test_default_is_builtin_marker ( self, results ) :
        assert results.loglikelihood_name == _default_loglikelihood_name

    def test_dump_stores_marker ( self, results ) :
        assert results.dump()[ 'loglikelihood_name' ] \
            == _default_loglikelihood_name

    def test_roundtrip_preserves_marker ( self, results ) :
        d = results.dump()
        d[ 'loglikelihood_name' ] = 'my_likelihoods.student_t'
        ret = Results.load( d )
        assert ret.loglikelihood_name == 'my_likelihoods.student_t'

    def test_file_without_marker_loads_as_builtin ( self, results ) :
        # files written before this release have no such key
        d = results.dump()
        d.pop( 'loglikelihood_name' )
        ret = Results.load( d )
        assert ret.loglikelihood_name == _default_loglikelihood_name


# ===========================================================================
# bayes_factor
# ===========================================================================

class TestBayesFactor :

    def _mock_run ( self, logz, name = _default_loglikelihood_name ) :
        return SimpleNamespace( logz = logz, loglikelihood_name = name )

    def test_raw_floats ( self ) :
        assert bayes_factor( 3., 1. ) == pytest.approx( 2. )

    def test_results_like_objects_same_builtin_no_warning ( self ) :
        with _warnings.catch_warnings( record = True ) as record :
            _warnings.simplefilter( 'always' )
            log_bf = bayes_factor( self._mock_run( 3. ), self._mock_run( 1. ) )
        assert log_bf == pytest.approx( 2. )
        assert len( record ) == 0

    def test_mixed_raw_and_object ( self ) :
        assert bayes_factor( self._mock_run( 3. ), 1. ) == pytest.approx( 2. )

    def test_different_likelihoods_warn ( self ) :
        with pytest.warns( UserWarning, match = 'different log-likelihoods' ) :
            bayes_factor( self._mock_run( 3., name = 'my_mod.like_a' ),
                          self._mock_run( 1., name = 'my_mod.like_b' ) )

    def test_same_custom_likelihood_warns ( self ) :
        with pytest.warns( UserWarning, match = 'custom log-likelihood' ) :
            bayes_factor( self._mock_run( 3., name = 'my_mod.like_a' ),
                          self._mock_run( 1., name = 'my_mod.like_a' ) )

    def test_missing_evidence_raises ( self ) :
        with pytest.raises( ValueError, match = 'evidence' ) :
            bayes_factor( self._mock_run( None ), self._mock_run( 1. ) )

    def test_results_without_evidence_raises ( self, results ) :
        # a real Results from a run with no evidence estimate (emcee-like)
        with pytest.raises( ValueError, match = 'evidence' ) :
            bayes_factor( results, results )
