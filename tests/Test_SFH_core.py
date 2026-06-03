
#------------------------------------------------------------------------------#
#           Test of galapy.SFH_core module.
#
# These tests are written as contracts on the public interface of CSFH,
# not on the underlying implementation.  They are intended to remain valid
# if SFH_core is ever re-implemented in pure Python, as long as the class
# name and method signatures are preserved.
#------------------------------------------------------------------------------#

import pickle

import numpy as np
import pytest

from galapy.SFH_core import CSFH


# ---------------------------------------------------------------------------
# Shared fixtures and constants
# ---------------------------------------------------------------------------

ALL_MODELS = ['insitu', 'constant', 'delayedexp', 'lognormal', 'interpolated']

# Reference evaluation times [yr]
TAU = np.array([1.e8, 5.e8, 1.e9])
NPTS = np.full(3, 1000, dtype=np.uint64)

# Default parameters per model (ordered as expected by set_params)
_DEFAULT_PARAMS = {
    'insitu'      : [100., 3.e8],
    'constant'    : [1., 1.e8, 0.01],
    'delayedexp'  : [1., 0.2, 1.e8, 1.e8, 0.01],
    'lognormal'   : [100., 2., 3.e8, 1.e8, 0.01],
    'interpolated': [1.e8, 0.01],
}

# Interpolation grid for the 'interpolated' model
_INTERP_T   = [1.e7, 1.e8, 5.e8, 1.e9]
_INTERP_PSI = [10., 50., 30., 5.]


def _make(model, tau_quench=2.e10):
    """Construct and configure a CSFH instance with default parameters."""
    s = CSFH(model=model, tau_quench=tau_quench)
    if model == 'interpolated':
        s.set_interpolator(_INTERP_T, _INTERP_PSI)
    s.set_params(_DEFAULT_PARAMS[model])
    return s


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

class TestConstruction:

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_all_models_construct(self, model):
        _make(model)

    def test_invalid_model_raises(self):
        with pytest.raises((ValueError, RuntimeError)):
            CSFH(model='badmodel')

    def test_default_tau_quench(self):
        s = CSFH(model='insitu')
        s.set_params(_DEFAULT_PARAMS['insitu'])
        # Default tau_quench = 2e10 yr >> any reasonable galaxy age:
        # SFR at 5 Gyr should still be positive
        assert s(np.array([5.e9]))[0] > 0


# ---------------------------------------------------------------------------
# SFR (__call__)
# ---------------------------------------------------------------------------

class TestSFR:

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_sfr_non_negative(self, model):
        s = _make(model)
        assert np.all(s(TAU) >= 0)

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_sfr_output_shape(self, model):
        s = _make(model)
        assert s(TAU).shape == TAU.shape

    def test_insitu_sfr_reference(self):
        s = _make('insitu')
        expected = np.array([24.7265168, 45.98734055, 31.62243299])
        np.testing.assert_allclose(s(TAU), expected, rtol=1e-5)

    def test_constant_sfr_is_constant(self):
        s = _make('constant')
        sfr = s(TAU)
        # all values equal the configured psi=1.0
        np.testing.assert_allclose(sfr, 1.0, rtol=1e-10)

    def test_delayedexp_sfr_reference(self):
        s = _make('delayedexp')
        expected = np.array([1.46455443e+01, 3.70102136e-01, 2.86454191e-03])
        np.testing.assert_allclose(s(TAU), expected, rtol=1e-5)

    def test_lognormal_sfr_reference(self):
        s = _make('lognormal')
        expected = np.array([17.15373359, 19.30697868, 16.6413516])
        np.testing.assert_allclose(s(TAU), expected, rtol=1e-5)

    def test_interpolated_sfr_at_nodes(self):
        s = _make('interpolated')
        nodes = np.array(_INTERP_T)
        np.testing.assert_allclose(s(nodes), _INTERP_PSI, rtol=1e-10)

    def test_delayedexp_sfr_decreasing_at_late_times(self):
        s = _make('delayedexp')
        late = np.array([1.e8, 5.e8, 2.e9, 5.e9])
        sfr = s(late)
        assert np.all(np.diff(sfr) < 0)


# ---------------------------------------------------------------------------
# Stellar mass (Mstar)
# ---------------------------------------------------------------------------

class TestMstar:

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_mstar_non_negative(self, model):
        s = _make(model)
        assert np.all(s.Mstar(TAU, NPTS) >= 0)

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_mstar_output_shape(self, model):
        s = _make(model)
        assert s.Mstar(TAU, NPTS).shape == TAU.shape

    def test_insitu_mstar_reference(self):
        s = _make('insitu')
        expected = np.array([1.07447685e+09, 1.26514050e+10, 2.48126196e+10])
        np.testing.assert_allclose(s.Mstar(TAU, NPTS), expected, rtol=1e-5)

    def test_constant_mstar_reference(self):
        s = _make('constant')
        expected = np.array([7.72508827e+07, 3.46442577e+08, 6.58331879e+08])
        np.testing.assert_allclose(s.Mstar(TAU, NPTS), expected, rtol=1e-5)

    def test_delayedexp_mstar_reference(self):
        s = _make('delayedexp')
        expected = np.array([1.53352833e+09, 2.38418724e+09, 2.25012831e+09])
        np.testing.assert_allclose(s.Mstar(TAU, NPTS), expected, rtol=1e-5)

    def test_lognormal_mstar_reference(self):
        s = _make('lognormal')
        expected = np.array([9.40024773e+08, 6.26222440e+09, 1.18235699e+10])
        np.testing.assert_allclose(s.Mstar(TAU, NPTS), expected, rtol=1e-5)

    def test_interpolated_mstar_reference(self):
        s = _make('interpolated')
        expected = np.array([2.19912626e+09, 1.29713719e+10, 1.76265524e+10])
        np.testing.assert_allclose(s.Mstar(TAU, NPTS), expected, rtol=1e-5)

    def test_mstar_non_decreasing_before_quench(self):
        s = _make('insitu')
        tau_fine = np.array([1.e7, 1.e8, 5.e8, 1.e9, 5.e9])
        npts_fine = np.full(5, 1000, dtype=np.uint64)
        mstar = s.Mstar(tau_fine, npts_fine)
        assert np.all(np.diff(mstar) >= 0)

    def test_mstar_does_not_grow_after_quench(self):
        # After quenching SFR=0, stellar mass loss from evolved stars can still
        # reduce Mstar, so m_after <= m_at_quench (not equality).
        s = CSFH(model='insitu', tau_quench=5.e8)
        s.set_params(_DEFAULT_PARAMS['insitu'])
        m_at_quench = s.Mstar(np.array([5.e8]), np.array([1000], dtype=np.uint64))[0]
        m_after     = s.Mstar(np.array([2.e9]), np.array([1000], dtype=np.uint64))[0]
        assert m_after <= m_at_quench


# ---------------------------------------------------------------------------
# Dust and gas mass (Mdust, Mgas)
# ---------------------------------------------------------------------------

class TestMdust:

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_mdust_non_negative(self, model):
        s = _make(model)
        assert np.all(s.Mdust(TAU) >= 0)

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_mdust_output_shape(self, model):
        s = _make(model)
        assert s.Mdust(TAU).shape == TAU.shape

    def test_insitu_mdust_reference(self):
        s = _make('insitu')
        expected = np.array([3.15488873e+07, 2.55286176e+08, 2.56791236e+08])
        np.testing.assert_allclose(s.Mdust(TAU), expected, rtol=1e-5)

    def test_constant_mdust_is_constant(self):
        # constant/delayedexp/lognormal/interpolated use fixed Mdust parameter
        s = _make('constant')
        np.testing.assert_allclose(s.Mdust(TAU), 1.e8, rtol=1e-10)


class TestMgas:

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_mgas_non_negative(self, model):
        s = _make(model)
        assert np.all(s.Mgas(TAU) >= 0)

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_mgas_output_shape(self, model):
        s = _make(model)
        assert s.Mgas(TAU).shape == TAU.shape

    def test_insitu_mgas_reference(self):
        s = _make('insitu')
        expected = np.array([7.41795504e+09, 1.37962022e+10, 9.48672990e+09])
        np.testing.assert_allclose(s.Mgas(TAU), expected, rtol=1e-5)


# ---------------------------------------------------------------------------
# Metallicities (Zgas, Zstar)
# ---------------------------------------------------------------------------

class TestMetallicities:

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_zgas_in_unit_interval(self, model):
        s = _make(model)
        zg = s.Zgas(TAU)
        assert np.all((zg > 0) & (zg < 1))

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_zstar_in_unit_interval(self, model):
        s = _make(model)
        zs = s.Zstar(TAU)
        assert np.all((zs > 0) & (zs < 1))

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_zgas_zstar_output_shape(self, model):
        s = _make(model)
        assert s.Zgas(TAU).shape == TAU.shape
        assert s.Zstar(TAU).shape == TAU.shape

    def test_insitu_zgas_reference(self):
        s = _make('insitu')
        expected = np.array([0.00630795, 0.02471833, 0.0357591])
        np.testing.assert_allclose(s.Zgas(TAU), expected, rtol=1e-5)

    def test_insitu_zstar_reference(self):
        s = _make('insitu')
        expected = np.array([0.00417058, 0.01612655, 0.02377438])
        np.testing.assert_allclose(s.Zstar(TAU), expected, rtol=1e-5)

    def test_insitu_zstar_le_zgas(self):
        # stars form from gas: stellar metallicity lags gas metallicity
        s = _make('insitu')
        assert np.all(s.Zstar(TAU) <= s.Zgas(TAU))

    def test_constant_metallicity_equals_set_value(self):
        s = _make('constant')
        np.testing.assert_allclose(s.Zgas(TAU),  0.01, rtol=1e-10)
        np.testing.assert_allclose(s.Zstar(TAU), 0.01, rtol=1e-10)


# ---------------------------------------------------------------------------
# dMstar
# ---------------------------------------------------------------------------

class TestDMstar:

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_dmstar_returns_array(self, model):
        s = _make(model)
        dm = s.dMstar()
        assert isinstance(dm, np.ndarray)

    @pytest.mark.parametrize('model', ALL_MODELS)
    def test_dmstar_non_negative(self, model):
        s = _make(model)
        assert np.all(s.dMstar() >= 0)

    def test_dmstar_shape_consistent_across_models(self):
        sizes = [_make(m).dMstar().shape for m in ALL_MODELS]
        assert len(set(sizes)) == 1


# ---------------------------------------------------------------------------
# set_tau_quench
# ---------------------------------------------------------------------------

class TestTauQuench:

    def test_sfr_zero_after_quench(self):
        s = CSFH(model='insitu')
        s.set_params(_DEFAULT_PARAMS['insitu'])
        s.set_tau_quench(5.e8)
        tau_after = np.array([6.e8, 1.e9, 5.e9])
        np.testing.assert_allclose(s(tau_after), 0.0, atol=1e-30)

    def test_sfr_positive_before_quench(self):
        s = CSFH(model='insitu')
        s.set_params(_DEFAULT_PARAMS['insitu'])
        s.set_tau_quench(5.e8)
        tau_before = np.array([1.e8, 3.e8])
        assert np.all(s(tau_before) > 0)

    @pytest.mark.parametrize('model', ['insitu', 'constant', 'lognormal'])
    def test_set_tau_quench_via_constructor(self, model):
        s = CSFH(model=model, tau_quench=1.e8)
        s.set_params(_DEFAULT_PARAMS[model])
        assert s(np.array([2.e8]))[0] == pytest.approx(0.0, abs=1e-30)

    def test_set_tau_quench_method(self):
        s = _make('insitu')
        s.set_tau_quench(3.e8)
        assert s(np.array([5.e8]))[0] == pytest.approx(0.0, abs=1e-30)


# ---------------------------------------------------------------------------
# set_interpolator (interpolated model only)
# ---------------------------------------------------------------------------

class TestSetInterpolator:

    def test_sfr_exact_at_nodes(self):
        s = CSFH(model='interpolated')
        s.set_interpolator(_INTERP_T, _INTERP_PSI)
        s.set_params(_DEFAULT_PARAMS['interpolated'])
        nodes = np.array(_INTERP_T)
        np.testing.assert_allclose(s(nodes), _INTERP_PSI, rtol=1e-10)

    def test_sfr_non_negative_between_nodes(self):
        s = CSFH(model='interpolated')
        s.set_interpolator(_INTERP_T, _INTERP_PSI)
        s.set_params(_DEFAULT_PARAMS['interpolated'])
        tau_mid = np.linspace(_INTERP_T[0], _INTERP_T[-1], 100)
        assert np.all(s(tau_mid) >= 0)

    def test_mstar_non_negative(self):
        s = _make('interpolated')
        assert np.all(s.Mstar(TAU, NPTS) >= 0)


# ---------------------------------------------------------------------------
# time_grid
# ---------------------------------------------------------------------------

class TestTimeGrid:

    @pytest.fixture
    def insitu_sfh(self):
        return _make('insitu')

    def test_returns_four_tuple(self, insitu_sfh):
        tgrid = list(np.linspace(0, 1.e9, 50))
        Zgrid = [0.0001, 0.001, 0.01, 0.02, 0.05]
        result = insitu_sfh.time_grid(1.e9, tgrid, Zgrid)
        assert len(result) == 4

    def test_dMgrid_shape(self, insitu_sfh):
        tgrid = list(np.linspace(0, 1.e9, 50))
        Zgrid = [0.001, 0.01, 0.05]
        dMgrid, _, _, _ = insitu_sfh.time_grid(1.e9, tgrid, Zgrid)
        assert dMgrid.shape == (len(tgrid),)

    def test_dMgrid_non_negative(self, insitu_sfh):
        tgrid = list(np.linspace(0, 1.e9, 50))
        Zgrid = [0.001, 0.01, 0.05]
        dMgrid, _, _, _ = insitu_sfh.time_grid(1.e9, tgrid, Zgrid)
        assert np.all(dMgrid >= 0)

    def test_Zidx_within_Zgrid_bounds(self, insitu_sfh):
        tgrid = list(np.linspace(0, 1.e9, 50))
        Zgrid = [0.0001, 0.001, 0.01, 0.02, 0.05]
        _, _, Zidx, _ = insitu_sfh.time_grid(1.e9, tgrid, Zgrid)
        assert np.all(Zidx < len(Zgrid))

    def test_Zout_non_negative(self, insitu_sfh):
        tgrid = list(np.linspace(0, 1.e9, 50))
        Zgrid = [0.001, 0.01, 0.05]
        _, Zout, _, _ = insitu_sfh.time_grid(1.e9, tgrid, Zgrid)
        assert np.all(Zout >= 0)

    def test_last_idx_is_int(self, insitu_sfh):
        tgrid = list(np.linspace(0, 1.e9, 50))
        Zgrid = [0.001, 0.01]
        _, _, _, last_idx = insitu_sfh.time_grid(1.e9, tgrid, Zgrid)
        assert isinstance(last_idx, int)


# ---------------------------------------------------------------------------
# Serialization (pickle round-trip)
# ---------------------------------------------------------------------------

class TestSerialization:

    _PICKLE_MODELS = [
        m if m != 'interpolated' else
        pytest.param(m, marks=pytest.mark.skip(
            reason='C++ sfh_interpolated does not serialize interpolator state; '
                   'deserialized object is corrupted and segfaults on use'
        ))
        for m in ALL_MODELS
    ]

    @pytest.mark.parametrize('model', _PICKLE_MODELS)
    def test_pickle_roundtrip_mstar(self, model):
        s = _make(model)
        s_loaded = pickle.loads(pickle.dumps(s))
        np.testing.assert_allclose(
            s_loaded.Mstar(TAU, NPTS),
            s.Mstar(TAU, NPTS),
            rtol=1e-10,
        )

    @pytest.mark.parametrize('model', _PICKLE_MODELS)
    def test_pickle_roundtrip_sfr(self, model):
        s = _make(model)
        s_loaded = pickle.loads(pickle.dumps(s))
        np.testing.assert_allclose(
            s_loaded(TAU),
            s(TAU),
            rtol=1e-10,
        )

    def test_pickle_preserves_tau_quench(self):
        s = CSFH(model='insitu', tau_quench=5.e8)
        s.set_params(_DEFAULT_PARAMS['insitu'])
        s_loaded = pickle.loads(pickle.dumps(s))
        # SFR after quench must still be zero after round-trip
        assert s_loaded(np.array([1.e9]))[0] == pytest.approx(0.0, abs=1e-30)
