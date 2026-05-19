#------------------------------------------------------------------------------#
#                  Test of galapy.sampling.Run module.                         #
#------------------------------------------------------------------------------#
"""
Tests for the galapy fitting pipeline (Run module).

Covers:
- Run.initialize: builds global_dict with model, observation, handler, noise
- Free vs. fixed parameter registration
- AGN template parameter keys (agn.template.*)
- SFH model variants
- Run.loglikelihood and Run.logprob behaviour
- Upper-limit handling
- Noise (calibration error) model integration
- Prior transform utility
"""

import numpy as np
import pytest

from galapy.sampling import Run
from galapy.Galaxy import PhotoGXY
from galapy.sampling.Observation import Observation
from galapy.sampling.Statistics import transform_to_prior_unit_cube
from galapy.Handlers import ModelParameters

# ---------------------------------------------------------------------------
# Shared test data
# ---------------------------------------------------------------------------

BANDS  = ['GOODS.b', 'GOODS.v', 'GOODS.i', 'GOODS.z']
FLUXES = np.array([1.0e-4, 3.0e-4, 9.0e-4, 2.0e-3])
ERRORS = FLUXES * 0.1
UPLIMS = np.zeros(len(BANDS), dtype=bool)

# Minimal galaxy parameter set – insitu SFH, a handful of free params
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

AGN_PARAMS = {
    'agn.fAGN'        : ([-3., 3.], True),
    'agn.template.ct' : 40,
    'agn.template.al' : 0.,
    'agn.template.be' : -0.5,
    'agn.template.ta' : 6.,
    'agn.template.rm' : 60,
    'agn.template.ia' : 0.001,
}

NOISE_PARAMS = {'f_cal': ([-10., 1.], True)}


def _init(**kwargs):
    """Call Run.initialize with shared defaults; return global_dict."""
    params = dict(BASE_PARAMS)
    params.update(kwargs.pop('extra_params', {}))
    Run.initialize(
        BANDS, FLUXES, ERRORS, UPLIMS,
        BANDS, params,
        **kwargs
    )
    return Run.global_dict


def _random_par():
    """Draw a random parameter vector uniformly from the prior."""
    h = Run.global_dict['handler']
    return h.rng.uniform(*h.par_prior.T)


# ===========================================================================
# Initialization
# ===========================================================================

class TestInitialize:

    def test_global_dict_has_required_keys(self):
        gd = _init()
        assert set(gd.keys()) >= {'data', 'model', 'noise', 'handler'}

    def test_model_is_photogxy(self):
        gd = _init()
        assert isinstance(gd['model'], PhotoGXY)

    def test_observation_stores_fluxes(self):
        gd = _init()
        assert isinstance(gd['data'], Observation)
        np.testing.assert_allclose(np.sort(gd['data'].fluxes), np.sort(FLUXES))

    def test_observation_stores_errors(self):
        gd = _init()
        np.testing.assert_allclose(np.sort(gd['data'].errors), np.sort(ERRORS))

    def test_observation_stores_uplims(self):
        gd = _init()
        assert not np.any(gd['data'].uplims)

    def test_noise_is_none_by_default(self):
        gd = _init()
        assert gd['noise'] is None

    def test_handler_is_model_parameters(self):
        gd = _init()
        assert isinstance(gd['handler'], ModelParameters)

    def test_with_noise_model(self):
        gd = _init(noise_model='calibration_error', noise_params=NOISE_PARAMS)
        assert gd['noise'] is not None

    def test_with_agn(self):
        gd = _init(do_AGN=True, extra_params=AGN_PARAMS)
        assert 'AGN' in gd['model'].components

    def test_with_radio(self):
        gd = _init(do_Radio=True)
        assert 'nebular' in gd['model'].components

    def test_with_xray(self):
        gd = _init(do_Xray=True)
        assert 'XRB' in gd['model'].components

    def test_invalid_noise_model_warns(self):
        with pytest.warns(UserWarning, match='not valid'):
            gd = _init(noise_model='nonexistent_model', noise_params={})
        assert gd['noise'] is None

    def test_custom_filter(self):
        wl    = np.linspace(3000, 10000, 100)
        trans = np.exp(-0.5 * ((wl - 5000) / 500) ** 2)
        params = dict(BASE_PARAMS)
        params['redshift'] = 0.5
        Run.initialize(
            ['my_band'], np.array([5e-4]), np.array([5e-5]), np.array([False]),
            [], params,
            filter_kwargs={'my_band': {'wavelengths': wl, 'photons': trans}},
        )
        assert Run.global_dict['data'] is not None


# ===========================================================================
# Parameter registration
# ===========================================================================

class TestParameterRegistration:

    def test_free_params_are_identified(self):
        _init()
        free = list(Run.global_dict['handler'].par_free)
        assert 'galaxy.age' in free
        assert 'galaxy.sfh.psi_max' in free

    def test_fixed_param_not_in_free(self):
        _init()
        free = list(Run.global_dict['handler'].par_free)
        assert 'galaxy.redshift' not in free

    def test_prior_shape_matches_free_count(self):
        _init()
        h = Run.global_dict['handler']
        assert h.par_prior.shape == (len(h.par_free), 2)

    def test_agn_template_keys_accepted_without_warning(self):
        import warnings
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter('always')
            _init(do_AGN=True, extra_params=AGN_PARAMS)
        agn_template_warnings = [
            w for w in record
            if 'agn.template' in str(w.message)
        ]
        assert len(agn_template_warnings) == 0

    def test_agn_fAGN_is_free(self):
        _init(do_AGN=True, extra_params=AGN_PARAMS)
        free = list(Run.global_dict['handler'].par_free)
        assert 'galaxy.agn.fAGN' in free

    def test_agn_old_keys_trigger_warning(self):
        old_key_params = dict(BASE_PARAMS)
        old_key_params.update({'agn.fAGN': ([-3., 3.], True), 'agn.ct': 40})
        with pytest.warns(UserWarning, match='galaxy.agn.ct'):
            Run.initialize(
                BANDS, FLUXES, ERRORS, UPLIMS,
                BANDS, old_key_params,
                do_AGN=True,
            )

    def test_noise_free_param_registered(self):
        _init(noise_model='calibration_error', noise_params=NOISE_PARAMS)
        free = list(Run.global_dict['handler'].par_free)
        assert 'noise.f_cal' in free


# ===========================================================================
# SFH model variants
# ===========================================================================

@pytest.mark.parametrize('sfh_model,extra', [
    ('insitu', {}),
    ('delayedexp', {
        'sfh.psi_norm' : ([0., 4.], True),
        'sfh.k_shape'  : 1.0,
        'sfh.Mdust'    : 1e8,
        'sfh.Zgxy'     : 0.02,
    }),
    ('lognormal', {
        'sfh.psi_norm'   : ([0., 4.], True),
        'sfh.sigma_star' : 1.0,
        'sfh.Mdust'      : 1e8,
        'sfh.Zgxy'       : 0.02,
    }),
    ('constant', {
        'sfh.psi'   : ([0., 4.], True),
        'sfh.Mdust' : 1e8,
        'sfh.Zgxy'  : 0.02,
    }),
])
def test_sfh_model_initializes_and_evaluates(sfh_model, extra):
    params = dict(BASE_PARAMS)
    if sfh_model != 'insitu':
        params.pop('sfh.psi_max', None)
        params.pop('sfh.tau_star', None)
    params.update(extra)
    Run.initialize(BANDS, FLUXES, ERRORS, UPLIMS, BANDS, params, sfh_model=sfh_model)
    gd = Run.global_dict
    par = gd['handler'].par_prior.mean(axis=1)
    ll  = Run.loglikelihood(par)
    assert isinstance(float(ll), float)


# ===========================================================================
# Loglikelihood
# ===========================================================================

class TestLoglikelihood:

    def setup_method(self):
        _init()

    def test_returns_float(self):
        ll = Run.loglikelihood(_random_par())
        assert isinstance(float(ll), float)

    def test_finite_or_neginf(self):
        ll = Run.loglikelihood(_random_par())
        assert np.isfinite(ll) or ll == -np.inf

    def test_finite_at_prior_midpoint(self):
        par = Run.global_dict['handler'].par_prior.mean(axis=1)
        assert np.isfinite(Run.loglikelihood(par))

    def test_all_uplims(self):
        uplims_all = np.ones(len(BANDS), dtype=bool)
        Run.initialize(BANDS, FLUXES, ERRORS, uplims_all, BANDS, BASE_PARAMS)
        ll = Run.loglikelihood(Run.global_dict['handler'].par_prior.mean(axis=1))
        assert isinstance(float(ll), float)

    def test_mixed_uplims(self):
        uplims_mixed = np.array([False, True, False, True])
        Run.initialize(BANDS, FLUXES, ERRORS, uplims_mixed, BANDS, BASE_PARAMS)
        ll = Run.loglikelihood(Run.global_dict['handler'].par_prior.mean(axis=1))
        assert isinstance(float(ll), float)

    @pytest.mark.parametrize('method', ['simple', 'chi2', 'S12'])
    def test_uplim_methods(self, method):
        uplims_mixed = np.array([False, True, False, False])
        Run.initialize(BANDS, FLUXES, ERRORS, uplims_mixed, BANDS, BASE_PARAMS)
        par = Run.global_dict['handler'].par_prior.mean(axis=1)
        ll  = Run.loglikelihood(par, method_uplims=method)
        assert isinstance(float(ll), float)

    def test_with_noise_model(self):
        _init(noise_model='calibration_error', noise_params=NOISE_PARAMS)
        par = Run.global_dict['handler'].par_prior.mean(axis=1)
        assert np.isfinite(Run.loglikelihood(par))

    def test_agn_model(self):
        _init(do_AGN=True, extra_params=AGN_PARAMS)
        par = Run.global_dict['handler'].par_prior.mean(axis=1)
        ll  = Run.loglikelihood(par)
        assert isinstance(float(ll), float)


# ===========================================================================
# Logprob (prior-bounded loglikelihood)
# ===========================================================================

class TestLogprob:

    def setup_method(self):
        _init()

    def test_finite_or_neginf_inside_prior(self):
        lp = Run.logprob(_random_par())
        assert np.isfinite(lp) or lp == -np.inf

    def test_above_prior_is_neginf(self):
        par_out = Run.global_dict['handler'].par_prior[:, 1] + 1.0
        assert Run.logprob(par_out) == -np.inf

    def test_below_prior_is_neginf(self):
        par_out = Run.global_dict['handler'].par_prior[:, 0] - 1.0
        assert Run.logprob(par_out) == -np.inf

    def test_midpoint_equals_loglikelihood(self):
        par = Run.global_dict['handler'].par_prior.mean(axis=1)
        assert Run.logprob(par) == pytest.approx(Run.loglikelihood(par))


# ===========================================================================
# Prior transform utility
# ===========================================================================

class TestPriorTransform:

    def test_zero_maps_to_lower_bound(self):
        prior  = np.array([[2., 5.], [0., 1.]])
        result = transform_to_prior_unit_cube(np.zeros(2), prior)
        np.testing.assert_allclose(result, prior[:, 0])

    def test_one_maps_to_upper_bound(self):
        prior  = np.array([[2., 5.], [0., 1.]])
        result = transform_to_prior_unit_cube(np.ones(2), prior)
        np.testing.assert_allclose(result, prior[:, 1])

    def test_half_maps_to_midpoint(self):
        prior  = np.array([[0., 10.], [4., 8.]])
        result = transform_to_prior_unit_cube(np.full(2, 0.5), prior)
        np.testing.assert_allclose(result, prior.mean(axis=1))

    def test_output_within_prior_range(self):
        prior   = np.array([[1., 4.], [-2., 3.]])
        samples = np.random.default_rng(0).random((100, 2))
        for s in samples:
            r = transform_to_prior_unit_cube(s, prior)
            assert np.all(r >= prior[:, 0])
            assert np.all(r <= prior[:, 1])
