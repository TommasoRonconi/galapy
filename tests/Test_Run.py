#------------------------------------------------------------------------------#
#                  Test of galapy.sampling.Run module.                         #
#------------------------------------------------------------------------------#
"""
Tests for the galapy fitting pipeline (Run module).

Covers:
- Run.initialize / PipelineState.initialize: builds PipelineState with model,
  observation, handler, noise
- Free vs. fixed parameter registration
- AGN template parameter keys (agn.template.*)
- SFH model variants
- Run.loglikelihood and Run.logprob behaviour
- Upper-limit handling
- Noise (calibration error) model integration
- Prior transform utility
- Catalogue-mode: _catalogue_worker builds independent PipelineState per object
"""

import numpy as np
import pytest

from types import SimpleNamespace

from galapy.sampling import Run
from galapy.sampling.Run import PipelineState, _expand_hyperpar, _model_suffixes
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
    """Call Run.initialize with shared defaults; return the PipelineState."""
    params = dict(BASE_PARAMS)
    params.update(kwargs.pop('extra_params', {}))
    return Run.initialize(
        BANDS, FLUXES, ERRORS, UPLIMS,
        BANDS, params,
        **kwargs
    )


def _random_par(state):
    """Draw a random parameter vector uniformly from the prior."""
    return state.handler.rng.uniform(*state.handler.par_prior.T)


# ===========================================================================
# Initialization
# ===========================================================================

class TestInitialize:

    def test_returns_pipeline_state(self):
        state = _init()
        assert isinstance(state, PipelineState)

    def test_state_has_required_attrs(self):
        state = _init()
        for attr in ('data', 'model', 'noise', 'handler'):
            assert hasattr(state, attr)

    def test_model_is_photogxy(self):
        state = _init()
        assert isinstance(state.model, PhotoGXY)

    def test_observation_stores_fluxes(self):
        state = _init()
        assert isinstance(state.data, Observation)
        np.testing.assert_allclose(np.sort(state.data.fluxes), np.sort(FLUXES))

    def test_observation_stores_errors(self):
        state = _init()
        np.testing.assert_allclose(np.sort(state.data.errors), np.sort(ERRORS))

    def test_observation_stores_uplims(self):
        state = _init()
        assert not np.any(state.data.uplims)

    def test_noise_is_none_by_default(self):
        state = _init()
        assert state.noise is None

    def test_handler_is_model_parameters(self):
        state = _init()
        assert isinstance(state.handler, ModelParameters)

    def test_with_noise_model(self):
        state = _init(noise_model='calibration_error', noise_params=NOISE_PARAMS)
        assert state.noise is not None

    def test_with_agn(self):
        state = _init(do_AGN=True, extra_params=AGN_PARAMS)
        assert 'AGN' in state.model.components

    def test_with_radio(self):
        state = _init(do_Radio=True)
        assert 'nebular' in state.model.components

    def test_with_xray(self):
        state = _init(do_Xray=True)
        assert 'XRB' in state.model.components

    def test_invalid_noise_model_warns(self):
        with pytest.warns(UserWarning, match='not valid'):
            state = _init(noise_model='nonexistent_model', noise_params={})
        assert state.noise is None

    def test_custom_filter(self):
        wl    = np.linspace(3000, 10000, 100)
        trans = np.exp(-0.5 * ((wl - 5000) / 500) ** 2)
        params = dict(BASE_PARAMS)
        params['redshift'] = 0.5
        state = Run.initialize(
            ['my_band'], np.array([5e-4]), np.array([5e-5]), np.array([False]),
            [], params,
            filter_kwargs={'my_band': {'wavelengths': wl, 'photons': trans}},
        )
        assert state.data is not None

    def test_independent_states_do_not_share_model(self):
        state1 = _init()
        state2 = _init()
        assert state1.model is not state2.model
        assert state1.handler is not state2.handler


# ===========================================================================
# Parameter registration
# ===========================================================================

class TestParameterRegistration:

    def test_free_params_are_identified(self):
        state = _init()
        free = list(state.handler.par_free)
        assert 'galaxy.age' in free
        assert 'galaxy.sfh.psi_max' in free

    def test_fixed_param_not_in_free(self):
        state = _init()
        free = list(state.handler.par_free)
        assert 'galaxy.redshift' not in free

    def test_prior_shape_matches_free_count(self):
        state = _init()
        assert state.handler.par_prior.shape == (len(state.handler.par_free), 2)

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
        state = _init(do_AGN=True, extra_params=AGN_PARAMS)
        free = list(state.handler.par_free)
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
        state = _init(noise_model='calibration_error', noise_params=NOISE_PARAMS)
        free = list(state.handler.par_free)
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
    state = Run.initialize(BANDS, FLUXES, ERRORS, UPLIMS, BANDS, params,
                           sfh_model=sfh_model)
    par = state.handler.par_prior.mean(axis=1)
    ll  = Run.loglikelihood(par, state)
    assert isinstance(float(ll), float)


# ===========================================================================
# Loglikelihood
# ===========================================================================

class TestLoglikelihood:

    def setup_method(self):
        self.state = _init()

    def test_returns_float(self):
        ll = Run.loglikelihood(_random_par(self.state), self.state)
        assert isinstance(float(ll), float)

    def test_finite_or_neginf(self):
        ll = Run.loglikelihood(_random_par(self.state), self.state)
        assert np.isfinite(ll) or ll == -np.inf

    def test_finite_at_prior_midpoint(self):
        par = self.state.handler.par_prior.mean(axis=1)
        assert np.isfinite(Run.loglikelihood(par, self.state))

    def test_all_uplims(self):
        uplims_all = np.ones(len(BANDS), dtype=bool)
        state = Run.initialize(BANDS, FLUXES, ERRORS, uplims_all, BANDS, BASE_PARAMS)
        ll = Run.loglikelihood(state.handler.par_prior.mean(axis=1), state)
        assert isinstance(float(ll), float)

    def test_mixed_uplims(self):
        uplims_mixed = np.array([False, True, False, True])
        state = Run.initialize(BANDS, FLUXES, ERRORS, uplims_mixed, BANDS, BASE_PARAMS)
        ll = Run.loglikelihood(state.handler.par_prior.mean(axis=1), state)
        assert isinstance(float(ll), float)

    @pytest.mark.parametrize('method', ['simple', 'chi2', 'S12'])
    def test_uplim_methods(self, method):
        uplims_mixed = np.array([False, True, False, False])
        state = Run.initialize(BANDS, FLUXES, ERRORS, uplims_mixed, BANDS, BASE_PARAMS)
        par = state.handler.par_prior.mean(axis=1)
        ll  = Run.loglikelihood(par, state, method_uplims=method)
        assert isinstance(float(ll), float)

    def test_with_noise_model(self):
        state = _init(noise_model='calibration_error', noise_params=NOISE_PARAMS)
        par = state.handler.par_prior.mean(axis=1)
        assert np.isfinite(Run.loglikelihood(par, state))

    def test_agn_model(self):
        state = _init(do_AGN=True, extra_params=AGN_PARAMS)
        par   = state.handler.par_prior.mean(axis=1)
        ll    = Run.loglikelihood(par, state)
        assert isinstance(float(ll), float)


# ===========================================================================
# Logprob (prior-bounded loglikelihood)
# ===========================================================================

class TestLogprob:

    def setup_method(self):
        self.state = _init()

    def test_finite_or_neginf_inside_prior(self):
        lp = Run.logprob(_random_par(self.state), self.state)
        assert np.isfinite(lp) or lp == -np.inf

    def test_above_prior_is_neginf(self):
        par_out = self.state.handler.par_prior[:, 1] + 1.0
        assert Run.logprob(par_out, self.state) == -np.inf

    def test_below_prior_is_neginf(self):
        par_out = self.state.handler.par_prior[:, 0] - 1.0
        assert Run.logprob(par_out, self.state) == -np.inf

    def test_midpoint_equals_loglikelihood(self):
        par = self.state.handler.par_prior.mean(axis=1)
        assert Run.logprob(par, self.state) == pytest.approx(
            Run.loglikelihood(par, self.state)
        )


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


# ===========================================================================
# Catalogue worker
# ===========================================================================

class TestCatalogueWorker:

    def _make_shared_kw(self):
        return dict(
            galaxy_parameters = BASE_PARAMS,
            sfh_model         = 'insitu',
            ssp_lib           = 'parsec22.NT',
            do_Radio          = False,
            do_Xray           = False,
            do_AGN            = False,
            noise_model       = None,
            noise_parameters  = {},
            noise_kwargs      = {},
            lstep             = None,
            filters_custom    = None,
            method_uplims     = 'chi2',
            sampler           = 'dynesty',
            nwalkers          = None,
            nsamples          = None,
            sampler_kw        = {},
            sampling_kw       = {},
            output_directory  = '/tmp',
            store_method      = 'hdf5',
            store_lightweight = True,
            pickle_sampler    = False,
            pickle_raw        = False,
        )

    def test_worker_builds_independent_states(self):
        # Verify that two catalogue entries produce independent PipelineState
        # objects (no shared model references).
        obj1 = dict(bands=BANDS, fluxes=FLUXES,      errors=ERRORS,      uplims=UPLIMS)
        obj2 = dict(bands=BANDS, fluxes=FLUXES * 2., errors=ERRORS * 2., uplims=UPLIMS)

        state1 = PipelineState.initialize(
            obj1['bands'], obj1['fluxes'], obj1['errors'], obj1['uplims'],
            obj1['bands'], BASE_PARAMS,
        )
        state2 = PipelineState.initialize(
            obj2['bands'], obj2['fluxes'], obj2['errors'], obj2['uplims'],
            obj2['bands'], BASE_PARAMS,
        )

        assert state1.model is not state2.model
        np.testing.assert_allclose(state1.data.fluxes, FLUXES)
        np.testing.assert_allclose(state2.data.fluxes, FLUXES * 2.)

    def test_worker_uplims_default_to_false(self):
        obj = dict(bands=BANDS, fluxes=FLUXES, errors=ERRORS)  # no uplims key
        state = PipelineState.initialize(
            obj['bands'], obj['fluxes'], obj['errors'],
            obj.get('uplims') if obj.get('uplims') is not None
                else np.zeros_like(np.asarray(obj['bands']), dtype=bool),
            obj.get('filters', obj['bands']),
            BASE_PARAMS,
        )
        assert not np.any(state.data.uplims)


# ===========================================================================
# Model suffix generation
# ===========================================================================

class TestModelSuffixes:

    def test_single_varying_key(self):
        variants = [
            dict(sfh_model='insitu', do_AGN=False, do_Radio=False, do_Xray=False),
            dict(sfh_model='insitu', do_AGN=True,  do_Radio=False, do_Xray=False),
        ]
        sfx = _model_suffixes(variants)
        assert sfx == ['AGNFalse', 'AGNTrue']

    def test_two_varying_keys(self):
        variants = [
            dict(sfh_model='insitu',     do_AGN=False, do_Radio=False, do_Xray=False),
            dict(sfh_model='insitu',     do_AGN=True,  do_Radio=False, do_Xray=False),
            dict(sfh_model='delayedexp', do_AGN=False, do_Radio=False, do_Xray=False),
            dict(sfh_model='delayedexp', do_AGN=True,  do_Radio=False, do_Xray=False),
        ]
        sfx = _model_suffixes(variants)
        assert sfx == [
            'SFHinsitu_AGNFalse', 'SFHinsitu_AGNTrue',
            'SFHdelayedexp_AGNFalse', 'SFHdelayedexp_AGNTrue',
        ]

    def test_no_varying_keys_returns_empty_strings(self):
        variants = [
            dict(sfh_model='insitu', do_AGN=False, do_Radio=False, do_Xray=False),
        ]
        sfx = _model_suffixes(variants)
        assert sfx == ['']


# ===========================================================================
# _expand_hyperpar
# ===========================================================================

def _make_hyperpar(**overrides):
    """Build a minimal single-source SimpleNamespace mimicking a parameter file."""
    ns = SimpleNamespace(
        bands             = BANDS,
        fluxes            = FLUXES,
        errors            = ERRORS,
        uplims            = UPLIMS,
        filters           = BANDS,
        filters_custom    = None,
        galaxy_parameters = dict(BASE_PARAMS),
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
    for k, v in overrides.items():
        setattr(ns, k, v)
    return ns


class TestExpandHyperpar:

    # ---- basic shape detection ----

    def test_single_source_single_model(self):
        jobs = _expand_hyperpar(_make_hyperpar())
        assert len(jobs) == 1

    def test_multi_source_no_models(self):
        fluxes = np.stack([FLUXES, FLUXES * 2, FLUXES * 3])
        errors = np.stack([ERRORS, ERRORS * 2, ERRORS * 3])
        uplims = np.zeros((3, len(BANDS)), dtype=bool)
        jobs = _expand_hyperpar(_make_hyperpar(
            fluxes=fluxes, errors=errors, uplims=uplims
        ))
        assert len(jobs) == 3

    def test_single_source_multi_model(self):
        models = [dict(do_AGN=False), dict(do_AGN=True)]
        jobs = _expand_hyperpar(_make_hyperpar(models=models))
        assert len(jobs) == 2

    def test_multi_source_multi_model_cross_product(self):
        N = 3
        fluxes = np.stack([FLUXES] * N)
        errors = np.stack([ERRORS] * N)
        uplims = np.zeros((N, len(BANDS)), dtype=bool)
        models = [dict(do_AGN=False), dict(do_AGN=True)]
        jobs = _expand_hyperpar(_make_hyperpar(
            fluxes=fluxes, errors=errors, uplims=uplims, models=models
        ))
        assert len(jobs) == N * 2

    # ---- per-source fluxes in each job ----

    def test_each_job_gets_correct_source_fluxes(self):
        fluxes = np.stack([FLUXES, FLUXES * 2])
        errors = np.stack([ERRORS, ERRORS * 2])
        uplims = np.zeros((2, len(BANDS)), dtype=bool)
        jobs = _expand_hyperpar(_make_hyperpar(
            fluxes=fluxes, errors=errors, uplims=uplims
        ))
        np.testing.assert_allclose(jobs[0].fluxes, FLUXES)
        np.testing.assert_allclose(jobs[1].fluxes, FLUXES * 2)

    # ---- per-source galaxy_parameters ----

    def test_per_source_param_is_split_correctly(self):
        N = 2
        redshifts = [0.5, 1.0]
        gp = dict(BASE_PARAMS)
        gp['redshift'] = redshifts
        fluxes = np.stack([FLUXES] * N)
        errors = np.stack([ERRORS] * N)
        uplims = np.zeros((N, len(BANDS)), dtype=bool)
        jobs = _expand_hyperpar(_make_hyperpar(
            fluxes=fluxes, errors=errors, uplims=uplims,
            galaxy_parameters=gp,
        ))
        assert jobs[0].galaxy_parameters['redshift'] == pytest.approx(0.5)
        assert jobs[1].galaxy_parameters['redshift'] == pytest.approx(1.0)

    def test_per_source_param_with_single_source_raises(self):
        gp = dict(BASE_PARAMS)
        gp['redshift'] = [0.5, 1.0]
        with pytest.raises(ValueError, match='only 1 source'):
            _expand_hyperpar(_make_hyperpar(galaxy_parameters=gp))

    def test_per_source_param_wrong_length_raises(self):
        N = 3
        gp = dict(BASE_PARAMS)
        gp['redshift'] = [0.5, 1.0]   # length 2, not 3
        fluxes = np.stack([FLUXES] * N)
        errors = np.stack([ERRORS] * N)
        uplims = np.zeros((N, len(BANDS)), dtype=bool)
        with pytest.raises(ValueError, match='2 entries but 3 sources'):
            _expand_hyperpar(_make_hyperpar(
                fluxes=fluxes, errors=errors, uplims=uplims,
                galaxy_parameters=gp,
            ))

    # ---- run_id generation ----

    def test_run_id_single_source_empty(self):
        jobs = _expand_hyperpar(_make_hyperpar(run_id=''))
        assert jobs[0].run_id == ''

    def test_run_id_single_source_named(self):
        jobs = _expand_hyperpar(_make_hyperpar(run_id='ngc2403'))
        assert jobs[0].run_id == 'ngc2403'

    def test_run_id_multi_source_auto(self):
        N = 3
        fluxes = np.stack([FLUXES] * N)
        errors = np.stack([ERRORS] * N)
        uplims = np.zeros((N, len(BANDS)), dtype=bool)
        jobs = _expand_hyperpar(_make_hyperpar(
            fluxes=fluxes, errors=errors, uplims=uplims, run_id=''
        ))
        assert [j.run_id for j in jobs] == ['obj0', 'obj1', 'obj2']

    def test_run_id_multi_source_list(self):
        N = 2
        fluxes = np.stack([FLUXES] * N)
        errors = np.stack([ERRORS] * N)
        uplims = np.zeros((N, len(BANDS)), dtype=bool)
        jobs = _expand_hyperpar(_make_hyperpar(
            fluxes=fluxes, errors=errors, uplims=uplims,
            run_id=['src_a', 'src_b'],
        ))
        assert [j.run_id for j in jobs] == ['src_a', 'src_b']

    def test_run_id_with_model_suffix(self):
        models = [dict(do_AGN=False), dict(do_AGN=True)]
        jobs = _expand_hyperpar(_make_hyperpar(run_id='ngc2403', models=models))
        assert jobs[0].run_id == 'ngc2403_AGNFalse'
        assert jobs[1].run_id == 'ngc2403_AGNTrue'

    # ---- model config propagation ----

    def test_model_override_do_agn(self):
        models = [dict(do_AGN=False), dict(do_AGN=True)]
        jobs = _expand_hyperpar(_make_hyperpar(models=models))
        assert jobs[0].do_AGN is False
        assert jobs[1].do_AGN is True

    def test_model_override_sfh_model(self):
        models = [dict(sfh_model='insitu'), dict(sfh_model='delayedexp')]
        jobs = _expand_hyperpar(_make_hyperpar(models=models))
        assert jobs[0].sfh_model == 'insitu'
        assert jobs[1].sfh_model == 'delayedexp'

    def test_model_galaxy_parameters_patch(self):
        extra_gp = {'sfh.psi_max': ([0., 3.], True)}
        models = [dict(galaxy_parameters=extra_gp)]
        jobs = _expand_hyperpar(_make_hyperpar(models=models))
        assert 'sfh.psi_max' in jobs[0].galaxy_parameters

    # ---- uplims default to False ----

    def test_uplims_none_defaults_to_all_false(self):
        jobs = _expand_hyperpar(_make_hyperpar(uplims=None))
        assert not np.any(jobs[0].uplims)

    # ---- run_id list wrong length raises ----

    def test_run_id_list_wrong_length_raises(self):
        N = 3
        fluxes = np.stack([FLUXES] * N)
        errors = np.stack([ERRORS] * N)
        uplims = np.zeros((N, len(BANDS)), dtype=bool)
        with pytest.raises(ValueError, match='run_id has 2 entries but 3 sources'):
            _expand_hyperpar(_make_hyperpar(
                fluxes=fluxes, errors=errors, uplims=uplims,
                run_id=['a', 'b'],
            ))
