import os
os.environ["OMP_NUM_THREADS"]        = "1"
os.environ["MKL_NUM_THREADS"]        = "1"
os.environ["OPENBLAS_NUM_THREADS"]   = "1"
os.environ["BLAS_NUM_THREADS"]       = "1"
os.environ["NUMEXPR_NUM_THREADS"]    = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
import numpy

import inspect
import warnings
import argparse
import functools
import importlib.util
from types import SimpleNamespace

from galapy.PhotometricSystem import PMS
from galapy.Galaxy import PhotoGXY
from galapy.Noise import CalibrationError
from galapy.Handlers import ModelParameters
from galapy.sampling.Statistics import gaussian_loglikelihood
from galapy.sampling.Sampler import Sampler
from galapy.sampling.Observation import Observation
from galapy.sampling.Results import generate_output_base, dump_results

_default_sampling_kw = {
    'dynesty' : dict(
        dlogz_init=0.05, nlive_init=500, nlive_batch=100,
        maxiter_init=10000, maxiter_batch=1000, maxbatch=10,
        stop_kwargs = {'target_n_effective': int(5e6)}
    ),
    'emcee' : {},
    'nautilus' : { 'f_live' : 0.01, 'n_eff' : 8000,
                   'discard_exploration' : True, 'verbose' : True },
}    

################################################################################

class PipelineState :
    """Encapsulates all per-object fitting state.

    Replaces the module-level ``global_dict`` so that multiple independent
    fits can coexist in the same Python process (catalogue mode) without
    sharing mutable global state.

    The log-likelihood used by the run is carried here as well, so that every
    code path (serial, parallel, catalogue) picks it up from the same place
    without threading an extra argument through the whole call stack.
    """

    def __init__ ( self, data, model, noise, handler, loglikelihood = None ) :
        self.data    = data
        self.model   = model
        self.noise   = noise
        self.handler = handler
        # ``None`` selects the built-in likelihood; anything else is validated
        # here so that a malformed custom likelihood fails at construction
        # rather than inside the sampling loop.
        self.loglikelihood = _resolve_loglikelihood( loglikelihood )

    @classmethod
    def initialize (
            cls,
            bands, fluxes, errors, uplims, filter_args, params,
            sfh_model    = 'insitu', ssp_lib = 'parsec22.NT',
            do_Radio     = False, do_Xray = False, do_AGN = False,
            noise_model  = None, noise_params = {},
            gxy_kwargs   = {}, noise_kwargs = {}, filter_kwargs = {},
            loglikelihood = None
    ) :

        #########################################################################
        # Build photometric system

        pms = PMS( *filter_args, **filter_kwargs )

        #########################################################################
        # Build observation

        data = Observation( bands, fluxes, errors, uplims, pms )

        #########################################################################
        # Build photometric galaxy

        model = PhotoGXY( pms      = data.pms,
                          sfh      = { 'model' : sfh_model },
                          csp      = { 'ssp_lib' : ssp_lib },
                          do_Radio = do_Radio,
                          do_Xray  = do_Xray,
                          do_AGN   = do_AGN,
                          **gxy_kwargs )
        noise = None
        if noise_model is not None :
            if noise_model == 'calibration_error' :
                noise = CalibrationError( **noise_kwargs )
            else :
                warnings.warn( 'The noise model chosen is not valid, noise will be ignored' )

        #########################################################################
        # Build parameters handler

        sample_params = { '.'.join( ['galaxy', key] ) : value
                          for key, value in params.items() }
        if noise is not None :
            sample_params.update( { '.'.join( ['noise', key] ) : value
                                    for key, value in noise_params.items() } )
            handler = ModelParameters( model, noise,
                                       sample_params = sample_params )
        else :
            handler = ModelParameters( model, sample_params = sample_params )

        #########################################################################
        # Set fixed parameters and initial values for free parameters

        init = handler.return_nested()
        try :
            init['galaxy']['age'] = min(
                model.cosmo.age( init['galaxy']['redshift'] ) - 1.0,
                init['galaxy']['age']
            )
            model.set_parameters( **init['galaxy'] )
        except RuntimeError :
            pass
        if noise is not None : noise.set_parameters( **init['noise'] )

        return cls( data, model, noise, handler,
                    loglikelihood = loglikelihood )

################################################################################

def initialize ( bands, fluxes, errors, uplims, filter_args, params,
                 sfh_model = 'insitu', ssp_lib = 'parsec22.NT',
                 do_Radio = False, do_Xray = False, do_AGN = False,
                 noise_model = None, noise_params = {},
                 gxy_kwargs = {}, noise_kwargs = {}, filter_kwargs = {},
                 loglikelihood = None ) :
    """Build and return a :class:`PipelineState` for one object.

    This is a convenience wrapper around :meth:`PipelineState.initialize`
    kept as a free function for API familiarity.

    Returns
    -------
    PipelineState
    """
    return PipelineState.initialize(
        bands, fluxes, errors, uplims, filter_args, params,
        sfh_model    = sfh_model,   ssp_lib      = ssp_lib,
        do_Radio     = do_Radio,    do_Xray      = do_Xray,
        do_AGN       = do_AGN,      noise_model  = noise_model,
        noise_params = noise_params, gxy_kwargs  = gxy_kwargs,
        noise_kwargs = noise_kwargs, filter_kwargs = filter_kwargs,
        loglikelihood = loglikelihood,
    )

################################################################################

def loglikelihood ( par, state, **kwargs ) :
    """The built-in Gaussian log-likelihood of galapy.

    This is also the reference implementation for user-provided likelihoods
    (see the ``loglikelihood`` hyper-parameter of the parameter file): any
    replacement must honour the same contract.

    Parameters
    ----------
    par : array-like
        1-D array with the current values of the free parameters, in the order
        given by ``state.handler.par_free``.
    state : PipelineState
        The state of the run, giving access to ``state.handler`` (parameter
        handler), ``state.model`` (the galaxy model), ``state.noise`` (the
        noise model, or ``None``) and ``state.data`` (the ``Observation``).
    **kwargs
        Extra keyword arguments forwarded by the sampler through ``logl_kw``
        (currently ``method_uplims``). A custom likelihood should always
        accept ``**kwargs``.

    Returns
    -------
    float
        A finite scalar, or ``-numpy.inf`` for parameter sets that are
        rejected — either because the model refuses them (``RuntimeError``
        from ``set_parameters``) or because the likelihood is not finite.
        Returning anything else (an array, a ``nan``) will corrupt the run.
    """

    nested = state.handler.return_nested( par )
    with numpy.errstate( all = 'ignore' ) :
        try :
            state.model.set_parameters( **nested['galaxy'] )
            flux_model = numpy.asarray( state.model.photoSED() )
        except RuntimeError :
            return -numpy.inf

        if state.noise is not None :
            state.noise.set_parameters( **nested['noise'] )
            errors      = state.noise.apply( state.data.errors, flux_model )
            noise_llike = numpy.log( 2 * numpy.pi * errors**2 ).sum()
            if not numpy.isfinite( noise_llike ) :
                return -numpy.inf
            llike = ( gaussian_loglikelihood( data   = state.data.fluxes,
                                              error  = errors,
                                              model  = flux_model,
                                              uplims = state.data.uplims,
                                              **kwargs ) - 0.5 * noise_llike )
        else :
            llike = gaussian_loglikelihood( data   = state.data.fluxes,
                                            error  = state.data.errors,
                                            model  = flux_model,
                                            uplims = state.data.uplims,
                                            **kwargs )

    return llike if numpy.isfinite( llike ) else -numpy.inf

################################################################################

def loglikelihood_name ( func ) :
    """Fully-qualified name identifying a log-likelihood callable.

    Returns the ``module.qualified_name`` string, which is exactly what
    ``pickle`` stores when serialising a function by reference. It is used as
    a provenance marker in the results file so that two runs can be checked
    for likelihood compatibility before their evidences are compared.
    """
    # A functools.partial binds per-object data, not behaviour: unwrap it
    # (repeatedly, partials can nest) so that the marker identifies the
    # underlying function. The bound values are not part of the marker.
    while isinstance( func, functools.partial ) :
        func = func.func
    module   = getattr( func, '__module__',  None ) or '<unknown>'
    qualname = getattr( func, '__qualname__', None )
    if qualname is None :
        # Callable instances carry no ``__qualname__`` of their own: fall back
        # to the class, which is stable across runs (``repr`` would embed the
        # memory address, making two identical runs look incompatible).
        qualname = type( func ).__qualname__
    return f'{module}.{qualname}'

def _warn_if_not_picklable ( func ) :
    """Warn when ``func`` will not survive being sent to a worker process.

    Parallel runs ship the likelihood to the workers by pickling it, which
    stores a *reference* (module plus qualified name) rather than the code
    itself. Lambdas, nested functions and anything defined directly inside the
    parameter file — imported by path as the throw-away module
    ``hyper_parameters`` — have no importable path and cannot be rebuilt by a
    spawned worker.
    """
    # A partial is itself pickled by value, storing a reference to the wrapped
    # callable plus the bound arguments: what must be importable is the
    # innermost function, so the checks below apply to that.
    while isinstance( func, functools.partial ) :
        func = func.func

    qualname = getattr( func, '__qualname__', '' ) or ''
    module   = getattr( func, '__module__',   None )

    if qualname == '<lambda>' :
        reason = 'it is a lambda'
    elif '<locals>' in qualname :
        reason = 'it is defined inside another function (a closure)'
    elif module in ( None, 'hyper_parameters' ) :
        reason = ( f'it is defined in the "{module}" namespace, which a worker '
                   'process cannot import' )
    else :
        return

    warnings.warn(
        f'The custom loglikelihood "{loglikelihood_name(func)}" is most likely '
        f'not picklable ({reason}). It will work in serial runs '
        '(galapy-fit --serial) but is expected to fail as soon as the sampling '
        'is parallelised. Define it in a regular, importable module (installed '
        'or reachable through PYTHONPATH) and import it in the parameter file.'
    )

def _resolve_loglikelihood ( func ) :
    """Validate a user-provided log-likelihood and return the callable to use.

    Parameters
    ----------
    func : callable or None
        Value of the ``loglikelihood`` hyper-parameter. ``None`` selects the
        built-in :func:`loglikelihood`.

    Returns
    -------
    callable
        A callable honouring the ``(par, state, **kwargs) -> float`` contract.

    Raises
    ------
    TypeError
        If ``func`` is neither ``None`` nor a callable, or if it cannot be
        called with the two mandatory positional arguments ``(par, state)``.
    """
    if func is None :
        return loglikelihood

    if not callable( func ) :
        raise TypeError(
            'The "loglikelihood" hyper-parameter must be either None (use the '
            'built-in Gaussian likelihood) or a callable with signature '
            f'(par, state, **kwargs), got an object of type '
            f'"{type(func).__name__}".'
        )

    try :
        signature = inspect.signature( func )
    except ( TypeError, ValueError ) :
        # C-implemented callables may not expose a signature: nothing to check.
        signature = None

    if signature is not None :
        try :
            signature.bind( 'par', 'state' )
        except TypeError as err :
            raise TypeError(
                f'The custom loglikelihood "{loglikelihood_name(func)}" cannot '
                'be called as loglikelihood(par, state). It must accept the '
                'two mandatory positional arguments `par` (the free-parameter '
                'vector) and `state` (the PipelineState of the run). '
                f'Original error: {err}'
            ) from err

        if not any( p.kind is inspect.Parameter.VAR_KEYWORD
                    for p in signature.parameters.values() ) :
            warnings.warn(
                f'The custom loglikelihood "{loglikelihood_name(func)}" does '
                'not accept **kwargs. Sampling keyword arguments (currently '
                '"method_uplims") are forwarded to the likelihood and will '
                'raise a TypeError unless they are among its named arguments.'
            )

    _warn_if_not_picklable( func )

    return func

def _nautilus_loglikelihood ( logl, state, par, **kwargs ) :
    """Adapter keeping the likelihood contract positional under nautilus.

    nautilus wraps the likelihood with ``functools.partial( func, *args,
    **kwargs )``, so positional extras are PREPENDED to the sampled vector and
    the state cannot follow ``par`` positionally. Pre-binding ``logl`` and
    ``state`` here (with :func:`functools.partial`, see :func:`sample`) hands
    nautilus a callable of ``par`` alone, while the user's likelihood keeps
    receiving ``( par, state, **kwargs )`` positionally, exactly as with
    dynesty and emcee. Module-level on purpose: the resulting partial is
    picklable by reference for nautilus' internal pool.
    """
    return logl( par, state, **kwargs )

################################################################################

def logprob ( par, state, **kwargs ) :

    pmin, pmax = state.handler.par_prior.T
    if all( ( pmin < par ) & ( par < pmax ) ) :
        return getattr( state, 'loglikelihood', loglikelihood )( par, state,
                                                                 **kwargs )

    return -numpy.inf

################################################################################

def sample ( state, sampler = 'dynesty', nwalkers = None, nsamples = None,
             sampler_kw = {}, logl_kw = {}, run_sampling_kw = {},
             Ncpu = 1, pool = None ) :

    # Sampler keywords
    sampler_kw = dict( sampler_kw )

    # Sampling keywords
    sampling_kw = dict( _default_sampling_kw.get( sampler, {} ) )
    sampling_kw.update( run_sampling_kw )

    # Log-likelihood of the run: the custom one carried by the state, when
    # present, otherwise the built-in default. The emcee branch goes through
    # ``logprob``, which resolves it the same way.
    logl = getattr( state, 'loglikelihood', loglikelihood )

    if sampler == 'dynesty' :

        from galapy.sampling.Statistics import transform_to_prior_unit_cube

        sampler_kw.update(
            { 'ptform_args' : ( state.handler.par_prior, ),
              'logl_args'   : ( state, ),
              'logl_kwargs' : logl_kw,
              'queue_size'  : Ncpu,
            }
        )
        sampler = Sampler( loglikelihood   = logl,
                           ndim            = len( state.handler.par_free ),
                           sampler         = sampler,
                           prior_transform = transform_to_prior_unit_cube,
                           pool            = pool,
                           sampler_kw      = sampler_kw )

        sampler.run_sampling( sampling_kw = sampling_kw )

    elif sampler == 'emcee' :

        sampler_kw.update( { 'args'   : ( state, ),
                              'kwargs' : logl_kw } )
        sampler = Sampler( loglikelihood = logprob,
                           ndim          = len( state.handler.par_free ),
                           sampler       = sampler,
                           nwalkers      = nwalkers,
                           pool          = pool,
                           sampler_kw    = sampler_kw )

        pos_init = state.handler.rng.uniform(
            *state.handler.par_prior.T,
            size = ( nwalkers, len( state.handler.par_free ) )
        )

        sampler.run_sampling( pos_init, nsamples,
                              sampling_kw = sampling_kw )

    elif sampler == 'nautilus' :

        # Define prior-transform to map parameters in the unit-cube
        from galapy.sampling.Statistics import transform_to_prior_unit_cube

        # Build sampler — likelihood_kwargs and prior_kwargs go to the constructor.
        # NOTE: nautilus wraps both callables with
        #   functools.partial( func, *func_args, **func_kwargs )
        # so anything passed through prior_args/likelihood_args is PREPENDED to
        # the sampled vector, reversing the argument order. The prior limits
        # are therefore passed by keyword, while the state is pre-bound to the
        # likelihood through the _nautilus_loglikelihood adapter, which keeps
        # the ( par, state, **kwargs ) contract positional for the (possibly
        # custom) likelihood, exactly as with the other samplers.
        sampler_kw.update(
            { 'prior_kwargs' : { 'prior_limits' : state.handler.par_prior },
              'likelihood_kwargs' : logl_kw,
            }
        )
        sampler = Sampler( loglikelihood = functools.partial(
                               _nautilus_loglikelihood, logl, state ),
                           ndim = len( state.handler.par_free ),
                           sampler = sampler,
                           prior_transform = transform_to_prior_unit_cube,
                           pool = pool,
                           sampler_kw = sampler_kw )

        # Run sampling
        sampler.run_sampling( sampling_kw = sampling_kw )

    else :
        raise ValueError( f'The sampler chosen "{sampler}" is not valid. '
                          'Valid samplers are ["dynesty", "emcee", "nautilus"].' )

    return sampler

################################################################################

def store_results ( state, sampler,
                    out_dir = '.', name = '',
                    method = 'hdf5', lightweight = False,
                    pickle_sampler = False, pickle_raw = False,
                    store_quantities = None ) :

    outbase = generate_output_base( out_dir = out_dir, name = name )
    _ = dump_results( model      = state.model,
                      handler    = state.handler,
                      data       = state.data,
                      sampler    = sampler,
                      noise      = state.noise,
                      outbase    = outbase,
                      method     = method,
                      lightweight = lightweight,
                      derived    = store_quantities,
                      loglikelihood_name = loglikelihood_name(
                          getattr( state, 'loglikelihood', loglikelihood )
                      ) )
    sampler.save_results( outbase        = outbase,
                          pickle_sampler = pickle_sampler,
                          pickle_raw     = pickle_raw )
    return;

################################################################################

def _sample_serial ( state, which_sampler = 'dynesty',
                     nwalkers = None, nsamples = None,
                     sampler_kw = {}, logl_kw = {}, run_sampling_kw = {},
                     out_dir = '.', name = '',
                     store_method = 'hdf5', store_lightweight = False,
                     pickle_sampler = False, pickle_raw = False,
                     store_quantities = None ) :

    sampler = sample(
        state,
        sampler         = which_sampler,
        sampler_kw      = sampler_kw,
        logl_kw         = logl_kw,
        run_sampling_kw = run_sampling_kw,
        nwalkers        = nwalkers,
        nsamples        = nsamples,
    )

    store_results(
        state, sampler,
        out_dir    = out_dir,
        name       = name,
        method     = store_method,
        lightweight = store_lightweight,
        pickle_sampler = pickle_sampler,
        pickle_raw     = pickle_raw,
        store_quantities = store_quantities,
    )

    return;

################################################################################

def _sample_parallel ( state, which_sampler = 'dynesty',
                       nwalkers = None, nsamples = None,
                       sampler_kw = {}, logl_kw = {}, run_sampling_kw = {},
                       Ncpu = None,
                       out_dir = '.', name = '',
                       store_method = 'hdf5', store_lightweight = False,
                       pickle_sampler = False, pickle_raw = False,
                       store_quantities = None ) :
    import multiprocessing as mp

    if Ncpu is None :
        Ncpu = mp.cpu_count()

    if which_sampler == 'dynesty' :

        from galapy.sampling.Statistics import transform_to_prior_unit_cube
        from dynesty.pool import Pool as DynestyPool

        sampling_kw = dict( _default_sampling_kw.get( which_sampler, {} ) )
        sampling_kw.update( run_sampling_kw )

        with DynestyPool(
                Ncpu,
                getattr( state, 'loglikelihood', loglikelihood ),
                transform_to_prior_unit_cube,
                logl_args    = ( state, ),
                logl_kwargs  = logl_kw,
                ptform_args  = ( state.handler.par_prior, ),
        ) as pool :
            sampler = Sampler(
                loglikelihood   = pool.loglike,
                ndim            = len( state.handler.par_free ),
                sampler         = which_sampler,
                prior_transform = pool.prior_transform,
                pool            = pool,
                sampler_kw      = sampler_kw,
            )
            sampler.run_sampling( sampling_kw = sampling_kw )

    elif which_sampler == 'emcee' :

        import sys
        #_ctx = 'fork' if sys.platform.startswith( 'linux' ) else 'spawn'
        _ctx = 'spawn'
        with mp.get_context( _ctx ).Pool( Ncpu ) as pool :
            sampler = sample(
                state,
                sampler         = which_sampler,
                sampler_kw      = sampler_kw,
                logl_kw         = logl_kw,
                run_sampling_kw = run_sampling_kw,
                nwalkers        = nwalkers,
                nsamples        = nsamples,
                Ncpu            = Ncpu,
                pool            = pool,
            )

    elif which_sampler == 'nautilus' :

        # nautilus parallelises the likelihood calls through the pool it is
        # handed, so the pool is simply forwarded to sample(). What travels to
        # the workers is the partial built there, which carries the state
        # pre-bound to the likelihood through _nautilus_loglikelihood.
        _ctx = 'spawn'
        with mp.get_context( _ctx ).Pool( Ncpu ) as pool :
            sampler = sample(
                state,
                sampler         = which_sampler,
                sampler_kw      = sampler_kw,
                logl_kw         = logl_kw,
                run_sampling_kw = run_sampling_kw,
                nwalkers        = nwalkers,
                nsamples        = nsamples,
                Ncpu            = Ncpu,
                pool            = pool,
            )

    else :
        raise ValueError( f'The sampler chosen "{which_sampler}" is not valid. '
                          'Valid samplers are ["dynesty", "emcee", "nautilus"].' )

    store_results(
        state, sampler,
        out_dir    = out_dir,
        name       = name,
        method     = store_method,
        lightweight = store_lightweight,
        pickle_sampler = pickle_sampler,
        pickle_raw     = pickle_raw,
        store_quantities = store_quantities,
    )

    return;

################################################################################

def _model_suffixes ( resolved_variants ) :
    """Return one name suffix per model variant, encoding the varying keys.

    Only keys that differ across variants are included, e.g. 'AGNTrue',
    'SFHinsitu_AGNFalse'.
    """
    _KEYS = [
        ( 'sfh_model', 'SFH'   ),
        ( 'do_AGN',    'AGN'   ),
        ( 'do_Radio',  'Radio' ),
        ( 'do_Xray',   'Xray'  ),
    ]
    varying = [
        ( key, pfx ) for key, pfx in _KEYS
        if len( set( str( r[key] ) for r in resolved_variants ) ) > 1
    ]
    return [
        '_'.join( f'{pfx}{r[key]}' for key, pfx in varying )
        for r in resolved_variants
    ]

################################################################################

def _expand_hyperpar ( hyperpar ) :
    """Expand a (possibly multi-source, multi-model) hyperpar into a flat list
    of NxK single-job SimpleNamespace objects.

    Each returned object is a fully-resolved specification for one
    (source, model-variant) pair — equivalent to a single-object parameter
    file ready for :func:`_catalogue_worker`.

    Parameters
    ----------
    hyperpar : module-like namespace
        The parameter file loaded as a Python module by :func:`_run`.

    Returns
    -------
    list of SimpleNamespace
        Length NxK, where N is the number of sources (inferred from
        ``hyperpar.fluxes.shape[0]`` when 2-D, else 1) and K is
        ``len(hyperpar.models)`` when present, else 1.
    """
    # ------------------------------------------------------------------ #
    # Step 1 — detect N                                                  #
    # ------------------------------------------------------------------ #
    raw_fluxes = numpy.asarray( hyperpar.fluxes )
    if raw_fluxes.ndim == 1 :
        N         = 1
        fluxes_2d = raw_fluxes[ numpy.newaxis, : ]
        errors_2d = numpy.asarray( hyperpar.errors )[ numpy.newaxis, : ]
        raw_uplims = hyperpar.uplims
        uplims_2d = (
            numpy.zeros( ( 1, raw_fluxes.shape[0] ), dtype = bool )
            if raw_uplims is None
            else numpy.asarray( raw_uplims )[ numpy.newaxis, : ]
        )
    elif raw_fluxes.ndim == 2 :
        N         = raw_fluxes.shape[0]
        fluxes_2d = raw_fluxes
        errors_2d = numpy.asarray( hyperpar.errors )
        raw_uplims = hyperpar.uplims
        uplims_2d = (
            numpy.zeros_like( fluxes_2d, dtype = bool )
            if raw_uplims is None
            else numpy.asarray( raw_uplims )
        )
    else :
        raise ValueError(
            f'fluxes must be a 1-D or 2-D array-like, got shape {raw_fluxes.shape}.'
        )

    # ------------------------------------------------------------------ #
    # Step 2 — detect K                                                  #
    # ------------------------------------------------------------------ #
    raw_models     = getattr( hyperpar, 'models', None ) or []
    model_variants = list( raw_models ) if raw_models else [ {} ]
    K              = len( model_variants )

    # ------------------------------------------------------------------ #
    # Step 2b — the log-likelihood is a top-level, run-wide setting       #
    # ------------------------------------------------------------------ #
    # Validate once, here, so that a malformed custom likelihood is reported
    # in the main process before any job is built.
    custom_logl = getattr( hyperpar, 'loglikelihood', None )
    _ = _resolve_loglikelihood( custom_logl )

    # It must NOT be overridable per model variant: the likelihood is the term
    # carrying the data, and the whole point of running K variants on the same
    # source is to compare their evidences. Evidences computed with different
    # likelihoods do not form a Bayes factor about the models -- the ratio
    # would mostly measure which likelihood assigns more probability mass to
    # the dataset. This is a silent scientific error, hence a hard failure.
    for j, mv in enumerate( model_variants ) :
        if 'loglikelihood' in mv :
            raise ValueError(
                f'models[{j}] specifies a per-variant "loglikelihood". The '
                'log-likelihood must be identical across model variants: '
                'evidences obtained with different likelihoods are not '
                'comparable and their ratio is not a Bayes factor about the '
                'models. Set "loglikelihood" once at the top level of the '
                'parameter file. To genuinely compare two likelihoods, run '
                'them from separate parameter files.'
            )

    # ------------------------------------------------------------------ #
    # Step 3 — validate and split per-source galaxy_parameters overrides #
    # A plain list or 1-D numpy array (not a tuple) signals per-source   #
    # fixed values.  Tuples are prior specs and are left untouched.      #
    # ------------------------------------------------------------------ #
    gp_base    = dict( hyperpar.galaxy_parameters )
    per_source = {}
    for key, val in list( gp_base.items() ) :
        if isinstance( val, ( list, numpy.ndarray ) ) and not isinstance( val, tuple ) :
            arr = numpy.asarray( val )
            if arr.ndim == 1 :
                if N == 1 :
                    raise ValueError(
                        f"galaxy_parameters['{key}'] is a list but only 1 source "
                        f"was detected from fluxes shape. "
                        f"Use a scalar for a single-source fixed parameter."
                    )
                if len( arr ) != N :
                    raise ValueError(
                        f"galaxy_parameters['{key}'] has {len(arr)} entries but "
                        f"{N} sources were detected from fluxes shape."
                    )
                per_source[ key ] = arr
                del gp_base[ key ]

    # ------------------------------------------------------------------ #
    # Step 4 — run_id list                                               #
    # ------------------------------------------------------------------ #
    raw_run_id = getattr( hyperpar, 'run_id', '' ) or ''
    if isinstance( raw_run_id, ( list, numpy.ndarray ) ) :
        run_ids = [ str( r ) for r in raw_run_id ]
        if len( run_ids ) != N :
            raise ValueError(
                f'run_id has {len(run_ids)} entries but {N} sources detected.'
            )
    elif raw_run_id :
        run_ids = (
            [ str( raw_run_id ) ] if N == 1
            else [ f'{raw_run_id}{i:d}' for i in range( N ) ]
        )
    else :
        run_ids = [ '' ] if N == 1 else [ f'obj{i:d}' for i in range( N ) ]

    # ------------------------------------------------------------------ #
    # Step 5 — model suffixes (only when K > 1)                          #
    # ------------------------------------------------------------------ #
    if K > 1 :
        resolved = [
            dict(
                sfh_model = mv.get( 'sfh_model', hyperpar.sfh_model ),
                do_AGN    = mv.get( 'do_AGN',    hyperpar.do_AGN    ),
                do_Radio  = mv.get( 'do_Radio',  hyperpar.do_Radio  ),
                do_Xray   = mv.get( 'do_Xray',   hyperpar.do_Xray   ),
            )
            for mv in model_variants
        ]
        model_sfx = _model_suffixes( resolved )
    else :
        model_sfx = [ '' ]

    # ------------------------------------------------------------------ #
    # Step 6 — build NxK job list                                        #
    # ------------------------------------------------------------------ #
    filters      = getattr( hyperpar, 'filters',        hyperpar.bands )
    filters_cust = getattr( hyperpar, 'filters_custom', None )

    jobs = []
    for i in range( N ) :
        src_gp = dict( gp_base )
        for key, arr in per_source.items() :
            src_gp[ key ] = float( arr[ i ] )

        for j, mv in enumerate( model_variants ) :
            job_gp = dict( src_gp )
            job_gp.update( mv.get( 'galaxy_parameters', {} ) )

            sfx        = model_sfx[ j ]
            job_run_id = f'{run_ids[i]}_{sfx}' if sfx else run_ids[ i ]

            jobs.append( SimpleNamespace(
                bands             = hyperpar.bands,
                fluxes            = fluxes_2d[ i ],
                errors            = errors_2d[ i ],
                uplims            = uplims_2d[ i ],
                filters           = filters,
                filters_custom    = filters_cust,
                galaxy_parameters = job_gp,
                sfh_model         = mv.get( 'sfh_model', hyperpar.sfh_model ),
                ssp_lib           = mv.get( 'ssp_lib',   getattr( hyperpar, 'ssp_lib', 'parsec22.NT' ) ),
                do_AGN            = mv.get( 'do_AGN',    hyperpar.do_AGN   ),
                do_Radio          = mv.get( 'do_Radio',  hyperpar.do_Radio ),
                do_Xray           = mv.get( 'do_Xray',   hyperpar.do_Xray  ),
                noise_model       = mv.get( 'noise_model',      hyperpar.noise_model      ),
                noise_parameters  = mv.get( 'noise_parameters', hyperpar.noise_parameters ),
                noise_kwargs      = mv.get( 'noise_kwargs',     hyperpar.noise_kwargs     ),
                lstep             = hyperpar.lstep,
                method_uplims     = hyperpar.method_uplims,
                loglikelihood     = custom_logl,
                sampler           = hyperpar.sampler,
                nwalkers          = getattr( hyperpar, 'nwalkers', None ),
                nsamples          = getattr( hyperpar, 'nsamples', None ),
                sampler_kw        = hyperpar.sampler_kw,
                sampling_kw       = hyperpar.sampling_kw,
                output_directory  = hyperpar.output_directory,
                run_id            = job_run_id,
                store_method      = hyperpar.store_method,
                store_lightweight = hyperpar.store_lightweight,
                store_quantities  = getattr( hyperpar, 'store_quantities', None ),
                pickle_sampler    = hyperpar.pickle_sampler,
                pickle_raw        = hyperpar.pickle_raw,
            ) )

    return jobs

################################################################################

def _catalogue_worker ( job ) :
    """Worker for one (source, model-variant) job.

    Receives a fully-resolved :class:`~types.SimpleNamespace` produced by
    :func:`_expand_hyperpar`.  Builds its own :class:`PipelineState` and
    runs a fit.  Designed to be called inside a forked/spawned subprocess.

    If ``job.cpu_set`` is present and ``os.sched_setaffinity`` is available
    (Linux), the process is pinned to those CPUs before any computation so
    that forked grandchildren (the inner sampler pool) inherit the same mask.
    If ``job.cpus_per_job > 1`` the inner sampler is run in parallel across
    those CPUs; otherwise it runs serially.
    """
    if hasattr( os, 'sched_setaffinity' ) :
        cpu_set = getattr( job, 'cpu_set', None )
        if cpu_set is not None :
            os.sched_setaffinity( 0, cpu_set )

    for _var in ( 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                  'BLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
                  'VECLIB_MAXIMUM_THREADS' ) :
        os.environ[ _var ] = '1'

    state = PipelineState.initialize(
        job.bands,
        job.fluxes,
        job.errors,
        job.uplims,
        job.filters,
        job.galaxy_parameters,
        sfh_model    = job.sfh_model,
        ssp_lib      = job.ssp_lib,
        do_Radio     = job.do_Radio,
        do_Xray      = job.do_Xray,
        do_AGN       = job.do_AGN,
        noise_model  = job.noise_model,
        noise_params = job.noise_parameters,
        gxy_kwargs   = { 'lstep' : job.lstep },
        noise_kwargs = job.noise_kwargs,
        filter_kwargs = job.filters_custom or {},
        loglikelihood = getattr( job, 'loglikelihood', None ),
    )

    cpus_per_job = getattr( job, 'cpus_per_job', 1 )
    _run_fn      = _sample_parallel if cpus_per_job > 1 else _sample_serial
    _extra       = { 'Ncpu' : cpus_per_job } if cpus_per_job > 1 else {}

    _run_fn(
        state,
        which_sampler     = job.sampler,
        nwalkers          = job.nwalkers,
        nsamples          = job.nsamples,
        sampler_kw        = job.sampler_kw,
        logl_kw           = { 'method_uplims' : job.method_uplims },
        run_sampling_kw   = job.sampling_kw,
        out_dir           = job.output_directory,
        name              = job.run_id,
        store_method      = job.store_method,
        store_lightweight = job.store_lightweight,
        store_quantities  = getattr( job, 'store_quantities', None ),
        pickle_sampler    = job.pickle_sampler,
        pickle_raw        = job.pickle_raw,
        **_extra,
    )

    return job.run_id


def _sample_catalogue ( jobs, Ncpu = None ) :
    """Distribute NxK jobs across MPI ranks (future multi-node entry point).

    Not called by :func:`_run` in the current single-node implementation.
    Reserved for the forthcoming MPI path where each rank receives a slice of
    jobs and runs them sequentially with the full per-rank CPU budget, matching
    the strategy in :func:`_run`.

    The ``ProcessPoolExecutor`` + CPU-binding infrastructure below is retained
    as a functional fallback for testing and as scaffolding for the MPI
    integration.
    """
    import multiprocessing as mp
    import sys

    if hasattr( os, 'sched_getaffinity' ) :
        avail = sorted( os.sched_getaffinity( 0 ) )
    else :
        avail = list( range( mp.cpu_count() ) )

    if Ncpu is not None :
        avail = avail[ :Ncpu ]

    n_parallel   = min( len( jobs ), len( avail ) )
    cpus_per_job = len( avail ) // n_parallel

    for i, job in enumerate( jobs ) :
        slot         = i % n_parallel
        job.cpu_set      = set( avail[ slot * cpus_per_job :
                                       ( slot + 1 ) * cpus_per_job ] )
        job.cpus_per_job = cpus_per_job

    _binding = ( 'active' if hasattr( os, 'sched_setaffinity' )
                 else 'not available on this platform' )
    _inner   = ( f'parallel ({cpus_per_job} CPUs/worker)'
                 if cpus_per_job > 1 else 'serial' )
    print( f'galapy-fit  [{len(jobs)} job(s)]  '
           f'{n_parallel} worker(s) x {cpus_per_job} CPU(s)  |  '
           f'inner sampler: {_inner}  |  CPU binding: {_binding}',
           flush = True )

    from concurrent.futures import ProcessPoolExecutor

    # ctx_name = 'fork' if sys.platform.startswith( 'linux' ) else 'forkserver'
    ctx_name = 'spawn'
    ctx      = mp.get_context( ctx_name )

    # ProcessPoolExecutor workers are non-daemon, so each catalogue worker
    # can itself spawn the inner sampler pool without hitting the
    # "daemonic processes are not allowed to have children" restriction.
    with ProcessPoolExecutor( max_workers = n_parallel,
                              mp_context  = ctx ) as executor :
        list( executor.map( _catalogue_worker, jobs ) )

    return;

################################################################################

def _run () :

    ####################################################################
    # Read command-line arguments:

    parser = argparse.ArgumentParser( description = 'options' )
    parser.add_argument( 'parameter_file',
                         default = "",
                         help = 'Path to parameter file' )
    parser.add_argument( '-s', '--serial',
                         dest    = 'serial',
                         action  = 'store_true',
                         help    = 'Run the program serially.' )
    parser.add_argument( '-mp', '--multiprocessing',
                         dest    = 'Ncpu',
                         type    = int,
                         default = None,
                         help    = ( 'flag for shared-memory parallel run: '
                                     'provide number of processes. '
                                     'If none is passed defaults to all the available CPUs.' ) )
    args = parser.parse_args()

    spec    = importlib.util.spec_from_file_location( "hyper_parameters",
                                                      args.parameter_file )
    hyperpar = importlib.util.module_from_spec( spec )
    spec.loader.exec_module( hyperpar )

    ####################################################################
    # Expand to NxK jobs

    jobs = _expand_hyperpar( hyperpar )

    ####################################################################
    # Sequential outer loop — each job receives the full CPU budget.
    # For multi-node cluster runs _sample_catalogue (MPI) will be used
    # instead; that path is not wired here yet.

    if len( jobs ) > 1 :
        _ncpu  = args.Ncpu if args.Ncpu is not None else os.cpu_count()
        _inner = 'serial' if args.serial else f'parallel ({_ncpu} CPUs)'
        print( f'galapy-fit  [{len(jobs)} job(s)]  sequential  |  '
               f'inner sampler: {_inner}',
               flush = True )

    for job in jobs :

        state = initialize(
            job.bands,
            job.fluxes,
            job.errors,
            job.uplims,
            job.filters,
            job.galaxy_parameters,
            sfh_model    = job.sfh_model,
            ssp_lib      = job.ssp_lib,
            do_Radio     = job.do_Radio,
            do_Xray      = job.do_Xray,
            do_AGN       = job.do_AGN,
            noise_model  = job.noise_model,
            noise_params = job.noise_parameters,
            gxy_kwargs   = { 'lstep' : job.lstep },
            noise_kwargs = job.noise_kwargs,
            filter_kwargs = job.filters_custom if job.filters_custom is not None else {},
            loglikelihood = getattr( job, 'loglikelihood', None ),
        )

        if args.serial :
            _sample_serial(
                state,
                which_sampler     = job.sampler,
                nwalkers          = job.nwalkers,
                nsamples          = job.nsamples,
                sampler_kw        = job.sampler_kw,
                logl_kw           = { 'method_uplims' : job.method_uplims },
                run_sampling_kw   = job.sampling_kw,
                out_dir           = job.output_directory,
                name              = job.run_id,
                store_method      = job.store_method,
                store_lightweight = job.store_lightweight,
                store_quantities  = getattr( job, 'store_quantities', None ),
                pickle_sampler    = job.pickle_sampler,
                pickle_raw        = job.pickle_raw,
            )
        else :
            _sample_parallel(
                state,
                which_sampler     = job.sampler,
                nwalkers          = job.nwalkers,
                nsamples          = job.nsamples,
                sampler_kw        = job.sampler_kw,
                logl_kw           = { 'method_uplims' : job.method_uplims },
                run_sampling_kw   = job.sampling_kw,
                Ncpu              = args.Ncpu,
                out_dir           = job.output_directory,
                name              = job.run_id,
                store_method      = job.store_method,
                store_lightweight = job.store_lightweight,
                store_quantities  = getattr( job, 'store_quantities', None ),
                pickle_sampler    = job.pickle_sampler,
                pickle_raw        = job.pickle_raw,
            )

    return;

################################################################################
# This commented out here because reasons ...
# if __name__ == '__main__' :
#     run()
################################################################################

_obs_params = {
    'single' : """
# The observed dataset is expressed as 4 array-likes containing respectively:
# - bands: a list of strings corresponding to the unique identifiers also used
#          in the 'filters' variable below
# - fluxes: measures of fluxes (or upper limits) at the corresponding bands
#           listed in variable 'bands'.
#           The input measure should be given in units of milli-Jansky
# - errors: 1-sigma error on the measures of the fluxes (or upper limits) listed
#           in the variable 'fluxes'
# - uplims: sequence of booleans identifying whether the values listed
#           in argument ``fluxes`` should be considered non-detection (``True``)
#           or a detection (``False``)
bands  = None
fluxes = None
errors = None
uplims = None

# The photometric filters system used in the observation.
# This parameter should be an iterable containing names
# of filters already present in the database, e.g.
#
# filters = ['GOODS.b', 'GOODS.i', 'GOODS.v', 'GOODS.z']
#
# NOTE that, if the bands listed in variable ``bands``
# are all present in the database and have the same names
# of the filters listed with ``galapy.PhotometricSystem.print_filters()``,
# this variable can be set as ``filters = bands``
filters = list(bands)

# Eventual custom photometric filters.
# This parameter should be a nested dictionary with user-defined transmissions.
# Such transmissions are passed as properly formatted dictionaries
#
# filters_custom = { 'filter_1' : { 'wavelengths' : array-like,
#                                   'photons' : array-like },
#                    'filter_2' : { 'wavelengths' : array-like,
#                                   'photons' : array-like },
#                    ...
#                    'filter_N' : { 'wavelengths' : array-like,
#                                   'photons' : array-like } }
#
# As the keywords in the lower level dictionaries suggest,
# the two arrays provided, must define the
# wavelength grid and the corresponding transmission in photon units
# (and thus they must have the same size).
# Note that the chosen keyword in the higher level dictionary
# (e.g. 'filter_1') will be used as unique identifier of
# the custom transmission.
filters_custom = None

# Method for treatment of the upper limits, when present.
# Available methods are:
# - 'simple' : a heaviside function which is 0. when the model's flux is
#              smaller than the upper-limit and +infty otherwise
# - 'chi2' : the distance between model and data is expressed as a normal chi-squared
# - 'S12' : modified version of the chi-squared integrating a gaussian, with mean = data-flux
#           and std = data-error, up to the value of the model's flux (Sawicki, 2012)
method_uplims = 'chi2'
""",
    'catalogue' : """
# Multi-source catalogue run.
#
# bands: list of M filter-name strings, shared across all N sources.
bands = None

# fluxes, errors, uplims: 2-D array-likes of shape (N, M).
#   fluxes : measured flux in milli-Jansky
#   errors : 1-sigma flux error in milli-Jansky
#   uplims : bool, True = upper limit / non-detection (default: all False)
#
# For a single source these may also be 1-D arrays of length M.
#
# Example loading from a numpy archive:
#   import numpy
#   _data  = numpy.load('survey.npz')
#   fluxes = _data['fluxes']   # shape (N, M)
#   errors = _data['errors']   # shape (N, M)
#   uplims = _data['uplims']   # shape (N, M), dtype bool
fluxes = None
errors = None
uplims = None

# The photometric filter system used in the observations.
# For a homogeneous catalogue this can simply mirror 'bands'.
filters = list(bands)

# Eventual custom photometric filters (shared across all sources).
filters_custom = None

# Method for treatment of upper limits.
# Available methods: 'simple', 'chi2', 'S12'
method_uplims = 'chi2'

# Model variants to test on each source (optional).
# Each dict may override any of: sfh_model, do_AGN, do_Radio, do_Xray,
# ssp_lib, noise_model, noise_parameters, noise_kwargs, and/or
# galaxy_parameters (a patch dict merged on top of the shared parameters
# defined below).
# If absent or empty, a single model (the shared configuration) is run.
#
# N x K jobs are launched -- one per (source, variant) pair.
# Output filenames get an auto-generated suffix encoding the varying keys,
# e.g. 'obj0_AGNFalse', 'obj0_AGNTrue', 'obj1_SFHinsitu_AGNFalse', ...
#
# Example -- test with and without AGN for every source:
#   models = [
#       dict(do_AGN=False),
#       dict(do_AGN=True),
#   ]
#
# Example -- compare two SFH topologies.
# Place topology-invariant parameters (age, redshift, ISM, ...) in the
# shared galaxy_parameters block below.  Put SFH-specific parameters
# only in each variant's galaxy_parameters; otherwise the handler warns
# about parameters not present in the other variant's model:
#   models = [
#       dict(sfh_model='insitu',
#            galaxy_parameters={
#                'sfh.psi_max'  : ([0., 4.], True),
#                'sfh.tau_star' : ([6., 11.], True),
#            }),
#       dict(sfh_model='delayedexp',
#            galaxy_parameters={
#                'sfh.psi_norm' : ([0., 4.], True),
#                'sfh.k_shape'  : ([0., 5.], False),
#                'sfh.tau_star' : ([6., 11.], True),
#            }),
#   ]
models = []
""",
}

_run_id_params = {
    'single' : """
# An identification name, it will be pre-pended to all files stored in the output directory
# (if the string is empty the current date+time will be used)
run_id = ''
""",
    'catalogue' : """
# Run IDs for the N sources (list of N strings).
# Each entry will be the prefix of the output files for that source.
# If left empty, sources are labelled automatically ('obj0', 'obj1', ...).
# When K > 1 model variants are specified, an auto-generated suffix is
# appended to distinguish the runs, e.g. 'obj0_AGNFalse', 'obj0_AGNTrue'.
run_id = []
""",
}

default_parameter_file = """
##################################################
# Parameters for building the observation to fit #
##################################################
{2:s}
########################################
# Parameters defining the galaxy model #
########################################

# Cosmological model to use for computing distances and ages.
# Pre-computed cosmologies present in the database:
# - WMAP7
# - WMAP9
# - Planck15
# - Planck18
# To use another cosmological model the user has to provide instead
# a dictionary with the following key-value couples:
# - 'redshift' : an iterable containing monothonically increasing
#                values of redshift
# - 'luminosity_distance' : an iterable with the luminosity distances
#                           corresponding to the redshift values in
#                           the first iterable
# - 'age' : an iterable with the age of the universe corresponding to
#           the redshift values in the first iterable
cosmo     = 'Planck18'

# The Star-Formation History model to use for building the galaxy model.
# Available models are:
# - 'insitu'
# - 'constant'
# - 'delayedexp'
# - 'lognormal'
sfh_model = '{0:s}'

# The SSP library to use for building the unattenuated stellar emission component
ssp_lib   = 'parsec22.NT'

# Whether to provide X-ray emission support (True) or not (False).
do_Xray = False

# Whether to provide radio-emission support (True) or not (False).
do_Radio = False

# Whether to build a galaxy containing an AGN (True) or not (False).
do_AGN = False

# Sub-sampling of the wavelength grid.
# If lstep is an integer it will consider a wavelength grid entry every lstep values.
# If lstep is a sequence of integers or a mask, only the wavelength grid entries
# corresponding to the indices provided will be considered.
# If None, it will consider the whole wavelength grid (safest choice)
lstep = None

# Eventual noise model to add. The default is ``None``, i.e. no noise will be added.
# Valid choices for the noise are:
# - 'calibration_error' : accounts for eventual systematics in the calibration of
#                         the errors in the observations.
#
# If this is set to ``None`` all the other hyper-parameters in this file
# related to noise will be ignored.
noise_model = 'calibration_error'

# Eventual keyword arguments to be passed to the noise model of choice
# (leave empty for no keyword arguments)
noise_kwargs = {{}}

#############################################
# Here define the fixed and free parameters #
#############################################

# The dictionary provided here will be used to set the fixed parameters to
# the given value or to flag parameters as 'free' with some associated prior.
# All the parameters that are not used in the specified galaxy model will be ignored.
# (e.g. if the galaxy model has been built with ``do_AGN = False`` all the eventual
#  AGN-related parameters provided will be ignored)
#
# - To set the parameter 'fixed_parameter' as FIXED provide a single value:
#   parameters = {{
#       ...
#       'fixed_parameter' : some_float_value,
#       ...
#   }}
#
# - To set the parameter 'free_parameter' as FREE provide a tuple
#   parameters = {{
#       ...
#       'free_parameter' : ( a_list, a_bool ),
#       ...
#   }}
#   The list ``a_list`` contains the minimum and maximum value of the
#   UNIFORM prior from which to draw samples for the 'free_parameter'.
#   The boolean ``a_bool`` states whether the prior has to be considered
#   * logarithmic: ``a_bool = True``, therefore samples will be drawn from the interval
#                  10**min(a_list) < 'free_parameter' < 10**max(a_list)
#   * linear: ``a_bool = False``, therefore samples will be drawn from the interval
#             min(a_list) < 'free_parameter' < max(a_list)

galaxy_parameters = {{

    ##########
    # Galaxy #
    ##########

    'age'      : ( [6., 11.], True ),
    'redshift' : ( [0., 10.], False ),

    ##########################
    # Star Formation History #
    ##########################

    'sfh.tau_quench' : ( [6., 11.], True ),
{1:s}
    ########################
    # Inter-Stellar Medium #
    ########################

    'ism.f_MC' : ( [0., 1.], False ),

    ##################
    # Molecular Clouds

    'ism.norm_MC' : 100.,
    'ism.N_MC'    : ( [0., 5.], True ),
    'ism.R_MC'    : ( [0., 5.], True ),
    'ism.tau_esc' : ( [4., 8.], True ),
    'ism.dMClow'  : 1.3,
    'ism.dMCupp'  : 1.6,

    ##############
    # Diffuse Dust

    'ism.norm_DD' : 1.0,
    'ism.Rdust'   : ( [0., 5.], True ),
    'ism.f_PAH'   : ( [0., 1.], False ),
    'ism.dDDlow'  : 0.7,
    'ism.dDDupp'  : 2.0,

    ###############
    # Synchrotron #
    ###############

    'syn.alpha_syn'   : 0.75,
    'syn.nu_self_syn' : 0.2,

    #####################
    # Nebular Free-Free #
    #####################

    'nff.Zgas' : 0.02,
    'nff.Zi'   : 1.,

    ###########################
    # Active Galactic Nucleus #
    ###########################

    'agn.fAGN' : ( [-3., 3.], True ),

    # Template-selecting parameters: these select one of the 24000 pre-computed
    # Fritz+2006 templates by snapping to the nearest discrete grid value.
    # Sampling over these parameters is discouraged.
    'agn.template.ct' : 40,
    'agn.template.al' : 0.,
    'agn.template.be' : -0.5,
    'agn.template.ta' : 6.,
    'agn.template.rm' : 60,
    'agn.template.ia' : 0.001,

}}

noise_parameters = {{

    #########
    # Noise #
    #########

    ###################
    # Calibration Error

    'f_cal' : ( [-10., 1.], True ),
}}

##############################
# Parameters for the sampler #
##############################

# Choose the sampler. Valid options are:
# 'emcee'    : Affine Invariant MCMC ensemble sampler
# 'dynesty'  : Dynamic Nested Sampler
# 'nautilus' : Neural Network-Boosted Nested Sampler
sampler = 'dynesty'

# EMCEE SAMPLER-SPECIFIC MANDATORY PARAMETERS
# - set the number of walkers (``nwalkers``)
# - set the chain length (``nsamples``)
nwalkers = 64
nsamples = 4096

# Sampler keyword arguments.
# These are the parameters passed to the constructor of the chosen sampler.
# Keys must match the chosen sampler; unrecognised keys will cause an error
# from the underlying library. See the relevant documentation:
# - emcee    (EnsembleSampler)     : https://emcee.readthedocs.io/en/stable/user/api/#emcee.EnsembleSampler
# - dynesty  (DynamicNestedSampler): https://dynesty.readthedocs.io/en/latest/api.html#dynesty.DynamicNestedSampler
# - nautilus (Sampler)             : https://nautilus-sampler.readthedocs.io/en/latest/api.html#nautilus.Sampler
# Example for dynesty -> decrease the number of random-walk steps:
# sampler_kw = {{'walks':25}}  # default is 'walks' : 50
# ( 'walks' : >= 50 is recommended for 15 < ndim < 25 )
sampler_kw = {{}}

# Sampling keyword arguments.
# These are the parameters passed to the method that runs the sampling.
# Keys must match the chosen sampler; unrecognised keys will cause an error
# from the underlying library. See the relevant documentation:
# - emcee    (run_mcmc)  : https://emcee.readthedocs.io/en/stable/user/api/#emcee.EnsembleSampler.run_mcmc
# - dynesty  (run_nested): https://dynesty.readthedocs.io/en/latest/api.html#dynesty.DynamicNestedSampler.run_nested
# - nautilus (run)       : https://nautilus-sampler.readthedocs.io/en/latest/api.html#nautilus.Sampler.run
sampling_kw = {{}}

# Custom log-likelihood (advanced).
#
# When None (the default) galapy uses its built-in Gaussian log-likelihood,
# i.e. galapy.sampling.Run.loglikelihood.
#
# To use your own, set this to a callable with the signature
#
#   def my_loglikelihood ( par, state, **kwargs ) :
#       ...
#       return llike
#
# where
# - par    : 1-D array with the current values of the free parameters, ordered
#            as state.handler.par_free
# - state  : the galapy.sampling.Run.PipelineState of the run, giving access to
#            state.handler (parameter handler), state.model (the galaxy model),
#            state.noise (the noise model, or None) and state.data (the
#            Observation being fitted)
# - kwargs : extra keyword arguments forwarded by the sampler (currently
#            'method_uplims'); always accept **kwargs
#
# The function MUST return a scalar, and MUST return -numpy.inf both for
# parameter sets the model rejects (set_parameters raises RuntimeError) and
# whenever the result is not finite. Evaluate the model inside a
# numpy.errstate(all='ignore') context, as the built-in likelihood does.
#
# IMPORTANT: in parallel runs the likelihood is sent to the worker processes by
# reference, so it has to live in a module the workers can import. Define it in
# a separate .py file that is installed or reachable through PYTHONPATH and
# import it here; lambdas, nested functions, and functions defined directly in
# this parameter file will NOT work in parallel. Binding per-object data with
# functools.partial HERE is fine instead: a partial is pickled by value, so
# only the wrapped function needs to be importable.
#
# NOTE: this is a run-wide setting. It cannot be overridden per entry of the
# 'models' list, because evidences computed with different likelihoods are not
# comparable and their ratio is not a Bayes factor about the models.
#
# Example (photometry + external stellar-mass estimate for this object):
#   from functools import partial
#   from my_likelihoods import mstar_loglikelihood
#   loglikelihood = partial( mstar_loglikelihood,
#                            logmstar_obs = 10.65, logmstar_err = 0.15 )
#
# See the "Custom likelihood" how-to in the documentation for a complete,
# copy-pasteable template.
loglikelihood = None

# Output directory (note that if the directory does not exist it will be created)
output_directory = ''

{3:s}
# The method used for storing results.
# Possible choices are:
# - 'hdf5' : uses HDF5 format to save a dictionary containing all the informations
#            to build the Results class used for the analysis and visualization of the
#            sampling-run results. (This is the safest choice, also in terms of security)
# - 'pickle' : uses the standard python pickle format to directly save an instance of
#              the Results class used for the analysis and visualization of the
#              sampling-run results. (Recommended only for local usage, pickled objects
#              are not safe for distribution)
store_method = 'hdf5'

# Only available if the output format chosen is HDF5 (see parameter store_method).
# If True it stores only the chains, weights and loglikelihood values obtained by the
# sampling run (along with information to re-build all the models used in running.
# If False it will save all the infos on the derived quantities as well.
# Selecting lightweight storage, it does not change the way users will ultimatelly access
# the Results object, what will change is the time spent for loading it.
# With lightweight storage the size of the output file is smaller, and the Results object
# is computed (instantiated) when loading the file (so this will take up to some minutes).
# With lightweight storage off, all the quantities are computed at the end of the sampling
# run and are ready to use but the output file can reach a size of up to some GiB.
store_lightweight = False

# Which derived physical quantities to compute from the posterior and store in
# the results file. Choose any subset of the built-in quantities:
#   'SED'   : spectral energy distribution        [flux, mJy ; array]
#   'Mstar' : stellar mass                         [Msun]
#   'Mdust' : dust mass                            [Msun]
#   'Mgas'  : gas mass                             [Msun]
#   'Zstar' : stellar metallicity                  [absolute, metal mass fraction]
#   'Zgas'  : gas metallicity                      [absolute, metal mass fraction]
#   'SFR'   : star-formation rate                  [Msun / yr]
#   'TMC'   : molecular-cloud temperature          [K]
#   'TDD'   : diffuse-dust temperature             [K]
# - None (default) stores all of them.
# - A shorter list makes the output file smaller and post-processing faster.
# - 'SED' is ALWAYS stored, whether or not you list it.
# Example: store_quantities = ['Mstar', 'SFR', 'Mdust']
#
# NOTE: this selects among the built-in quantities only. To store your *own*
# quantities (a colour, a band luminosity, the average attenuation, ...) add
# them after the run with Results.add_property() -- see the
# "Custom derived quantities" how-to in the documentation.
store_quantities = None

# Whether to pickle the sampler raw results.
# (might be useful for analyzing run statistics)
pickle_raw = False

# Whether to pickle the sampler at the end-of-run state.
# (might be useful for extending the run)
pickle_sampler = False

#############################
# ... and that's all folks! #
#############################
"""

def _generate_parameter_file () :

    ####################################################################
    # Read command-line arguments:

    parser = argparse.ArgumentParser(
        description = (
            'Writes the parameter file that has to be modified by users '
            'to meet the needs of their sampling.'
        )
    )
    parser.add_argument( '--name', '-n',
                         dest = 'name',
                         type = str,
                         default = 'galapy_hyper_parameters',
                         help = (
                             'provide here the name you want to give to ' +
                             'the parameter file, you can also choose ' +
                             'a path different to the current working directory: ' +
                             '/path/to/chosen_name\n' +
                             'NOTE THAT the extension ".py" will be appended to ' +
                             'the name chosen (this just to guarantee proper formatting ' +
                             'when opening the file in a text editor). ' +
                             'DEFAULT: ${PWD}/galapy_hyper_parameters.py'
                         ) )
    parser.add_argument( '--SFH_model', '-sfh',
                         dest = 'sfh_model',
                         type = str,
                         default = None,
                         help = (
                             'Choose a SFH model. ' +
                             'Available models are:' +
                             ' insitu, constant, delayedexp, lognormal. ' +
                             'DEFAULT: None'
                         ) )
    parser.add_argument( '--catalogue', '-cat',
                         dest = 'catalogue',
                         action = 'store_true',
                         help = (
                             'Generate a multi-source parameter file. '
                             'fluxes/errors/uplims are described as 2-D '
                             'arrays of shape (N, M); an optional ``models`` '
                             'list enables testing K model variants per source, '
                             'launching N x K parallel jobs in total.'
                         ) )
    args = parser.parse_args()

    ####################################################################

    sfh_models = set(['insitu', 'constant', 'delayedexp', 'lognormal'])
    sfh_params = {
        'insitu' : """
    # In-Situ model
    'sfh.psi_max' : ( [0., 4.], True ),
    'sfh.tau_star' : ( [6., 11.], True ),
        """,
        'constant' : """
    # Constant model
    'sfh.psi'   : ( [0., 4.], True ),
    'sfh.Mdust' : ( [6., 14.], True ),
    'sfh.Zgxy'  : ( [0., 1.], False ),
        """,
        'delayedexp' : """
    # Delayed-Exponential model
    'sfh.psi_norm' : ( [0., 4.], True ),
    'sfh.k_shape'  : ( [0., 5.], False ),
    'sfh.tau_star' : ( [6., 11.], True ),
    'sfh.Mdust' : ( [6., 14.], True ),
    'sfh.Zgxy'  : ( [0., 1.], False ),
        """,
        'lognormal' : """
    # Log-Normal model
    'sfh.psi_norm'   : ( [0., 4.], True ),
    'sfh.sigma_star' : ( [0., 5.], False ),
    'sfh.tau_star'   : ( [6., 11.], True ),
    'sfh.Mdust' : ( [6., 14.], True ),
    'sfh.Zgxy'  : ( [0., 1.], False ),
        """,
    }

    sfh_params_string = ''
    if args.sfh_model not in sfh_models and args.sfh_model is not None :
        raise RuntimeError(
            'The chosen model is not available. '
            'To see a list of the available choices call '
            'this function with the --help argument'
        )
    elif args.sfh_model is None :
        args.sfh_model = 'insitu'
        for sfh in sfh_models :
            sfh_params_string += sfh_params[sfh]
    else :
        sfh_params_string = sfh_params[args.sfh_model]

    obs_mode = 'catalogue' if args.catalogue else 'single'

    with open( args.name + '.py', 'w' ) as paramfile :
        paramfile.write( default_parameter_file.format(
            args.sfh_model, sfh_params_string,
            _obs_params[obs_mode], _run_id_params[obs_mode],
        ) )

    return;

################################################################################
