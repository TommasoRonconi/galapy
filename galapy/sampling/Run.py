import os
os.environ["OMP_NUM_THREADS"]        = "1"
os.environ["MKL_NUM_THREADS"]        = "1"
os.environ["OPENBLAS_NUM_THREADS"]   = "1"
os.environ["BLAS_NUM_THREADS"]       = "1"
os.environ["NUMEXPR_NUM_THREADS"]    = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
import numpy

import warnings
import argparse
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
}

################################################################################

class PipelineState :
    """Encapsulates all per-object fitting state.

    Replaces the module-level ``global_dict`` so that multiple independent
    fits can coexist in the same Python process (catalogue mode) without
    sharing mutable global state.
    """

    def __init__ ( self, data, model, noise, handler ) :
        self.data    = data
        self.model   = model
        self.noise   = noise
        self.handler = handler

    @classmethod
    def initialize (
            cls,
            bands, fluxes, errors, uplims, filter_args, params,
            sfh_model    = 'insitu', ssp_lib = 'parsec22.NT',
            do_Radio     = False, do_Xray = False, do_AGN = False,
            noise_model  = None, noise_params = {},
            gxy_kwargs   = {}, noise_kwargs = {}, filter_kwargs = {}
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

        return cls( data, model, noise, handler )

################################################################################

def initialize ( bands, fluxes, errors, uplims, filter_args, params,
                 sfh_model = 'insitu', ssp_lib = 'parsec22.NT',
                 do_Radio = False, do_Xray = False, do_AGN = False,
                 noise_model = None, noise_params = {},
                 gxy_kwargs = {}, noise_kwargs = {}, filter_kwargs = {} ) :
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
    )

################################################################################

def loglikelihood ( par, state, **kwargs ) :

    nested = state.handler.return_nested( par )
    try :
        state.model.set_parameters( **nested['galaxy'] )
        flux_model = numpy.asarray( state.model.photoSED() )
    except RuntimeError :
        return -numpy.inf

    if state.noise is not None :
        state.noise.set_parameters( **nested['noise'] )
        errors = state.noise.apply( state.data.errors, flux_model )
        noise_llike = numpy.log( 2 * numpy.pi * errors**2 ).sum()
        if numpy.isnan( noise_llike ) :
            return -numpy.inf
        return gaussian_loglikelihood( data   = state.data.fluxes,
                                       error  = errors,
                                       model  = flux_model,
                                       uplims = state.data.uplims,
                                       **kwargs ) - 0.5 * noise_llike
    return gaussian_loglikelihood( data   = state.data.fluxes,
                                   error  = state.data.errors,
                                   model  = flux_model,
                                   uplims = state.data.uplims,
                                   **kwargs )

################################################################################

def logprob ( par, state, **kwargs ) :

    pmin, pmax = state.handler.par_prior.T
    if all( ( pmin < par ) & ( par < pmax ) ) :
        return loglikelihood( par, state, **kwargs )

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

    if sampler == 'dynesty' :

        from galapy.sampling.Statistics import transform_to_prior_unit_cube

        sampler_kw.update(
            { 'ptform_args' : ( state.handler.par_prior, ),
              'logl_args'   : ( state, ),
              'logl_kwargs' : logl_kw,
              'queue_size'  : Ncpu,
            }
        )
        sampler = Sampler( loglikelihood   = loglikelihood,
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

    else :
        raise ValueError( f'The sampler chosen "{sampler}" is not valid. '
                          'Valid samplers are ["dynesty", "emcee"].' )

    return sampler

################################################################################

def store_results ( state, sampler,
                    out_dir = '.', name = '',
                    method = 'hdf5', lightweight = False,
                    pickle_sampler = False, pickle_raw = False ) :

    outbase = generate_output_base( out_dir = out_dir, name = name )
    _ = dump_results( model      = state.model,
                      handler    = state.handler,
                      data       = state.data,
                      sampler    = sampler,
                      noise      = state.noise,
                      outbase    = outbase,
                      method     = method,
                      lightweight = lightweight )
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
                     pickle_sampler = False, pickle_raw = False ) :

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
    )

    return;

################################################################################

def _sample_parallel ( state, which_sampler = 'dynesty',
                       nwalkers = None, nsamples = None,
                       sampler_kw = {}, logl_kw = {}, run_sampling_kw = {},
                       Ncpu = None,
                       out_dir = '.', name = '',
                       store_method = 'hdf5', store_lightweight = False,
                       pickle_sampler = False, pickle_raw = False ) :
    import multiprocessing as mp

    if Ncpu is None :
        Ncpu = mp.cpu_count()

    if which_sampler == 'dynesty' :

        from galapy.sampling.Statistics import transform_to_prior_unit_cube
        from dynesty.pool import Pool as DynestyPool

        sampling_kw = dict( _default_sampling_kw.get( which_sampler, {} ) )
        sampling_kw.update( run_sampling_kw )

        with DynestyPool(
                Ncpu, loglikelihood, transform_to_prior_unit_cube,
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

        with mp.get_context( 'fork' ).Pool( Ncpu ) as pool :
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
                          'Valid samplers are ["dynesty", "emcee"].' )

    store_results(
        state, sampler,
        out_dir    = out_dir,
        name       = name,
        method     = store_method,
        lightweight = store_lightweight,
        pickle_sampler = pickle_sampler,
        pickle_raw     = pickle_raw,
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
    of N×K single-job SimpleNamespace objects.

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
        Length N×K, where N is the number of sources (inferred from
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
                sampler           = hyperpar.sampler,
                nwalkers          = getattr( hyperpar, 'nwalkers', None ),
                nsamples          = getattr( hyperpar, 'nsamples', None ),
                sampler_kw        = hyperpar.sampler_kw,
                sampling_kw       = hyperpar.sampling_kw,
                output_directory  = hyperpar.output_directory,
                run_id            = job_run_id,
                store_method      = hyperpar.store_method,
                store_lightweight = hyperpar.store_lightweight,
                pickle_sampler    = hyperpar.pickle_sampler,
                pickle_raw        = hyperpar.pickle_raw,
            ) )

    return jobs

################################################################################

def _catalogue_worker ( job ) :
    """Worker for one (source, model-variant) job.

    Receives a fully-resolved :class:`~types.SimpleNamespace` produced by
    :func:`_expand_hyperpar`.  Builds its own :class:`PipelineState` and
    runs a serial fit.  Designed to be called inside a forked/spawned
    subprocess so that thread-count limits take effect before any numpy
    operations.
    """
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
    )

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
        pickle_sampler    = job.pickle_sampler,
        pickle_raw        = job.pickle_raw,
    )

    return job.run_id


def _sample_catalogue ( jobs, Ncpu = None ) :
    """Run one serial fit per job in parallel across worker processes."""
    import multiprocessing as mp
    import sys

    if Ncpu is None :
        Ncpu = min( mp.cpu_count(), len( jobs ) )

    ctx_name = 'fork' if sys.platform.startswith( 'linux' ) else 'forkserver'
    ctx      = mp.get_context( ctx_name )

    with ctx.Pool( Ncpu ) as pool :
        pool.map( _catalogue_worker, jobs )

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
    # Multi-job mode: catalogue-level parallelism

    if len( jobs ) > 1 :
        _sample_catalogue( jobs, Ncpu = args.Ncpu )
        return;

    ####################################################################
    # Single-job mode

    job   = jobs[ 0 ]
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
noise_model = None

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

# Choose the sampler, valid options are
# 'emcee' : Affine Invariant MCMC ensemble sampler
# 'dynesty' : Dynamic Nested Sampler
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
# - emcee   (EnsembleSampler): https://emcee.readthedocs.io/en/stable/user/api/#emcee.EnsembleSampler
# - dynesty (DynamicNestedSampler): https://dynesty.readthedocs.io/en/latest/api.html#dynesty.DynamicNestedSampler
# Example for dynesty — decrease the number of random-walk steps:
# sampler_kw = {{'walks':25}}  # default is 'walks' : 50
# ( 'walks' : >= 50 is recommended for 15 < ndim < 25 )
sampler_kw = {{}}

# Sampling keyword arguments.
# These are the parameters passed to the method that runs the sampling.
# Keys must match the chosen sampler; unrecognised keys will cause an error
# from the underlying library. See the relevant documentation:
# - emcee   (run_mcmc): https://emcee.readthedocs.io/en/stable/user/api/#emcee.EnsembleSampler.run_mcmc
# - dynesty (run_nested): https://dynesty.readthedocs.io/en/latest/api.html#dynesty.DynamicNestedSampler.run_nested
sampling_kw = {{}}

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
