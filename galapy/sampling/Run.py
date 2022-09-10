import numpy
import os
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import importlib.util

from galapy.Galaxy import PhotoGXY
from galapy.GalaxyParameters import GXYParameters
from galapy.sampling.Statistics import gaussian_loglikelihood
from galapy.sampling.Sampler import Sampler
from galapy.sampling.Observation import Observation
from galapy.internal.utils import set_nested

################################################################################

def initialize ( bands, fluxes, errors, uplims, filters, params,
                 ssp_lib = 'br22.NT', do_Radio = False, do_Xray = False, do_AGN = False,
                 model_kwargs = {} ) :

    #########################################################################
    # Build observation

    data = Observation( bands, fluxes, errors, uplims, filters )
    
    #########################################################################
    # Build photometric galaxy

    model = PhotoGXY( pms = data.pms,
                      csp = { 'ssp_lib' : ssp_lib },
                      do_Radio = do_Radio,
                      do_Xray = do_Xray,
                      do_AGN = do_AGN,
                      **model_kwargs )
    
    #########################################################################
    # Build parameters handler

    handler = GXYParameters( model, params )

    #########################################################################
    # Set fixed parameters and initial values for free parameters

    init = {}
    for klist, v in handler.parameters.items() :
        set_nested( init, klist.split('.'), v )
    model.set_parameters( **init )

    return data, model, handler

################################################################################

def loglikelihood ( par, data, model, handler ) :

    try :
        model.set_parameters( **handler.return_nested( par ) )
        flux_model = numpy.asarray( model.photoSED() )
    except RuntimeError :
        return -numpy.inf
        
    return gaussian_loglikelihood( data   = data.fluxes,
                                   error  = data.errors,
                                   model  = flux_model,
                                   uplims = data.uplims )

def loglikelihood_globals ( par ) :
    return loglikelihood( par, gxy_data, gxy_model, gxy_params )

################################################################################

# Include prior in likelihood (here uniform prior)
def logprob ( par, data, model, handler ) :
    
    pmin, pmax = handler.par_prior.T
    if all( ( pmin < par ) & ( par < pmax ) ) :
        return loglikelihood( par, data, model, handler )
    
    return -numpy.inf

def logprob_globals ( par ) :
    return logprob( par, gxy_data, gxy_model, gxy_params )

################################################################################
        
def sample_serial ( hyperpar ) :

    # Run the initializer, it populates the scope with global variables:
    gxy_data, gxy_model, gxy_params = initialize( bands = hyperpar.bands,
                                                  fluxes = hyperpar.fluxes,
                                                  errors = hyperpar.errors,
                                                  uplims = hyperpar.uplims,
                                                  filters = hyperpar.filters,
                                                  params = hyperpar.parameters,
                                                  ssp_lib = hyperpar.ssp_lib,
                                                  do_Radio = hyperpar.do_Radio,
                                                  do_Xray = hyperpar.do_Xray,
                                                  do_AGN = hyperpar.do_AGN,
                                                  model_kwargs = {
                                                      'lstep' : hyperpar.lstep
                                                  } )

    # It's run-time!
    sampler = None
    if hyperpar.sampler == 'dynesty' :
        
        # Define prior-transform to map parameters in the unit-cube
        from galapy.sampling.Statistics import transform_to_prior_unit_cube

        # Build sampler
        kw = hyperpar.sampler_kw
        kw.update( { 'logl_args' : ( gxy_data, gxy_model, gxy_params ),
                     'ptform_args' : (gxy_params.par_prior,) } ) 
        sampler = Sampler( loglikelihood = loglikelihood,
                           ndim = len( gxy_params.par_free ),
                           sampler = hyperpar.sampler,
                           prior_transform = transform_to_prior_unit_cube,
                           dynesty_sampler_kw = kw )

        # Run sampling
        sampler.run_sampling( dynesty_sampling_kw = hyperpar.sampling_kw )
        
    if hyperpar.sampler == 'emcee' :

        # Build sampler
        kw = hyperpar.sampler_kw
        kw.update( { 'args' : ( gxy_data, gxy_model, gxy_params ) } ) 
        sampler = Sampler( loglikelihood = logprob,
                           ndim = len( gxy_params.par_free ),
                           sampler = hyperpar.sampler,
                           nwalkers = hyperpar.nwalkers,
                           emcee_sampler_kw = kw )

        # Define initial position for walkers
        pos_init = gxy_params.rng.uniform(*gxy_params.par_prior.T,
                                          size=(hyperpar.nwalkers,
                                                len(gxy_params.par_free)))

        # Run sampling
        sampler.run_sampling( pos_init, hyperpar.nsamples,
                              emcee_sampling_kw = hyperpar.sampling_kw )

    # Store results:
    sampler.save_results( out_dir = hyperpar.output_directory,
                          name = hyperpar.run_id,
                          pickle_sampler = hyperpar.pickle_sampler )

    return;

################################################################################
        
def sample_parallel ( hyperpar, Ncpu = None ) :
    import multiprocessing as mp

    # Use all available CPUs if no specific number is provided.
    if Ncpu is None :
        Ncpu = mp.cpu_count()
        
    with mp.Pool( Ncpu ) as pool :

        # It's run-time!
        sampler = None
        if hyperpar.sampler == 'dynesty' :
        
            # Define prior-transform to map parameters in the unit-cube
            from galapy.sampling.Statistics import transform_to_prior_unit_cube

            # Build sampler
            kw = hyperpar.sampler_kw
            kw.update( { 'ptform_args' : ( gxy_params.par_prior, ),
                         'queue_size' : Ncpu } ) 
            sampler = Sampler( loglikelihood = loglikelihood_globals,
                               ndim = len( gxy_params.par_free ),
                               sampler = hyperpar.sampler,
                               prior_transform = transform_to_prior_unit_cube,
                               pool = pool,
                               dynesty_sampler_kw = kw )
            
            # Run sampling
            sampler.run_sampling( dynesty_sampling_kw = hyperpar.sampling_kw )
        
        if hyperpar.sampler == 'emcee' :

            # Build sampler
            kw = hyperpar.sampler_kw
            sampler = Sampler( loglikelihood = logprob_globals,
                               ndim = len( gxy_params.par_free ),
                               sampler = hyperpar.sampler,
                               nwalkers = hyperpar.nwalkers,
                               pool = pool,
                               emcee_sampler_kw = kw )
            
            # Define initial position for walkers
            pos_init = gxy_params.rng.uniform(*gxy_params.par_prior.T,
                                              size=(hyperpar.nwalkers,
                                                    len(gxy_params.par_free)))
            
            # Run sampling
            sampler.run_sampling( pos_init, hyperpar.nsamples,
                                  emcee_sampling_kw = hyperpar.sampling_kw )

    # Store results:
    sampler.save_results( out_dir = hyperpar.output_directory,
                          name = hyperpar.run_id,
                          pickle_sampler = hyperpar.pickle_sampler )
    
    return;

################################################################################

def run () :

    ####################################################################
    # Read command-line arguments:
    
    parser = argparse.ArgumentParser( description = 'options' )
    parser.add_argument( 'parameter_file',
                         default = "",
                         help = 'Path to parameter file')
    parser.add_argument( '-s', '--serial',
                         dest = 'serial',
                         action = 'store_true',
                         help = 'Run the program serially.')
    parser.add_argument( '-mp', '--multiprocessing',
                         dest = 'Ncpu',
                         type = int,
                         default = None,
                         help = ( 'flag for shared-memory parallel run: ' +
                                  'provide number of processes. ' +
                                  'If none is passed defaults to all the available CPUs.') )
    args = parser.parse_args()
    
    spec = importlib.util.spec_from_file_location( "hyper_parameters",
                                                   args.parameter_file )
    hyperpar = importlib.util.module_from_spec( spec )
    spec.loader.exec_module( hyperpar )
    
    ####################################################################

    if args.serial :
        sample_serial( hyperpar )
    else :
        global gxy_data
        global gxy_model
        global gxy_params
        
        init_args = ( hyperpar.bands,
                      hyperpar.fluxes,
                      hyperpar.errors,
                      hyperpar.uplims,
                      hyperpar.filters,
                      hyperpar.parameters,
                      hyperpar.ssp_lib,
                      hyperpar.do_Radio,
                      hyperpar.do_Xray,
                      hyperpar.do_AGN,
                      { 'lstep' : hyperpar.lstep } )
        gxy_data, gxy_model, gxy_params = initialize( *init_args )
        sample_parallel( hyperpar, Ncpu = args.Ncpu )

    return;    

################################################################################

if __name__ == '__main__' :    

    run()

################################################################################
