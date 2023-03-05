import numpy
import os
import warnings
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import importlib.util

from galapy.Galaxy import PhotoGXY
from galapy.Noise import CalibrationError
from galapy.Handlers import ModelParameters
from galapy.sampling.Statistics import gaussian_loglikelihood
from galapy.sampling.Sampler import Sampler
from galapy.sampling.Observation import Observation
from galapy.sampling.Results import generate_output_base, dump_results
from galapy.internal.utils import set_nested

################################################################################

def initialize ( bands, fluxes, errors, uplims, filters, params,
                 sfh_model = 'insitu', ssp_lib = 'br22.NT',
                 do_Radio = False, do_Xray = False, do_AGN = False,
                 noise_model = None, noise_params = {},
                 gxy_kwargs = {}, noise_kwargs = {} ) :

    #########################################################################
    # Build observation

    data = Observation( bands, fluxes, errors, uplims, filters )
    
    #########################################################################
    # Build photometric galaxy

    model = PhotoGXY( pms = data.pms,
                      sfh = { 'model' : sfh_model },
                      csp = { 'ssp_lib' : ssp_lib },
                      do_Radio = do_Radio,
                      do_Xray = do_Xray,
                      do_AGN = do_AGN,
                      **gxy_kwargs )
    noise = None
    if noise_model is not None :
        if noise_model == 'calibration_error' :
            noise = CalibrationError( **noise_kwargs )
        else :
            warnings.warn( 'The noise model chosen is not valid, noise will be ignored' )
    
    #########################################################################
    # Build parameters handler

    sample_params = { '.'.join( ['galaxy', key] ) : value for key, value in params.items() }
    if noise is not None :
        sample_params.update( { '.'.join( ['noise', key] ) : value
                                for key, value in noise_params.items() } )
    handler = ModelParameters( model, noise,
                               sample_params = sample_params )

    #########################################################################
    # Set fixed parameters and initial values for free parameters

    init = {}
    for klist, v in handler.parameters.items() :
        set_nested( init, klist.split('.'), v )
    model.set_parameters( **init['galaxy'] )
    if noise is not None : noise.set_parameters( **init['noise'] )

    return data, model, noise, handler

################################################################################

def loglikelihood ( par, data, model, noise, handler, **kwargs ) :

    nested = handler.return_nested( par )
    try :
        model.set_parameters( **nested['galaxy'] )
        flux_model = numpy.asarray( model.photoSED() )
    except RuntimeError :
        return -numpy.inf

    if model.age > model.UA :
        return -numpy.inf

    if noise is not None :
        noise.set_parameters( **nested['noise'] )
        errors = noise.apply( data.errors, flux_model )
        noise_llike = numpy.log( 2 * numpy.pi * errors**2 ).sum()
        if numpy.isnan( noise_llike ) :
            return - numpy.inf
        return gaussian_loglikelihood( data   = data.fluxes,
                                       error  = errors,
                                       model  = flux_model,
                                       uplims = data.uplims,
                                       **kwargs ) - 0.5 * noise_llike
    return gaussian_loglikelihood( data   = data.fluxes,
                                   error  = data.errors,
                                   model  = flux_model,
                                   uplims = data.uplims,
                                   **kwargs )

def loglikelihood_globals ( par, **kwargs ) :
    return loglikelihood( par, *global_args, **kwargs )

################################################################################

# Include prior in likelihood (here uniform prior)
def logprob ( par, data, model, noise, handler, **kwargs ) :
    
    pmin, pmax = handler.par_prior.T
    if all( ( pmin < par ) & ( par < pmax ) ) :
        return loglikelihood( par, data, model, noise, handler, **kwargs )
    
    return -numpy.inf

def logprob_globals ( par, **kwargs ) :
    return logprob( par, *global_args, **kwargs )

################################################################################
        
def sample_serial ( hyperpar ) :

    # Run the initializer, it populates the scope with global variables:
    gxy_data, gxy_model, gxy_noise, gxy_params = initialize(
        bands        = hyperpar.bands,
        fluxes       = hyperpar.fluxes,
        errors       = hyperpar.errors,
        uplims       = hyperpar.uplims,
        filters      = hyperpar.filters,
        params       = hyperpar.galaxy_parameters,
        sfh_model    = hyperpar.sfh_model,
        ssp_lib      = hyperpar.ssp_lib,
        do_Radio     = hyperpar.do_Radio,
        do_Xray      = hyperpar.do_Xray,
        do_AGN       = hyperpar.do_AGN,
        noise_model  = hyperpar.noise_model,
        noise_params = hyperpar.noise_parameters,
        gxy_kwargs = {
            'lstep' : hyperpar.lstep
        },
        noise_kwargs = hyperpar.noise_kwargs
    )

    # It's run-time!
    sampler = None
    if hyperpar.sampler == 'dynesty' :
        
        # Define prior-transform to map parameters in the unit-cube
        from galapy.sampling.Statistics import transform_to_prior_unit_cube

        # Build sampler
        kw = hyperpar.sampler_kw
        kw.update(
            { 'logl_args' : ( gxy_data, gxy_model, gxy_noise, gxy_params ),
              'logl_kwargs' : { 'method_uplims' : hyperpar.method_uplims },
              'ptform_args' : (gxy_params.par_prior,)
            }
        ) 
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
        kw.update(
            { 'args' : ( gxy_data, gxy_model, gxy_noise, gxy_params ),
              'kwargs' : { 'method_uplims' : hyperpar.method_uplims }
            }
        ) 
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
    outbase = generate_output_base( out_dir = hyperpar.output_directory,
                                    name = hyperpar.run_id )
    _ = dump_results( model = gxy_model, handler = gxy_params, data = gxy_data,
                      sampler = sampler, noise = gxy_noise,outbase = outbase )
    sampler.save_results( outbase = outbase,
                          pickle_sampler = hyperpar.pickle_sampler,
                          pickle_raw = hyperpar.pickle_raw )

    return;

################################################################################
        
def sample_parallel ( hyperpar, Ncpu = None ) :
    import multiprocessing as mp

    
    # Use all available CPUs if no specific number is provided.
    if Ncpu is None :
        Ncpu = mp.cpu_count()
        
    with mp.Pool( Ncpu ) as pool :

        gxy_data, gxy_model, gxy_noise, gxy_params = global_args
        # It's run-time!
        sampler = None
        if hyperpar.sampler == 'dynesty' :
        
            # Define prior-transform to map parameters in the unit-cube
            from galapy.sampling.Statistics import transform_to_prior_unit_cube

            # Build sampler
            kw = hyperpar.sampler_kw
            kw.update(
                { 'ptform_args' : ( gxy_params.par_prior, ),
                  'logl_kwargs' : { 'method_uplims' : hyperpar.method_uplims },
                  'queue_size' : Ncpu
                }
            ) 
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
            kw.update( { 'kwargs' : { 'method_uplims' : hyperpar.method_uplims } } )
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
    outbase = generate_output_base( out_dir = hyperpar.output_directory,
                                    name = hyperpar.run_id )
    _ = dump_results( model = gxy_model, handler = gxy_params, data = gxy_data,
                      sampler = sampler, noise = gxy_noise, outbase = outbase )
    sampler.save_results( outbase = outbase,
                          pickle_sampler = hyperpar.pickle_sampler,
                          pickle_raw = hyperpar.pickle_raw )
    
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
        global global_args
        # global gxy_data
        # global gxy_model
        # global gxy_noise
        # global gxy_params
    
        init_args = (
            hyperpar.bands,
            hyperpar.fluxes,
            hyperpar.errors,
            hyperpar.uplims,
            hyperpar.filters,
            hyperpar.galaxy_parameters,
        )
        init_kwargs = dict(
            sfh_model    = hyperpar.sfh_model,
            ssp_lib      = hyperpar.ssp_lib,
            do_Radio     = hyperpar.do_Radio,
            do_Xray      = hyperpar.do_Xray,
            do_AGN       = hyperpar.do_AGN,
            noise_model  = hyperpar.noise_model,
            noise_params = hyperpar.noise_parameters,
            gxy_kwargs = {
                'lstep' : hyperpar.lstep
            },
            noise_kwargs = hyperpar.noise_kwargs,
        )
        # gxy_data, gxy_model, gxy_params = initialize( *init_args, **init_kwargs )
        global_args = initialize( *init_args, **init_kwargs )
        sample_parallel( hyperpar, Ncpu = args.Ncpu )

    return;    

################################################################################
# This commented out here because reasons ...
# if __name__ == '__main__' :    
#     run()
################################################################################

default_parameter_file = f"""
##################################################
# Parameters for building the observation to fit #
##################################################

# The photometric filters system used in the observation.
# This parameter can be either an iterable containing names
# of filters already present in the database, e.g.
#
# filters = ['b_goods', 'i_goods', 'v_goods', 'z_goods']
#
# Alternatively it can be a nested dictionary with user-defined transmissions.
# Such transmissions are passed as properly formatted dictionaries
# 
# filters = {{ 'filter_1' : {{ 'wavelenghts' : array-like,
#                            'photons' : array-like }},
#             'filter_2' : {{ 'wavelenghts' : array-like,
#                            'photons' : array-like }},
#             ...
#             'filter_N' : {{ 'wavelenghts' : array-like,
#                            'photons' : array-like }} }}
#
# As the keywords in the lower level dictionaries suggest,
# the two arrays provided, must define the
# wavelenght grid and the corresponding transmission in photon units
# (and thus they must have the same size).
# Note that the chosen keyword in the higher level dictionary
# (e.g. 'filter_1') will be used as unique identifier of
# the custom transmission.
filters = None

# The observed dataset is expressed as 4 array-likes containing respectively:
# - bands: a list of strings corresponding to the unique identifiers also used
#          in the 'filters' variable above (i.e. if you have provided as a list
#          of filters already present in the database you can set
#          bands = filters
# - fluxes: measures of fluxes (or upper limits) at the corresponding bands
#           listed in variable 'bands'.
#           The input measure should be given in units of milli-Jansky
# - errors: 1-sigma error on the measures of the fluxes (or upper limits) listed
#           in the variable 'fluxes'
# - uplims: sequence of booleans identifying whether the values listed 
#           in argument `fluxes` should be considered non-detection (`True`)
#           or a detection (`False`)
bands  = None
fluxes = None
errors = None
uplims = None

# Method for treatment of the upper limits, when present.
# Available methods are:
# - 'simple' : a heaviside function which is 0. when the model's flux is
#              smaller than the upper-limit and +infty otherwise
# - 'chi2' : the distance between model and data is expressed as a normal chi-squared
# - 'S12' : modified version of the chi-squared integrating a gaussian, with mean = data-flux
#           and std = data-error, up to the value of the model's flux (Sawicki, 2012)
method_uplims = 'chi2'

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
sfh_model = 'insitu'

# The SSP library to use for building the unattenuated stellar emission component
ssp_lib   = 'br22.NT'

# Whether to provide X-ray emission support (True) or not (False). 
do_Xray = False

# Whether to provide radio-emission support (True) or not (False).
do_Radio = False

# Whether to build a galaxy containing an AGN (True) or not (False).
do_AGN = False

# Sub-sampling of the wavelenght grid.
# If lstep is an integer it will consider a wavelenght grid entry every lstep values.
# If lstep is a sequence of integers or a mask, only the wavelenght grid entries
# corresponding to the indices provided will be considered.
# If None, it will consider the whole wavelenght grid (safest choice)
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
# (e.g. if the galaxy model has been built with `do_AGN = False` all the eventual
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
#   The list `a_list` contains the minimum and maximum value of the
#   UNIFORM prior from which to draw samples for the 'free_parameter'.
#   The boolean `a_bool` states whether the prior has to be considered
#   * logarithmic: `a_bool = True`, therefore samples will be drawn from the interval
#                  10**min(a_list) < 'free_parameter' < 10**max(a_list)
#   * linear: `a_bool = False`, therefore samples will be drawn from the interval
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

    # In-Situ model
    'sfh.psi_max' : ( [0., 4.], True ),
    'sfh.tau_star' : ( [6., 11.], True ),

    # Constant model
    'sfh.psi'   : ( [0., 4.], True ),
    'sfh.Mdust' : ( [6., 14.], True ),
    'sfh.Zgxy'  : ( [0., 1.], False ),
    
    # Delayed-Exponential model
    'sfh.psi_norm' : ( [0., 4.], True ),
    'sfh.k_shape'  : ( [0., 5.], False ),
    'sfh.tau_star' : ( [6., 11.], True ),
    'sfh.Mdust' : ( [6., 14.], True ),
    'sfh.Zgxy'  : ( [0., 1.], False ),
    
    # Log-Normal model
    'sfh.psi_norm'   : ( [0., 4.], True ),
    'sfh.sigma_star' : ( [0., 5.], False ),
    'sfh.tau_star'   : ( [6., 11.], True ),
    'sfh.Mdust' : ( [6., 14.], True ),
    'sfh.Zgxy'  : ( [0., 1.], False ),

    ########################
    # Inter-Stellar Medium #
    ########################

    'ism.f_MC' : ( [0., 1.], False ),

    ##################
    # Molecular Clouds
    
    'ism.norm_MC' : ( [-1., 4.], True ),
    'ism.N_MC'    : ( [0., 5.], True ),
    'ism.R_MC'    : ( [0., 5.], True ),
    'ism.tau_esc' : ( [4., 8.], True ),
    'ism.dMClow'  : ( [0., 5.], False ),
    'ism.dMCupp'  : ( [0., 5.], False ),
    
    ##############
    # Diffuse Dust
    
    'ism.norm_DD' : ( [-1., 4.], True ),
    'ism.Rdust'   : ( [0., 5.], True ),
    'ism.f_PAH'   : ( [0., 1.], False ),
    'ism.dDDlow'  : ( [0., 5.], False ),
    'ism.dDDupp'  : ( [0., 5.], False ),

    ###############
    # Synchrotorn #
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
    
    'agn.ct' : 40,
    'agn.al' : 0.,
    'agn.be' : -0.5,
    'agn.ta' : 6.,
    'agn.rm' : 60,
    'agn.ia' : 0.001,
    
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
sampler = 'emcee'

# Sampler keyword arguments.
# These are the parameters that will be passed to the constructor of
# the chosen sampler. See relative documentation for further informations:
# - emcee :
# - dynesty :
sampler_kw = {{}}

# Sampling keyword arguments
# These are the parameters that will be passed to the routine running the
# parameter-space sampling. See relative documentation for further informations:
# - emcee :
# - dynesty :
sampling_kw = {{}}

# Output directory (note that if the directory does not exist it will be created)
output_directory = ''

# An identification name, it will be pre-pended to all files stored in the output directory
# (if the string is empty the current date+time will be used)
run_id = ''

# Whether to pickle the sampler raw results.
# (might be useful for analyzing run statistics)
pickle_raw = True

# Whether to pickle the sampler at the end-of-run state.
# (might be useful for extending the run)
pickle_sampler = False

# EMCEE SAMPLER-SPECIFIC MANDATORY PARAMETERS
# - set the number of walkers (`nwalkers`)
# - set the chain lenght (`nsamples`)
nwalkers = 32
nsamples = 16

#############################
# ... and that's all folks! #
#############################
"""

def generate_parameter_file () :

    ####################################################################
    # Read command-line arguments:
    
    parser = argparse.ArgumentParser( description = 'options' )
    parser.add_argument( '--name', '-n',
                         dest = 'name',
                         type = str,
                         default = 'galapy_hyper_parameters',
                         help = (
                             'provide here the name you want to give to ' +
                             'the parameter file, you can also choose ' +
                             'a path different to the current working directory: ' +
                             '/path/to/chosen_name' +
                             'NOTE THAT the extension ".py" will be appended to ' +
                             'the name chosen (this just to guarantee proper formatting ' +
                             'when opening the file in a text editor). ' +
                             'DEFAULT: ${PWD}/galapy_hyper_parameters.py'
                         ) )
    args = parser.parse_args()
    
    ####################################################################

    with open( args.name + '.py', 'w' ) as paramfile :
        paramfile.write( default_parameter_file )

    return;

################################################################################
