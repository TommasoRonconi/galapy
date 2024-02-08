import numpy
import os
import warnings
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import importlib.util

from galapy.PhotometricSystem import PMS
from galapy.Galaxy import PhotoGXY
from galapy.Noise import CalibrationError
from galapy.Handlers import ModelParameters
from galapy.sampling.Statistics import gaussian_loglikelihood
from galapy.sampling.Sampler import Sampler
from galapy.sampling.Observation import Observation
from galapy.sampling.Results import generate_output_base, dump_results
from galapy.internal.utils import set_nested

_default_sampling_kw = {
    'dynesty' : dict(
        dlogz_init=0.05, nlive_init=500, nlive_batch=100,
        maxiter_init=10000, maxiter_batch=1000, maxbatch=10,
        stop_kwargs = {'target_n_effective': int(5e6)}
    ),
    'emcee' : {},
}    

global_dict = None

################################################################################

def initialize ( bands, fluxes, errors, uplims, filter_args, params,
                 sfh_model = 'insitu', ssp_lib = 'parsec22.NT',
                 do_Radio = False, do_Xray = False, do_AGN = False,
                 noise_model = None, noise_params = {},
                 gxy_kwargs = {}, noise_kwargs = {}, filter_kwargs = {} ) :

    global global_dict 
    
    #########################################################################
    # Build photometric system

    pms = PMS( *filter_args, **filter_kwargs )
    
    #########################################################################
    # Build observation

    data = Observation( bands, fluxes, errors, uplims, pms )
    
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
    else :
        handler = ModelParameters( model, sample_params = sample_params )

    #########################################################################
    # Set fixed parameters and initial values for free parameters

    #init = {}
    #for klist, v in handler.parameters.items() :
    #    set_nested( init, klist.split('.'), v )
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

    global_dict = { 'data' : data, 'model' : model, 'noise' : noise, 'handler' : handler }
    return;

################################################################################

def loglikelihood ( par, **kwargs ) :
    
    data    = global_dict['data']
    model   = global_dict['model']
    noise   = global_dict['noise']
    handler = global_dict['handler']

    nested = handler.return_nested( par )
    try :
        model.set_parameters( **nested['galaxy'] )
        flux_model = numpy.asarray( model.photoSED() )
    except RuntimeError :
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

################################################################################

# Include prior in likelihood (here uniform prior)
def logprob ( par, **kwargs ) :
    
    pmin, pmax = global_dict['handler'].par_prior.T
    if all( ( pmin < par ) & ( par < pmax ) ) :
        return loglikelihood( par, **kwargs )
    
    return -numpy.inf

################################################################################

def sample ( sampler = 'emcee', nwalkers = None, nsamples = None, 
             sampler_kw = {}, logl_kw = {}, run_sampling_kw = {},
             Ncpu = 1, pool = None ) :
    
    # Sampler keywords
    sampler_kw = dict( sampler_kw )
    
    # Sampling keywords
    sampling_kw = dict(_default_sampling_kw[sampler])
    sampling_kw.update(run_sampling_kw)
    
    # It's run-time!
    if sampler == 'dynesty' :
        
        # Define prior-transform to map parameters in the unit-cube
        from galapy.sampling.Statistics import transform_to_prior_unit_cube
    
        # Build sampler
        sampler_kw.update(
            { 'ptform_args' : ( global_dict['handler'].par_prior, ),
              'logl_kwargs' : logl_kw,
              'queue_size' : Ncpu
            }
        ) 
        sampler = Sampler( loglikelihood = loglikelihood,
                           ndim = len( global_dict['handler'].par_free ),
                           sampler = sampler,
                           prior_transform = transform_to_prior_unit_cube,
                           pool = pool,
                           dynesty_sampler_kw = sampler_kw )

        # Run sampling
        sampler.run_sampling( dynesty_sampling_kw = sampling_kw )
    
    if sampler == 'emcee' :
        
        # Build sampler
        sampler_kw.update( { 'kwargs' : logl_kw } )
        sampler = Sampler( loglikelihood = logprob,
                           ndim = len( global_dict['handler'].par_free ),
                           sampler = sampler,
                           nwalkers = nwalkers,
                           pool = pool,
                           emcee_sampler_kw = sampler_kw )

        # Define initial position for walkers
        pos_init = global_dict['handler'].rng.uniform(*global_dict['handler'].par_prior.T,
                                          size=(nwalkers,
                                                len(global_dict['handler'].par_free)))

        # Run sampling
        sampler.run_sampling( pos_init, nsamples,
                              emcee_sampling_kw = sampling_kw )
    
    # return the sampler for storing
    return sampler

################################################################################

def store_results ( sampler, 
                    out_dir = '.', name = '',
                    method = 'hdf5', lightweight = False,
                    pickle_sampler = False, pickle_raw = False ) :
        
    # Store results:
    outbase = generate_output_base( out_dir = out_dir,
                                    name = name )
    _ = dump_results( model = global_dict['model'],
                      handler = global_dict['handler'],
                      data = global_dict['data'],
                      sampler = sampler,
                      noise = global_dict['noise'],
                      outbase = outbase,
                      method = method,
                      lightweight = lightweight
                     )
    sampler.save_results( outbase = outbase,
                          pickle_sampler = pickle_sampler,
                          pickle_raw = pickle_raw )
    return;

################################################################################

def _sample_serial ( which_sampler = 'emcee', nwalkers = None, nsamples = None, 
                     sampler_kw = {}, logl_kw = {}, run_sampling_kw = {}, 
                     out_dir = '.', name = '',
                     store_method = 'hdf5', store_lightweight = False, 
                     pickle_sampler = False, pickle_raw = False  ) :

    # Run the sampling
    sampler = sample(
        sampler = which_sampler, 
        sampler_kw = sampler_kw, 
        logl_kw = logl_kw,
        run_sampling_kw = run_sampling_kw,
        nwalkers = nwalkers,
        nsamples = nsamples
    )
    
    # Store the results
    store_results(
        sampler, 
        out_dir = out_dir,
        name = name,
        method = store_method,
        lightweight = store_lightweight, 
        pickle_sampler = pickle_sampler,
        pickle_raw = pickle_raw
    )
    
    return;

################################################################################

def _sample_parallel ( which_sampler = 'emcee', nwalkers = None, nsamples = None, 
                       sampler_kw = {}, logl_kw = {}, run_sampling_kw = {},
                       Ncpu = None, 
                       out_dir = '.', name = '', 
                       store_method = 'hdf5', store_lightweight = False,
                       pickle_sampler = False, pickle_raw = False  ) :
    import multiprocessing as mp
    
    # Use all available CPUs if no specific number is provided.
    if Ncpu is None :
        Ncpu = mp.cpu_count()

    # Run the sampling
    sampler = None
    with mp.get_context("fork").Pool( Ncpu ) as pool :
        sampler = sample(
            sampler = which_sampler, 
            sampler_kw = sampler_kw, 
            logl_kw = logl_kw,
            run_sampling_kw = run_sampling_kw,
            nwalkers = nwalkers,
            nsamples = nsamples,
            Ncpu = Ncpu, pool = pool
        )

    # Store the results
    store_results(
        sampler, 
        out_dir = out_dir,
        name = name, 
        method = store_method,
        lightweight = store_lightweight, 
        pickle_sampler = pickle_sampler,
        pickle_raw = pickle_raw
    )
    
    return;

################################################################################

def _run () :

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

    # Initialization
    init_args = (
        hyperpar.bands,
        hyperpar.fluxes,
        hyperpar.errors,
        hyperpar.uplims if hyperpar.uplims is not None
        else numpy.zeros_like(hyperpar.bands, dtype=bool),
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
        filter_kwargs = hyperpar.filters_custom if hyperpar.filters_custom is not None else {},
    )
    initialize( *init_args, **init_kwargs )

    # Run
    if args.serial :
        _sample_serial(
            which_sampler = hyperpar.sampler,
            nwalkers = hyperpar.nwalkers,
            nsamples = hyperpar.nsamples, 
            sampler_kw = hyperpar.sampler_kw, 
            logl_kw = { 'method_uplims' : hyperpar.method_uplims }, 
            run_sampling_kw = hyperpar.sampling_kw,
            out_dir = hyperpar.output_directory,
            name = hyperpar.run_id,
            store_method = hyperpar.store_method,
            store_lightweight = hyperpar.store_lightweight,
            pickle_sampler = hyperpar.pickle_sampler,
            pickle_raw = hyperpar.pickle_raw
        )
    else :
        _sample_parallel(
            which_sampler = hyperpar.sampler,
            nwalkers = hyperpar.nwalkers,
            nsamples = hyperpar.nsamples, 
            sampler_kw = hyperpar.sampler_kw, 
            logl_kw = { 'method_uplims' : hyperpar.method_uplims }, 
            run_sampling_kw = hyperpar.sampling_kw,
            Ncpu = args.Ncpu, 
            out_dir = hyperpar.output_directory,
            name = hyperpar.run_id,
            store_method = hyperpar.store_method,
            store_lightweight = hyperpar.store_lightweight, 
            pickle_sampler = hyperpar.pickle_sampler,
            pickle_raw = hyperpar.pickle_raw
        )

    return;    

################################################################################
# This commented out here because reasons ...
# if __name__ == '__main__' :    
#     run()
################################################################################

default_parameter_file = """
##################################################
# Parameters for building the observation to fit #
##################################################

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
# filters = {{ 'filter_1' : {{ 'wavelengths' : array-like,
#                            'photons' : array-like }},
#             'filter_2' : {{ 'wavelengths' : array-like,
#                            'photons' : array-like }},
#             ...
#             'filter_N' : {{ 'wavelengths' : array-like,
#                            'photons' : array-like }} }}
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

# EMCEE SAMPLER-SPECIFIC MANDATORY PARAMETERS
# - set the number of walkers (``nwalkers``)
# - set the chain length (``nsamples``)
nwalkers = 64
nsamples = 4096

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

    with open( args.name + '.py', 'w' ) as paramfile :
        paramfile.write( default_parameter_file.format(args.sfh_model, sfh_params_string) )

    return;

################################################################################
