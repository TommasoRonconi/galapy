""" Implements the class used for storing sampling results.
"""

#############################################################################################
# External imports

import warnings
import numpy
import pickle

#############################################################################################
# Internal imports

from galapy.Galaxy import GXY
from galapy.Noise import Noise
from galapy.Handlers import ModelParameters, GXYParameters, NoiseParameters
from galapy.sampling.Observation import Observation
from galapy.sampling.Sampler import Sampler
from galapy.internal.utils import now_string, func_scalar_or_array, quantile_weighted

#############################################################################################

def generate_output_base ( out_dir = '', name = '' ) :
    """
        
    Parameters
    ----------
    out_dir : string
        Position in the filesystem where the results will be stored.
        Default to the directory where the command has been called.
    name : string 
        A string identifying the run that will be saved. 
        By default it will use a string with the current date+time 
    """
    import os
    
    # If no output directory is passed, set it to the current working directory
    if len( out_dir ) == 0 or out_dir is None :
        out_dir = os.getcwd()
        
    # First check whether the required output directory exists,
    # if not it will be created in the correct position of the file-system
    if not os.path.isdir( out_dir ) :
        try :
            # creates multi-level subdirs. similarly to the *Nix command `mkdir -p`
            # (while os.mkdir() only allows to create the highest level directory)
            os.makedirs( out_dir ) 
        except OSError:
            print ( f"Creation of the directory {out_dir} failed" )
            
    # If no name for the current run has been provided, set it to
    # a string with the current date+time: 'year+month+day+hour+minute'
    if len(name) == 0 :
        name = now_string()

    # Set the string with the output base name
    outbase = os.path.join( out_dir, name )

    return outbase

def dump_results ( model, handler, data, sampler, noise = None, outbase = '' ) :
    from time import time

    if len(outbase) == 0 :
        outbase = nowstring()

    if not isinstance( sampler, Sampler ) :
        raise ValueError( "Argument ``sampler`` should be an instance of type Sampler" )
    sample_res, sample_logl, sample_weights = sampler.return_samples_logl_weights()

    if not isinstance( model, GXY ) :
        raise ValueError( "Argument ``model`` should be an instance of type GXY" )

    if not isinstance( handler, ModelParameters ) :
        raise ValueError( "Argument ``handler`` should be an instance of type ModelParameters" )

    if not isinstance( data, Observation ) :
        raise ValueError( "Argument ``data`` should be an instance of type Observation" )

    if noise is not None :
        if not isinstance( noise, Noise ) :
            raise ValueError( "Argument ``noise`` should be an instance of type Noise" )
        
    print( 'Now processing the sampling results, this might require some time ...' )
    tstart = time()
    results = Results( model, handler,
                       sample_res, sample_logl,
                       sample_weights = sample_weights,
                       data = data, noise = noise,
                       sampler_name = sampler.which_sampler )
    ndur = time() - tstart
    print( f'... done in {ndur} seconds.' )
    
    # Pickle the Results instance
    with open( '_'.join( [ outbase, sampler.which_sampler,
                           'results.galapy.pickle' ] ), 'wb' ) as f :
        pickle.dump( results, f )

    print( f'Results stored in files with prefix: {outbase}' )

    return outbase 

#############################################################################################

class Results () :
    def __init__ ( self, model, handler, sample_res, sample_logl,
                   sample_weights = None, data = None, noise = None,
                   sampler_name = 'dynesty' ) :
        """ A class for storing the results of a sampling run.
        
        Parameters
        ----------
        model : galapy.Galaxy.GXY
            An instance of type GXY (or derived). It stores the model's architecture
            used for running sampling.
        handler : galapy.Handlers.ModelParameters 
            An instance ot type ModelParameters with the parameterisation used in
            the sampling run.
        sample_res : numpy.ndarray
            the matrix containing all the samples of the run
        sample_logl : numpy.ndarray
            1D array with the loglikelihood values corresponding to the samples
        sample_weights : numpy.ndarray
            (Optional) 1D array with the weight of each sample in the run. 
            Default is ``None``, in which case the array will be padded with ones
            (i.e. all samples have the same weight)
        data : galapy.sampling.Observation.Observation
            (Optional) An instance of type ``Observation`` with the fluxes 
            measurements used in the sampling run
        noise : galapy.Noise.Noise
            (Optional) An instance of type ``Noise`` with
            the eventual noise model used in the sampling run
        sampler_name : str
            Which sampler has been used in the sampling run (i.e. 'dynesty' or 'emcee', 
            default is 'dynesty')
        """
        
        # Store the model architecture
        if not isinstance(model, GXY) :
            raise AttributeError( 
                'Attribute "model" should be an instance of type ``GXY``'
            )
        self._mod = pickle.dumps(model)
        
        # Store the parameters specs
        if not isinstance(handler, ModelParameters) :
            raise AttributeError( 
                'Attribute "handler" should be an instance of type ``ModelParameters``'
            )
        self._han = pickle.dumps(handler)

        if sample_weights is None :
            sample_weights = numpy.ones_like( sample_logl )
        
        # Store the observation
        self.Ndof = 1
        self._obs = None
        if data is not None :
            if not isinstance(data,Observation) :
                raise AttributeError( 
                    'Attribute "data" should be an instance of type ``Observation``'
                )
            self._obs = pickle.dumps(data)
            self.Ndof = len(data.pms) - len(handler.par_free)

        # Store the noise model    
        self._noise = None
        if noise is not None :
            if not isinstance( noise, Noise ) :
                raise AttributeError( 
                    'Attribute "noise" should be an instance of type ``Noise``'
                )
            self._noise = pickle.dumps(noise)

        # Store the sampler's name and specs
        self.sampler = sampler_name
        self.ndim = len( handler.par_free )
        
        self.size = len(sample_res)
        if self.size != len(sample_logl) or self.size != len(sample_weights):
            raise RuntimeError( 
                'Arguments sample_res, sample_logl and sample_weights should have same lenght'
            )
        self.params = []
        self.logl    = sample_logl
        self.samples = sample_res
        self.weights = sample_weights
        self.wnot0 = ( self.weights > 0. )
        self.SED    = numpy.empty(shape=(self.size, 
                                         *model.wl().shape))
        self.Mstar  = numpy.empty(shape=(self.size,))
        self.Mdust  = numpy.empty(shape=(self.size,))
        self.Mgas   = numpy.empty(shape=(self.size,))
        self.Zstar  = numpy.empty(shape=(self.size,))
        self.Zgas   = numpy.empty(shape=(self.size,))
        self.SFR    = numpy.empty(shape=(self.size,))
        self.TMC    = numpy.empty(shape=(self.size,))
        self.TDD    = numpy.empty(shape=(self.size,))
        
        for i, par in enumerate(sample_res) :
            self.params += [handler.return_nested(par)['galaxy']]
            model.set_parameters( **self.params[-1] )
            age = model.age
            self.SED[i]   = model.get_SED()
            self.Mstar[i] = model.sfh.Mstar(age)
            self.Mdust[i] = model.sfh.Mdust(age)
            self.Mgas[i]  = model.sfh.Mgas(age)
            self.Zstar[i] = model.sfh.Zstar(age)
            self.Zgas[i]  = model.sfh.Zgas(age)
            self.SFR[i]   = model.sfh(age)
            self.TMC[i]   = model.ism.mc.T
            self.TDD[i]   = model.ism.dd.T
            
    def get_stored_quantities ( self ) :
        """ Returns a list with all the quantities stored in the instance
        """
        return list( self.__dict__.keys() )

    def get_residuals ( self, which_model = 'bestfit', standardised = True ) :
        """
        """
        _gxy = pickle.loads( self._mod )
        _obs = pickle.loads( self._obs )
        _han = pickle.loads( self._han )
        if which_model == 'bestfit' :
            params = self.get_bestfit( 'samples' )
        elif which_model == 'mean' :
            params = self.get_mean( 'samples' )
        elif which_model == 'median' :
            params = self.get_quantile( 'samples' )
        else :
            warnings.warn( 'Choice invalid, falling back to default choice (="bestfit")' )
            params = self.get_bestfit( 'samples' )
            
        nested = _han.return_nested(params)
        _gxy.set_parameters( **nested['galaxy'] )
        # if isinstance( _han, GXYParameters ) :
        #     _gxy.set_parameters( **nested )
        # else :
        #     _gxy.set_parameters( **nested['galaxy'] )
        
        model = _gxy.photoSED()
        data  = _obs.fluxes
        if standardised :
            error = _obs.errors
            if hasattr(self, '_noise') and self._noise is not None :
                _noi = pickle.loads( self._noise )
                _noi.set_parameters( **nested['noise'] )
                error = _noi.apply( error, model )
            return ( data - model ) / error
        return data - model

    def get_chi2 ( self, which_model = 'bestfit', reduced = True ) :
        """
        """
        chi = self.get_residuals( which_model, standardised = True )
        redfact = 1.
        if reduced :
            redfact /= self.Ndof
        return numpy.sum( chi**2 ) * redfact
    
    def get_bestfit ( self, key ) :
        """ Returns the bestfit value of the stored quantity corresponding 
        to the input key
        
        Parameters
        ----------
        key : str or str-sequence 
            if a string, it should name one of the stored quantities 
            if a list of strings, all the strings in the list will be
            matched.  
        
        Returns
        -------
        : scalar or ndarray
            depending on the input ``key``.
        """
        idmax = self.logl.argmax()
        return func_scalar_or_array(
            var = key,
            function = lambda k : self.__dict__[k][idmax]
        )

    def get_mean ( self, key ) :
        """ Returns the weighted mean value of the stored quantity
        corresponding to the input key
        
        Parameters
        ----------
        key : str or str-sequence 
            if a string, it should name one of the stored quantities 
            if a list of strings, all the strings in the list will be
            matched.  
        
        Returns
        -------
        : scalar or ndarray
            depending on the input ``key``.
        """
        return func_scalar_or_array(
            var = key,
            function = lambda k : numpy.average( self.__dict__[k][self.wnot0],
                                                 weights = self.weights[self.wnot0],
                                                 axis = 0 )
        )

    def get_std ( self, key ) :
        """ Returns the weighted standard deviation value of the stored 
        quantity corresponding to the input key
        
        Parameters
        ----------
        key : str or str-sequence 
            if a string, it should name one of the stored quantities 
            if a list of strings, all the strings in the list will be
            matched.  
        
        Returns
        -------
        : scalar or ndarray
            depending on the input ``key``.
        """
        return func_scalar_or_array(
            var = key,
            function = lambda k : numpy.sqrt(
                numpy.average( ( self.__dict__[k][self.wnot0] - self.get_mean( k ) )**2,
                               weights = self.weights[self.wnot0], axis = 0 )
            )
        )

    def get_quantile ( self, key, quantile = 0.5 ) :
        """ Returns the weighted standard deviation value of the stored 
        quantity corresponding to the input key
        
        Parameters
        ----------
        key : str or str-sequence 
            if a string, it should name one of the stored quantities 
            if a list of strings, all the strings in the list will be
            matched.  
        
        Returns
        -------
        : scalar or ndarray
            depending on the input ``key``.
        """
        return func_scalar_or_array(
            var = key,
            function = lambda k : quantile_weighted( self.__dict__[ k ][self.wnot0],
                                                     quantile,
                                                     weights = self.weights[self.wnot0],
                                                     axis = 0 )
        )
    
    def get_model ( self ) :
        """ Instantiates a model corresponding to the one used for sampling.
        
        Returns
        -------
        : galapy.Galaxy.GXY
        """
        return pickle.loads( self._mod )
    
    def get_sampling_params ( self ) :
        """ Instantiates a handler corresponding to the one used for sampling.
        
        Returns
        -------
        : galapy.GalaxyParameters.GXYParameters
        """
        return pickle.loads( self._han )
    
    def get_observation ( self ) :
        """ Instantiates a dataset corresponding to the one used for sampling.
        
        Returns
        -------
        : galapy.sampling.Observation.Observation
        """
        if self._obs is not None :
            return pickle.loads( self._obs )
        else :
            warnings.warn( "The current instance has no Observation stored, passing None" )
            return None
    
    def get_noise ( self ) :
        """ Instantiates a noise model corresponding to the one used for sampling.
        
        Returns
        -------
        : galapy.sampling.Noise.Noise
        """
        if hasattr(self, '_noise') and self._noise is not None :
            return pickle.loads( self._noise )
        else :
            warnings.warn( "The current instance has no Noise stored, passing None" )
            return None

#############################################################################################
