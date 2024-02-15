""" Implements the class used for storing sampling results.
"""

#############################################################################################
# External imports

import warnings
import numpy
import pickle
from collections.abc import MutableMapping as MM

#############################################################################################
# Internal imports

import galapy
from galapy.Galaxy import GXY, PhotoGXY
from galapy.Noise import Noise, CalibrationError
from galapy.Handlers import ModelParameters, GXYParameters, NoiseParameters
from galapy.sampling.Observation import Observation
from galapy.sampling.Sampler import Sampler
from galapy.internal.utils import now_string, func_scalar_or_array, quantile_weighted, get_credible_interval, find_nearest
from galapy.io.hdf5 import write_to_hdf5, load_from_hdf5

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

def dump_results ( model, handler, data, sampler,
                   noise = None, outbase = '',
                   method = 'hdf5', lightweight = False ) :
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

    if lightweight :
        if method in { 'hdf5', 'h5' } :
            outfile = '_'.join( [ outbase, sampler.which_sampler, 'results_light.galapy.hdf5' ] )
            write_to_hdf5(
                outfile,
                metadata = dict(
                    storage_method = 'light',
                    galapy_version = galapy.__version__,
                ),
                hard = True,
                results = {
                    'model' : model.dump(),
                    'handler' : handler.dump(),
                    'sample_res' : sample_res,
                    'sample_logl' : sample_logl,
                    'sample_weights' : sample_weights,
                    'data' : None if data is None else data.dump(),
                    'noise' : None if noise is None else noise.dump(),
                    'sampler_name' : sampler.which_sampler,
                }
            )
            return outbase
        else :
            warnings.warn(
                "Lightweight (`lightweight=True`) dumping only available with `method = 'hdf5'`. "
                "Falling back to HDF5 output."
            )
                
    print( 'Now processing the sampling results, this might require some time ...' )
    tstart = time()
    results = Results( model, handler,
                       sample_res, sample_logl,
                       sample_weights = sample_weights,
                       data = data, noise = noise,
                       sampler_name = sampler.which_sampler )
    ndur = time() - tstart
    print( f'... done in {ndur} seconds.' )

    if method == 'pickle' :
        # Pickle the Results instance
        with open( '_'.join( [ outbase, sampler.which_sampler,
                               'results.galapy.pickle' ] ), 'wb' ) as pfw :
            pickle.dump( results, pfw )
    if method in { 'hdf5', 'h5' } :
        outfile = '_'.join( [ outbase, sampler.which_sampler, 'results.galapy.hdf5' ] )
        write_to_hdf5(
            outfile,
            metadata = dict(
                storage_method = 'heavy',
                galapy_version = galapy.__version__,
            ),
            hard = True,
            results = results.dump()
        )

    print( f'Results stored in files with prefix: {outbase}' )

    return outbase

def load_results ( infile, method = None, lightweight = None ) :

    if method is None :
        method = infile.split('.')[-1]

    res = None
    
    if method in { 'h5', 'hdf5' } :
        res_dict = load_from_hdf5( infile )
        if lightweight is None :
            lightweight = res_dict['metadata']['storage_method'] == 'light'
        if lightweight :
            from time import time
            print( 'Now processing the sampling results, this might require some time ...' )
            tstart = time()
            res = Results(
                model = (
                    PhotoGXY.load(res_dict['results']['model'])
                    if 'pms_kwargs' in res_dict['results']['model']
                    else GXY.load(res_dict['results']['model'])
                ),
                handler = ModelParameters.load( res_dict['results']['handler'] ),
                sample_res = res_dict['results']['sample_res'],
                sample_logl = res_dict['results']['sample_logl'],
                sample_weights = (
                    res_dict['results']['sample_weights']
                    if 'sample_weights' in res_dict['results']
                    else None
                ),
                data = (
                    Observation.load( res_dict['results']['data'] )
                    if 'data' in res_dict['results']
                    else None
                ),
                noise = (
                    CalibrationError.load( res_dict['results']['noise'] )
                    if res_dict['results']['noise'] is not None
                    else None
                ),
                sampler_name = res_dict['results']['sampler_name']
            )
            ndur = time() - tstart
            print( f'... done in {ndur} seconds.' )
        else :
            res = Results.load( res_dict['results'] )
        return res

    if method == 'pickle' :
        with open( infile, 'rb' ) as pfr :
            res = pickle.load( pfr )
        return res

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
        if not isinstance(model, GXY) and not isinstance( model, MM ) :
            raise AttributeError( 
                'Attribute "model" should be an instance of type ``GXY``'
            )
        if isinstance( model, MM ) :
            self._mod = model
            model = self.get_model()
        else :
            self._mod = model.dump()
        
        # Store the parameters specs
        if not isinstance(handler, ModelParameters) and not isinstance( handler, MM ) :
            raise AttributeError( 
                'Attribute "handler" should be an instance of type ``ModelParameters``'
            )
        if isinstance( handler, MM ) :
            self._han = handler
            handler = self.get_handler()
        else :
            self._han = handler.dump()

        if sample_weights is None :
            sample_weights = numpy.ones_like( sample_logl )
        
        # Store the observation
        self.Ndof = 1
        self._obs = None
        if data is not None :
            if not isinstance(data,Observation) and not isinstance( data, MM ) :
                raise AttributeError( 
                    'Attribute "data" should be an instance of type ``Observation``'
                )
            if isinstance( data, MM ) :
                self._obs = data
                data = self.get_observation()
            else :
                self._obs = data.dump()
            self.Ndof = len(data.pms) - len(handler.par_free)

        # Store the noise model    
        self._noise = None
        if noise is not None :
            if not isinstance( noise, Noise ) and not isinstance( noise, MM ) :
                raise AttributeError( 
                    'Attribute "noise" should be an instance of type ``Noise``'
                )
            if isinstance( noise, MM ) :
                self._noise = noise
                noise = self.get_noise()
            else :
                self._noise = noise.dump()

        # Store the sampler's name and specs
        self.sampler = sampler_name
        self.ndim = len( handler.par_free )
        
        self.size = len(sample_res)
        if self.size != len(sample_logl) or self.size != len(sample_weights):
            raise RuntimeError( 
                'Arguments sample_res, sample_logl and sample_weights should have same length'
            )
        self.params = []
        self.logl    = numpy.asarray( sample_logl )
        self.samples = numpy.asarray( sample_res )
        self.weights = numpy.asarray( sample_weights )
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

            try :
                model.set_parameters( **self.params[-1] )
            except RuntimeError :
                self.SED[i]   = -numpy.inf * numpy.ones_like(model.wl())
                self.Mstar[i] = -numpy.inf
                self.Mdust[i] = -numpy.inf
                self.Mgas[i]  = -numpy.inf
                self.Zstar[i] = -numpy.inf
                self.Zgas[i]  = -numpy.inf
                self.SFR[i]   = -numpy.inf
                self.TMC[i]   = -numpy.inf
                self.TDD[i]   = -numpy.inf
                continue
            
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

    def dump ( self ) :
        return dict(
            # Models' architecture
            model = self._mod,
            handler = self._han,
            data = self._obs,
            noise = self._noise,
            # Sampling run hyperparameters
            sampler_name = self.sampler,
            size = self.size, Ndof = self.Ndof,
            # Sampling run stored quantities
            logl    = self.logl,
            samples = self.samples,
            weights = self.weights,
            wnot0   = self.wnot0,
            # Derived quantities
            SED     = self.SED,
            Mstar   = self.Mstar,
            Mdust   = self.Mdust,
            Mgas    = self.Mgas,
            Zstar   = self.Zstar,
            Zgas    = self.Zgas,
            SFR     = self.SFR,
            TMC     = self.TMC,
            TDD     = self.TDD,
        )

    @classmethod
    def load ( cls, dictionary ) :

        # build object
        ret = cls(
            model   = dict( dictionary['model'] ),
            handler = dict( dictionary['handler'] ),
            data    = dict( dictionary['data'] ),
            noise   = ( dict( dictionary['noise'] )
                        if dictionary['noise'] is not None
                        else None ),
            sampler_name = dictionary['sampler_name'],
            sample_res     = [],
            sample_logl    = [],
            sample_weights = []
        )
        # Sampling run stored quantities
        ret.logl    = dictionary['logl']
        ret.samples = dictionary['samples']
        ret.weights = dictionary['weights']
        ret.wnot0   = dictionary['wnot0']

        # Additional hyperparameters
        ret.Ndof    = dictionary['Ndof']
        ret.size    = dictionary['size']

        # Derived quantities
        ret.SED     = dictionary['SED']
        ret.Mstar   = dictionary['Mstar']
        ret.Mdust   = dictionary['Mdust']
        ret.Mgas    = dictionary['Mgas']
        ret.Zstar   = dictionary['Zstar']
        ret.Zgas    = dictionary['Zgas']
        ret.SFR     = dictionary['SFR']
        ret.TMC     = dictionary['TMC']
        ret.TDD     = dictionary['TDD']

        # Compute parameters' dictionaries
        handler = ret.get_handler()
        for i, par in enumerate(ret.samples) :
            ret.params += [handler.return_nested(par)['galaxy']]

        return ret
        
    def get_stored_quantities ( self ) :
        """ Returns a list with all the quantities stored in the instance
        """
        return list( self.__dict__.keys() )

    def get_residuals ( self, which_model = 'bestfit', standardised = True ) :
        """
        """
        _gxy = self.get_model() 
        _obs = self.get_observation()
        _han = self.get_handler()
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
        
        model = _gxy.photoSED()
        data  = _obs.fluxes
        if standardised :
            error = _obs.errors
            if hasattr(self, '_noise') and self._noise is not None :
                _noi = self.get_noise()
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

    def get_median ( self, key ) :
        """Returns the median of a stored quantity.
        This is a shortcut for ``Results.get_quantile( key, quantile=0.5 )``
        """
        return self.get_quantile( key, quantile=0.5 )

    def get_credible_interval ( self, key, percent = 0.68, centre = 'bestfit' ) :
        """Returns the credible interval around some position enclosing a user-defined
        probability integral.
        Automatically accounts for upper or lower limits on some parameter.
        If an interval can be defined it returns the lower and upper distances from 
        the centre of the interval.
        If only upper/lower limits can be defined, it returns the position of the given limit.

        Parameters
        ----------
        key : str
            the name of one of the stored quantities.
        percent : float
            (Optional, default = 0.68) probability enclosed by the interval
        centre : str
            (Optional, default = 'bestfit') one among ('bestfit', 'mean', 'median'),
            where to define the centre of the interval.
            Note that, if this is set to 'median', and percent=0.5, this is equivalent
            to call ``Results.get_quantile`` with the argument ``quantile=(0.25, 0.75)``
        
        Returns
        -------
        low, upp : tuple
            If the interval is completely defined by the samples, returns the credible 
            limits around the centre so that the integral of the posterior between
            ``centre-low`` and ``centre+upp`` is equal to ``percent``.
            In the case only upper/lower limits can be defined it returns either
            ``(-numpy.inf, upp)`` for an upper limit, so that the integral between ``-inf``
            and ``upp`` is equal to ``percent``, or ``(low, +numpy.inf)`` for a lower limit,
            so that the integral between ``-inf`` and ``low`` is equal to ``1-percent``.
        """
        
        samples = self.__dict__[key][self.wnot0]

        if centre == 'bestfit' :
            idcentre = self.logl[self.wnot0].argmax()
        elif centre == 'mean' :
            idcentre = find_nearest( samples, self.get_mean( key ) )
        elif centre == 'median' :
            idcentre = find_nearest( samples, self.get_median( key ) )
        else :
            raise RuntimeError(
                "``centre`` argument should be one among ('bestfit', 'mean', 'median')"
            )

        # get limits from the weighted distribution
        low, upp = get_credible_interval(
            samples, idcentre, percent, self.weights[self.wnot0]
        )

        # if upper-limit is found, return upper-limited indefinite interval
        if low is None :
            return -numpy.inf, upp
        # if lower-limit is found, return lower-limited indefinite interval
        if upp is None :
            return low, +numpy.inf
        # return lower and upper credibility distances from centre
        return samples[idcentre]-low, upp-samples[idcentre]
    
    def get_model ( self ) :
        """ Instantiates a model corresponding to the one used for sampling.
        
        Returns
        -------
        : galapy.Galaxy.GXY
        """
        if 'pms_kwargs' in self._mod.keys() :
            return PhotoGXY.load( self._mod )
        else :
            return GXY.load( self._mod )

    def get_handler ( self ) :
        """ Instantiates a handler corresponding to the one used for sampling.
        
        Returns
        -------
        : galapy.GalaxyParameters.ModelParameters
        """
        return ModelParameters.load( self._han )
    
    def get_observation ( self ) :
        """ Instantiates a dataset corresponding to the one used for sampling.
        
        Returns
        -------
        : galapy.sampling.Observation.Observation
        """
        if self._obs is not None :
            return Observation.load( self._obs )
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
            return CalibrationError.load( self._noise )
        else :
            warnings.warn( "The current instance has no Noise stored, passing None" )
            return None
    
    def get_sampling_params ( self ) :
        """ Instantiates a handler corresponding to the one used for sampling.
        
        Returns
        -------
        : galapy.GalaxyParameters.GXYParameters
        """
        return self.get_handler()

#############################################################################################
