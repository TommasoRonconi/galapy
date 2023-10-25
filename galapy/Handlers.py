########################################################################################

# External imports
import numpy
import warnings

# Internal imports
from galapy.Galaxy import GXY
from galapy.Noise import Noise
from galapy.internal.abc import Model
from galapy.internal.utils import set_nested, unwrap_items

########################################################################################

class ModelParameters () :

    def __init__ ( self, *args, sample_params = {}, rng_kw = {} ) :

        # Set the default random engine 
        self.rng = numpy.random.default_rng(**rng_kw)
        
        # Extract parameterisation from the models
        self.parameters = {}
        _i = 0
        for model in args :
            if not isinstance( model, Model ) :
                raise RuntimeError(
                    f'The provided model ``{model}`` is not valid. '
                    'Valid models are instances of class galapy.internal.abc.Model'
                )
            if isinstance( model, GXY ) :
                model_key = ['galaxy']
            elif isinstance( model, Noise ) :
                model_key = ['noise'] 
            else :
                model_key = [f'custom{_i}'] 
                _i += 1
            self.parameters.update(
                { '.'.join(model_key + k) : v for k, v in unwrap_items(model.params) }
            )
            
        # Read fixed and free parameters from dictionary
        self.par_log   = []
        self.par_free  = []
        self.par_prior = []
        for key, value in sample_params.items() :
            if key not in self.parameters.keys() :
                warnings.warn(
                    f'Parameter "{key}" is not present in the model '
                    'requested and will be ignored.')
                continue
            if isinstance(value, tuple) :
                prior, log = value
                self.par_free.append( key )
                self.par_prior.append( prior )
                self.par_log.append( log )
                v = self.rng.uniform( *prior, size=None )
                self.parameters[key] = 10.**v if log else v
            else :
                try :
                    self.parameters[key] = float(value)
                except ValueError as ve :
                    warnings.warn( f'The value ({value}) provided for the fixed parameter {key} '
                                   'is not valid and will be ignored. '
                                   f'Using the default value {self.parameters[key]} instead.' )
        self.par_log   = numpy.asarray(self.par_log)
        self.par_free  = numpy.asarray(self.par_free)
        self.par_prior = numpy.asarray(self.par_prior)
                    
    def return_nested ( self, par = None ) :
        """ From a list of parameters returns a nested dictionary in 
        the proper format. Here 'proper' means that such dictionary has
        the hierarchy necessary to pass it as keyword arguments dictionary
        to an object of type galapy.Galaxy.GXY()
        
        Parameters
        ----------
        par : array-like or iterable
          a list of values to be assigned to the parameters flagged as `free`
        
        Returns
        -------
        : dict
        A nested dictionary with the proper format to be passed as keyword 
        arguments to the constructor or to function `set_parameters` of 
        a `galapy.Galaxy.GXY()` object.
        
        Examples
        --------
        Given a galaxy model of type GXY:

        >>> from galapy.Galaxy import GXY
        >>> gxy = GXY()

        and a dictionary with the fixed and free parameters as follows
        
        >>> p = { 'age' : 1.e+8, 'sfh.psi_max' : ( [0., 5.], True ) }

        we can build a parameter handler:
 
        >>> from galapy.GalaxyParameters import GXYParameters
        >>> par = GXYParameters( gxy, p )
        >>> par.par_free
        numpy.array( [ 'sfh.psi_max' ] )
        
        Now, assuming we want to set a new value to the 'sfh.psi_max' parameter,
        >>> nested = par.return_nested( [ 123. ] )
        >>> nested 
        { 'sfh' : { 'psi_max' : 123. } }
        
        This dictionary can be used to set the parameters of the galaxy model
        >>> gxy.set_parameters( **nested )
        """

        # If no argument is passed use internal parameters
        ret = {}
        if par is None :
            ret = {}
            for klist, v in self.parameters.items() :
                set_nested( ret, klist.split('.'), v )
            return ret
        
        # Otherwise use free parameters
        par = numpy.asarray( par ) 
        if len(par) != len(self.par_free) :
            raise RuntimeError( f'Provided {len(par)} but there are exactly '
                                f'{len(self.par_free)} free parameters to set.' )
        val = numpy.zeros_like(par)
        val[self.par_log] = 10.**par[self.par_log]
        val[~self.par_log] = par[~self.par_log]
        for v, klist in zip( val, self.par_free ) :
            set_nested( ret, klist.split('.'), v )
        return ret


########################################################################################

class GXYParameters ( ModelParameters ) :
    def __init__ ( self, gxy_model, sample_params, rng_kw = {} ) :

        super().__init__( gxy_model,
                          sample_params = { '.'.join(['galaxy',k]) : v
                                            for k,v in sample_params.items() },
                          rng_kw = rng_kw )

    def return_nested ( self, par ) :
        return super().return_nested( par )['galaxy']
    

########################################################################################

class NoiseParameters ( ModelParameters ) :
    def __init__ ( self, noise_model, sample_params, rng_kw = {} ) :

        super().__init__( noise_model,
                          sample_params = { '.'.join(['noise',k]) : v
                                            for k,v in sample_params.items() },
                          rng_kw = rng_kw )

    def return_nested ( self, par ) :
        return super().return_nested( par )['noise']

########################################################################################
