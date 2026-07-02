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
                   method = 'hdf5', lightweight = False,
                   derived = None ) :
    from time import time

    if len(outbase) == 0 :
        outbase = nowstring()

    if not isinstance( sampler, Sampler ) :
        raise ValueError( "Argument ``sampler`` should be an instance of type Sampler" )
    sample_res, sample_logl, sample_weights = sampler.return_samples_logl_weights()
    logz, logzerr = sampler.log_evidence()

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
                    'logz'    : logz,
                    'logzerr' : logzerr,
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
                       sample_weights   = sample_weights,
                       data             = data, noise = noise,
                       sampler_name     = sampler.which_sampler,
                       log_evidence     = logz,
                       log_evidence_err = logzerr,
                       derived          = derived )
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
                sampler_name     = res_dict['results']['sampler_name'],
                log_evidence     = res_dict['results'].get( 'logz',    None ),
                log_evidence_err = res_dict['results'].get( 'logzerr', None ),
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

    # Registry of the derived quantities computed for every sample at
    # construction time. Each entry maps a stored-attribute name to a callable
    # ``f(model)`` returning that quantity for a model whose parameters have
    # already been set to the sample's values. ``model.age`` is a plain
    # attribute (set inside ``set_parameters``), so functions that need the age
    # recover it from the model instead of receiving it as a separate argument.
    _default_properties = {
        'SED'   : lambda model : model.get_SED(),
        'Mstar' : lambda model : model.sfh.Mstar( model.age ),
        'Mdust' : lambda model : model.sfh.Mdust( model.age ),
        'Mgas'  : lambda model : model.sfh.Mgas( model.age ),
        'Zstar' : lambda model : model.sfh.Zstar( model.age ),
        'Zgas'  : lambda model : model.sfh.Zgas( model.age ),
        'SFR'   : lambda model : model.sfh( model.age ),
        'TMC'   : lambda model : model.ism.mc.T,
        'TDD'   : lambda model : model.ism.dd.T,
    }

    def __init__ ( self, model, handler, sample_res, sample_logl,
                   sample_weights = None, data = None, noise = None,
                   sampler_name = 'dynesty',
                   log_evidence = None, log_evidence_err = None,
                   derived = None ) :
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
        log_evidence : float
            (Optional) the logarithm of the evidence accumulated during the run (this
            number is not available for runs performed with the 'emcee' sampler)
        log_evidence_err : float
            (Optional) the error on the logarithm of the evidence accumulated during
            the run (this number is only available for runs performed with the 'dynesty' sampler)
        derived : sequence of str
            (Optional, default = ``None``) the subset of default derived
            quantities (keys of ``_default_properties``) to compute and store.
            ``None`` selects them all. ``'SED'`` is always computed and stored
            regardless of this argument; listing it explicitly is allowed but
            redundant (a warning is emitted).
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

        # Store the sampler's name, specs, and evidence estimate
        self.sampler = sampler_name
        self.ndim    = len( handler.par_free )
        self.logz    = log_evidence
        self.logzerr = log_evidence_err
        
        self.size = len(sample_res)
        if self.size != len(sample_logl) or self.size != len(sample_weights):
            raise RuntimeError( 
                'Arguments sample_res, sample_logl and sample_weights should have same length'
            )
        self.params  = []
        self.logl    = numpy.asarray( sample_logl )
        self.samples = numpy.asarray( sample_res )
        self.weights = numpy.asarray( sample_weights )
        self.wnot0   = ( self.weights > 0. )

        # Build the per-sample list of (free) nested parameter dictionaries.
        # The user-fixed parameters are not part of ``return_nested(par)``;
        # they are applied separately inside ``_compute_properties``.
        for par in sample_res :
            self.params += [ handler.return_nested( par )['galaxy'] ]

        # Compute and store the selected derived quantities. ``_derived`` keeps
        # track of which quantities are stored so that ``dump``/``load`` and
        # ``add_property`` stay in sync. ``SED`` is always computed and stored.
        self._derived = []
        self._compute_properties( self._select_properties( derived ),
                                  model = model, handler = handler )

    def _select_properties ( self, derived = None ) :
        """Build the ordered ``{name: callable}`` mapping of default quantities
        to compute, given an optional user selection.

        ``SED`` is always included (it is the primary product and the reference
        for the physical-validity gate in ``_compute_properties``). If the user
        lists it explicitly it is silently de-duplicated with a warning.

        Parameters
        ----------
        derived : sequence of str, optional
            The subset of ``_default_properties`` keys to store. ``None`` selects
            them all.

        Returns
        -------
        : dict
            Mapping ``{name: callable}`` drawn from ``_default_properties``,
            always containing ``'SED'``.
        """
        if derived is None :
            return dict( self._default_properties )

        names = list( derived )
        if 'SED' in names :
            warnings.warn( "'SED' is always computed and stored by default; "
                           "ignoring its presence in the requested `derived` list." )
            names = [ n for n in names if n != 'SED' ]

        unknown = [ n for n in names if n not in self._default_properties ]
        if len( unknown ) > 0 :
            raise KeyError(
                f"Unknown derived quantities {unknown}; valid choices are "
                f"{list( self._default_properties.keys() )}."
            )

        # 'SED' first, then the requested ones in the user-provided order.
        selected = { 'SED' : self._default_properties['SED'] }
        for n in names :
            selected[n] = self._default_properties[n]
        return selected

    def _compute_properties ( self, funcs, model = None, handler = None ) :
        """Evaluate derived-property functions over all the stored samples.

        For each sample the model parameters are set and every callable in
        ``funcs`` is evaluated as ``f(model)``; the returned values are stored
        as instance attributes named after the dictionary keys. The storage
        arrays are sized from the shape of each returned value, so scalar- and
        array-valued quantities (e.g. ``SED``) are handled uniformly.

        On a ``RuntimeError`` from ``set_parameters`` the sample is filled with
        ``-inf`` for every property. Historical behaviour is preserved: if the
        ``SED`` of a sample is non-finite, all the *other* quantities of that
        sample are set to ``-inf`` while the SED itself keeps its values.

        Parameters
        ----------
        funcs : dict
            Mapping ``{ name : callable }`` where each callable takes the model
            (with the current sample's parameters set) and returns the derived
            quantity for that sample.
        model : galapy.Galaxy.GXY, optional
            A model instance to reuse. If ``None`` a fresh one is built with
            ``get_model``.
        handler : galapy.Handlers.ModelParameters, optional
            The parameter handler to reuse. If ``None`` it is rebuilt with
            ``get_handler``.
        """
        if model is None :
            model = self.get_model()
        if handler is None :
            handler = self.get_handler()

        keys = list( funcs.keys() )
        _sentinel = -numpy.inf

        # Apply the user-fixed parameters once: return_nested() with no
        # argument returns *all* the stored parameters (fixed included), while
        # the per-sample dictionaries in self.params hold only the free ones.
        # Evaluating the functions on this valid state also lets us pre-size
        # the storage arrays, so array-valued quantities keep their shape even
        # if every sample below fails ``set_parameters``.
        out = None
        try :
            model.set_parameters( **handler.return_nested()['galaxy'] )
            with numpy.errstate( all = 'ignore' ) :
                probe = { k : funcs[k]( model ) for k in keys }
            out = { k : numpy.full( ( self.size, *numpy.shape( v ) ),
                                    _sentinel, dtype = float )
                    for k, v in probe.items() }
        except RuntimeError :
            pass   # fixed configuration invalid: fall back to lazy allocation

        for i, nested in enumerate( self.params ) :
            with numpy.errstate( all = 'ignore' ) :
                try :
                    model.set_parameters( **nested )
                except RuntimeError :
                    if out is not None :
                        for k in keys :
                            out[k][i] = _sentinel
                    continue
                vals = { k : funcs[k]( model ) for k in keys }

            # Allocate on the first successful sample if the probe above failed
            # (already-failed samples keep the sentinel value from ``full``).
            if out is None :
                out = { k : numpy.full( ( self.size, *numpy.shape( v ) ),
                                        _sentinel, dtype = float )
                        for k, v in vals.items() }

            for k, v in vals.items() :
                out[k][i] = v

            # Physical-validity gate: a non-finite SED invalidates all the other
            # quantities for this sample (the SED itself keeps its values). The
            # SED of the current batch is used when available (construction),
            # otherwise the canonical one stored at construction is used, so
            # quantities added later via ``add_property`` are gated identically.
            sed_ref = out['SED'] if 'SED' in out else getattr( self, 'SED', None )
            if sed_ref is not None and not numpy.isfinite( sed_ref[i] ).any() :
                for k in keys :
                    if k != 'SED' :
                        out[k][i] = _sentinel

        # No sample produced a value and the fixed configuration was invalid.
        if out is None :
            out = { k : numpy.full( ( self.size, ), _sentinel, dtype = float )
                    for k in keys }

        for k in keys :
            setattr( self, k, out[k] )
            if k not in self._derived :
                self._derived += [ k ]
        return

    def add_property ( self, func, name = None ) :
        """Compute and store one or more additional derived quantities.

        The new quantities are evaluated over all the stored samples exactly
        like the default ones (``SED``, ``Mstar``, ...) and become available to
        all the statistics methods (``get_mean``, ``get_quantile``, ...). They
        are also serialised by ``dump`` and restored by ``load``.

        Parameters
        ----------
        func : callable or dict
            Either a single callable ``f(model)`` or a dictionary
            ``{ name : callable }``. Each callable receives the model with the
            parameters of the current sample already set and returns the derived
            value for that sample (scalar or array). The model's age is
            available as ``model.age``.
        name : str, optional
            Name under which a single callable is stored. Ignored when ``func``
            is a dictionary. When ``None`` and ``func`` is a single callable an
            automatic name ``custom{N}`` is assigned.

        Returns
        -------
        : list
            The names of the quantities that have been added.

        Examples
        --------
        >>> res.add_property( lambda model : model.sfh.Mstar( model.age ),
        ...                   name = 'Mstar_check' )
        ['Mstar_check']
        >>> res.add_property( { 'Lbol' : lambda model : model.get_SED().sum() } )
        ['Lbol']
        """
        if not isinstance( func, MM ) :
            if name is not None :
                if not isinstance( name, str ) :
                    raise TypeError( 'Argument ``name`` should be a string.' )
                func = { name : func }
            else :
                icust = 0
                while hasattr( self, f'custom{icust:d}' ) :
                    icust += 1
                func = { f'custom{icust:d}' : func }
        for k, f in func.items() :
            if not hasattr( f, '__call__' ) :
                raise AttributeError(
                    'Argument ``func`` should be a callable or a dictionary of '
                    f'{{name:callable}}; "{k}" is not callable.'
                )
        self._compute_properties( func )
        return list( func.keys() )

    def dump ( self ) :
        ret = dict(
            # Models' architecture
            model = self._mod,
            handler = self._han,
            data = self._obs,
            noise = self._noise,
            # Sampling run hyperparameters
            sampler_name = self.sampler,
            size = self.size, Ndof = self.Ndof,
            logz    = self.logz,
            logzerr = self.logzerr,
            # Sampling run stored quantities
            logl    = self.logl,
            samples = self.samples,
            weights = self.weights,
            wnot0   = self.wnot0,
            # Names of the stored derived quantities (pipe-joined, mirroring
            # the convention used in galapy.Handlers.ModelParameters.dump)
            derived = '|'.join( self._derived ),
        )
        # Derived quantities (defaults plus any added via ``add_property``)
        for k in self._derived :
            ret[k] = getattr( self, k )
        return ret

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
        ret.logz    = dictionary.get( 'logz',    None )
        ret.logzerr = dictionary.get( 'logzerr', None )

        # Derived quantities. Files written before the introduction of the
        # ``derived`` key always stored exactly the nine default quantities.
        derived = dictionary.get( 'derived', None )
        if derived is None :
            derived = [ 'SED', 'Mstar', 'Mdust', 'Mgas',
                        'Zstar', 'Zgas', 'SFR', 'TMC', 'TDD' ]
        else :
            derived = derived.split('|')
        ret._derived = list( derived )
        for k in derived :
            setattr( ret, k, dictionary[k] )

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
