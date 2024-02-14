# External imports
# from abc import ABC, abstractmethod
from collections.abc import MutableMapping as MM
import numpy

###################################################################################

def trap_int ( xx, yy ) :
    r""" Trapezoid integration

    Parameters
    ----------
    xx : array-like
       x-domain grid
    yy : array-like
       y-domain grid
    
    Returns
    -------
    : float 
       the integral along the whole x-domain
    """
    xx = numpy.asarray(xx)
    yy = numpy.asarray(yy)
    return numpy.sum( 0.5 * ( yy[:-1] + yy[1:] ) * ( xx[1:] - xx[:-1] ), axis = 0 )

###################################################################################

class FlagVal () :
    def __init__ ( self, value, flag ) :
        self.value = value
        self.flag  = flag
    def __repr__ ( self ) :
        return f'{type(self).__name__}({self.value},"{self.flag}")'
    def __str__ ( self ) :
        return f'{self.flag}:{self.value}'

###################################################################################

def find_nearest ( array, value ) :
    r""" Function finding the indexes of array closer to the
    values in value.
    
    Parameters
    ----------
    array : array-like 
    
    value : array-like or scalar
    
    Return
    ------
    : numpy-array
      list of indexes of elements in array closer to 
      values in value.
    
    Warning
    -------
    Uniqueness of elements in the returned array is not guaranteed.
    """
    import numpy
    
    array = numpy.asarray( array )
    value = numpy.asarray( value )
    scalar_input = False
    if value.ndim == 0 :
        value = value[None] # makes 'value' 1D
        scalar_input = True
        
    idx = [ ( numpy.abs( array - _v ) ).argmin() for _v in value ] 

    if scalar_input :
        return numpy.squeeze( idx )
    return numpy.asarray( idx )

###################################################################################

def powerlaw_exp_cutoff ( El, gamma, Ecut ) :
    r""" Returns a powerlaw with exponential cut-off of the form

    .. math::
      
      y = E_\lambda^{-\gamma+3} e^{E_\lambda/E_\text{cut}}

    Parameters
    ----------
    El : scalar or array-like
      the x-values
    gamma : scalar
      the spectral index
    Ecut : scalar
      the cut-off x-scale
    
    Returns
    -------
    : scalar or array-like
    """
    return El**(-gamma+3) * numpy.exp(-El/Ecut) 

###################################################################################

def poly_N ( xx, coeff ) :
    r""" Method for computing N-order polynomia.
    The order of the polynomium is set by the length of the `coeff` list.
    
    The method computes

    .. math::
    
       y = \sum_0^N c_i * x^{N-i}
    
    Parameters
    ----------
    xx : scalar or array-like
      the x-values 
    coeff : list
      the values of the coefficients of the polynomium. The list must be
      ordered, with the first element corresponding to the coefficient multiplying
      the highest power of the x-variable and the last coefficient being the 
      scalar coefficient
    
    Returns
    -------
    : scalar or array-like
      the y-values
    """
    yy = 0.
    for _c in coeff :
        yy = yy * xx + _c
    return yy

###################################################################################

def unwrap_keys ( d ) :
    r""" Generator yielding the keys of a nested dictionary

    Parameters
    ----------
    d : dictionary
      A nested dictionary (whathever object inheriting from 'dict' 
      or collections.abc.MutableMapping or similar)
    
    Yields
    ------
    : list of keys    

    Examples
    --------
    >>> d = { 'a' : 1, 'b' : { 'c' : 2, 'd' : 3 } }
    >>> for i in unwrap_keys( d ) : print( i )
    ['a']
    ['b', 'c']
    ['b', 'd']
    """
    from collections.abc import MutableMapping
    for k, v in d.items() :
        if isinstance( v, MutableMapping ) :
            for t in unwrap_keys( v ) :
                yield [ k ] + t
        else :
            yield [ k ]

def unwrap_values ( d ) :
    r""" Generator yielding the values of a nested dictionary

    Parameters
    ----------
    d : dictionary
      A nested dictionary (whathever object inheriting from 'dict' 
      or collections.abc.MutableMapping or similar)
    
    Yields
    ------
    : values

    Examples
    --------
    >>> d = { 'a' : 1, 'b' : { 'c' : 2, 'd' : 3 } }
    >>> for i in unwrap_values( d ) : print( i )
    1
    2
    3
    """
    from collections.abc import MutableMapping
    for k, v in d.items() :
        if isinstance( v, MutableMapping ) :
            for t in unwrap_values( v ) :
                yield t
        else :
            yield v

def unwrap_items ( d ) :
    r""" Generator yielding the items of a nested dictionary

    Parameters
    ----------
    d : dictionary
      A nested dictionary (whathever object inheriting from 'dict' 
      or collections.abc.MutableMapping or similar)
    
    Yields
    ------
    : items
      a 2d tuple in which the first element is the list of nested keys
      and the second is the associated value

    Examples
    --------
    >>> d = { 'a' : 1, 'b' : { 'c' : 2, 'd' : 3 } }
    >>> for i in unwrap_items( d ) : print( i )
    (['a'], 1)
    (['b', 'c'], 2)
    (['b', 'd'], 3)
    """
    from collections.abc import MutableMapping
    for k, v in d.items() :
        if isinstance( v, MutableMapping ) :
            for tk, tv in unwrap_items( v ) :
                yield [ k ] + tk, tv
        else :
            yield [ k ], v

###################################################################################

def set_nested ( d, kl, v ) :
    r""" Recursive function for setting values in a nested dictionary.
    
    Parameters
    ----------
    d : dictionary
        the i/o dictionary
    kl : list
        a list of keywords to traverse the nested dictionary
        the first element is the highest position in the hierarchy
        of the nested dictionary
    v : whatever
        the value to set at corresponding nested sequence of keywords
    
    Returns
    -------
    : None
    
    Examples
    --------
    >>> d = {}
    >>> kl1 = ['first', 'second', 'third']
    >>> kl2 = ['a', 'b', 'c']
    >>> set_nested(d,kl1, 3.)
    >>> set_nested(d,kl2, True)
    >>> d
    {'first': {'second': {'third': 3.0}}, 'a': {'b': True}}
    
    If we want to add a new key-value pair to an already existing level:
    
    >>> kl3 = ['first', 'second', 'fourth']
    >>> set_nested(d, kl3, 4.)
    >>> d
    {'first': {'second': {'third': 3.0, 'fourth': 4.0}},'a': {'b': True}}
    """
    k = kl.pop(0)
    if len(kl) == 0 :
        d[k] = v
        return;
    try :
        set_nested( d[k], kl, v )
    except KeyError :
        d[k] = {}
        set_nested( d[k], kl, v )
    return;        

###################################################################################

def get_nested ( d, kl ) :
    if len( kl ) > 1 and isinstance( d, MM ) :
        return get_nested( d[ kl.pop( 0 ) ], kl )
    return d[ kl.pop( 0 ) ]

###################################################################################

def func_scalar_or_array ( var, function, *args, **kwargs ) :
    """ A function applying some function to a scalar or an ndarray.
    
    Parameters
    ----------
    var : scalar or ndarray
        The argument of the function
    function : callable
        A function taking a scalar as input
    args : sequence
        A list or tuple of additional arguments to be passed to function
    kwargs : dict
        A dictionary of additional keyword arguments to be passed to function
    
    Returns
    -------
    : scalar or ndarray
        Depends on the shape of ``var``: (i) if ``var`` is a scalar it is the result 
        of ``function(var)``; (ii) if ``var`` is an ndarray it is the result of 
        ``array([function(v) for v in var])``
    """
    var = numpy.asarray(var)
    scalar = False
    if var.ndim == 0 :
        var = var[None]
        scalar = True
    ret = numpy.array([function(v,*args,**kwargs) for v in var])
    if scalar :
        return ret[0]
    return ret

###################################################################################

def quantile_weighted ( values, quantiles, weights = None, 
                        values_sorted = False, old_style = False, 
                        axis=-1 ) :
    """ Very close to numpy.percentile, but supports weights.
    (partially stolen from stackoverflow: 
    https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy)

    Parameters
    ----------
    values : numpy.array 
        with data
    quantiles : array-like 
        quantiles needed (NOTE: quantiles should be in [0, 1]!)
    sample_weight : array-like 
        same length as ``values``
    values_sorted : bool
        if True, then will avoid sorting of initial array
    old_style : bool
        if True, will correct output to be consistent
        with numpy.percentile.
    axis : int
    
    Returns
    -------
    : numpy.array 
        computed quantiles.
    """
    values = numpy.array(values)
    quantiles = numpy.array(quantiles)
    if weights is None:
        weights = numpy.ones(len(values))
        
    weights = numpy.array(weights)
    if numpy.any(quantiles < 0) or numpy.any(quantiles > 1) :
        raise ValueError( 'Quantiles should be in [0, 1]' )        

    def get_weighted_1D( val, wgh ) :
        if not values_sorted :
            sorter = numpy.argsort(val)
            val = val[sorter]
            wgh = weights[sorter]
            
        weighted_quantiles = numpy.cumsum(wgh) - 0.5 * wgh
        if old_style : 
            # To be convenient with numpy.percentile
            weighted_quantiles -= weighted_quantiles[0]
            weighted_quantiles /= weighted_quantiles[-1]
        else :
            weighted_quantiles /= numpy.sum(wgh)
        return numpy.interp(quantiles, weighted_quantiles, val)
    
    if values.ndim > 1 :
        axes = numpy.delete(numpy.arange(values.ndim, dtype=numpy.int64), axis)
        res = numpy.empty( ( *[values.shape[a] for a in axes], *quantiles.shape ) )
        values = numpy.transpose(values, axes = numpy.append( axes, axis ))
        for i, val in enumerate(values) :
            res[i] = get_weighted_1D( val, weights )

        return res.T
    
    return get_weighted_1D( values, weights )

###################################################################################

def get_credible_interval(samples, idcentre, percent = 0.68, weights = None) :
    """Find lower and upper limits in sample array containing 
    given percent value around position.
    
    It splits the samples in 2 half-intervals:
    [lower, samples[idcentre]) and [samples[idcentre], upper) that are set to 
    contain 0.5*percent of all the weighted samples.
    If one of the two half-intervals does not have enough samples, it is considered
    un-defined. When this happens, the function consider it instead an upper/lower limit 
    on the sampled parameter and returns instead (None, upper-limit) or (lower-limit, None).
    
    Parameters
    ----------
    samples : ndarray or iterable
        the weighted samples over which to compute the estimator 
    idcentre : integer
        index in the samples iterable around which to compute the interval
    percent : float
        probability enclosed by the interval 
    weights :
        (Optional, default = None) weights of the samples, if not None
        this must have the same size of the samples argument
    
    Returns
    -------
    lower, upper : scalar floats
        [lower-, upper-) limits of the interval if the interval is defined.
        If it is not defined, the function will return either an upper limit or
        a lower limit, respectively.
    """
    
    if weights is None :
        weights = numpy.ones_like(samples)
    
    # sort the samples
    idsort = numpy.argsort(samples)
    
    # Normalize the weights to create a probability distribution
    weights = numpy.array(weights[idsort]) / numpy.sum(weights[idsort])

    # Calculate the cumulative distribution function (CDF)
    cdf = numpy.cumsum(weights)

    # Find the target value in the sorted array
    idx = numpy.nonzero(idsort == idcentre)

    # Find the lower bound based on the desired percentile
    lower_percentile = cdf[idx] - 0.5 * percent
    
    # account for upper limits
    if lower_percentile < 0.0 :
        upper_percentile = percent
        upper_idx = numpy.searchsorted(cdf, upper_percentile).item()
        return None, samples[idsort][upper_idx]
    
    # Find the upper bound based on the desired percentile
    upper_percentile = cdf[idx] + 0.5 * percent
    
    # account for lower limits
    if upper_percentile > 1.0 :
        lower_percentile = 1.0-percent
        lower_idx = numpy.searchsorted(cdf, lower_percentile).item()
        return samples[idsort][lower_idx], None

    # Find the lower and upper samples within the specified percentile range
    lower_idx = numpy.searchsorted(cdf, lower_percentile).item()
    upper_idx = numpy.searchsorted(cdf, upper_percentile).item()

    return samples[idsort][lower_idx], samples[idsort][upper_idx]

###################################################################################

def now_string () :
    """ Generates a string with the minute the function is called
    """
    from time import localtime as lt
    now = lt()
    return ( f'{now.tm_year:04}{now.tm_mon:02}{now.tm_mday:02}' +
             f'{now.tm_hour:02}{now.tm_min:02}' )

###################################################################################

def filter_strings ( inlist, fields ) :
    """Function for filtering a list of strings (uses ``fnmatch``)
    
    Parameters
    ----------
    inlist : sequence of str
        A list of strings to filter
    fields : str or sequence of str
        Either a single string or a sequence of strings. 
        Also accepts wildcards (e.g. ``which_params = 'sfh*'`` will show all the
        entries in ``inlist`` that contain the sub-string ``sfh``).
        
    Returns
    -------
    : str or sequence of str    
    
    Raises
    ------
    ValueError
    """
    
    import fnmatch
    
    if isinstance(fields, str) :
        fields = fnmatch.filter(inlist, fields)
    elif isinstance(fields, (list,tuple)) :
        fields = [item for sublist in [fnmatch.filter(inlist, s) for s in fields] for item in sublist]
    else :
        raise ValueError( 'fields argument should be either a string or a list of strings' )
        
    return fields
    
###################################################################################

def shorten_string ( string, maxlength = 10, splitchar = '.' ) :
    """Function for shortening a string, splitting it in components.
    
    Parameters
    ----------
    string : str
        Input string
    maxlength : int
        (default = 10) Maximum number of characters in the final string
    splitchar : str
        (default = '.') The character over which perform the splitting.
        
    Returns
    -------
    : str 
    The shortened string with length <= maxlength+len(splitchar)
    """
    
    components = string.split(splitchar)
    short_string = ''
    for component in numpy.flip(components) :
        if len(short_string) == 0 :
            short_string = component
        elif len(short_string) + len(component) + len(splitchar) < maxlength :
            short_string = splitchar.join([component, short_string])
        else :
            break
    return short_string
    
###################################################################################

def cat_to_dict ( infile, id_field = 'id', err_field = '_err', meta_fields = [], skip_fields = [] ) :
    """Converts an ASCII (e.g. Topcat-like) catalogue into a 2-levels dictionary.
    The 1st order dictionary contains a 2nd order dictionary for each entry in the catalogue.
    Each 2nd order dictionary contains data and meta-data about the named entry. 
    
    Parameters
    ----------
    infile : str
        Path to the ASCII file
    id_field : str
        (Optional, default='id') header name of the field containing the sources' ID
    err_field : str
        (Optional, default='_err') the sub-string identifying fields
        containing error measurements
    meta_fields : str or sequence of str
        (Optional, default empty list) which fields to consider meta-data
        Either a single string or a sequence of strings. 
        Also accepts wildcards (e.g. ``meta_fields = 'R*'`` will
        consider all the fields in header that start with ``R`` as meta-data).
    skip_fields : str or sequence of str
        (Optional, default empty list) which fields to skip
        Either a single string or a sequence of strings. 
        Also accepts wildcards (e.g. ``skip_fields = 'R*'`` will
        ignore all the fields in header that start with ``R``).
        
    
    Returns
    -------
    : dict
    """
    import re
    import warnings
    
    with open(infile, 'r') as f :
        lines = f.readlines()
    content = [ re.split( '\t+| +|,', line.rstrip() ) for line in lines ]
    
    gxydict = {}
    names = []
    for head, *gxys in zip(*content) :
        
        if head in filter_strings(content[0],skip_fields) :
            continue
        if head==id_field :
            names = gxys
            for name in names :
                gxydict[name] = {'bands':[], 'fluxes':[], 'errors':[]}
            continue
        for key, gxy in zip(names,gxys) :
            if numpy.any( [head==what 
                           for what in filter_strings(content[0], meta_fields)] ) :
                try :
                    gxydict[key][head] = float(gxy)
                except ValueError :
                    warnings.warn( f'{head} value {gxy} cannot be cast to float, skipping' )
                    continue
            elif err_field in head :
                try :
                    gxydict[key]['errors'].append(float(gxy))
                except ValueError :
                    continue
            else :
                try :
                    gxydict[key]['fluxes'].append(float(gxy))
                except ValueError :
                    continue
                gxydict[key]['bands'].append(head)

    for gxy in gxydict.values() :
        gxy['bands']  = numpy.asarray(gxy['bands'])
        gxy['fluxes'] = numpy.asarray(gxy['fluxes'])
        gxy['errors'] = numpy.asarray(gxy['errors'])
    
    return gxydict

###################################################################################

# def recurrent_return ( dd, keylist ) :
#     if len( keylist ) > 1 and isinstance( dd, MM ) :
#         return recurrent_return( dd[ keylist.pop( 0 ) ], keylist )
#     return dd[ keylist.pop( 0 ) ]

# def recurrent_update ( dd, **kwargs ) :
#     for k, v in kwargs.items() :
#         if isinstance(v,MM) :
#             rec_update( dd[k], **v )
#         elif k in dd.keys() :
#             dd[k] = v
#         else :
#             print( f'Key {k} not valid' )
#             continue
#     return;

###################################################################################

# class SmartAttribute ( object ) :
#     """ Descriptor object
    
#     Parameters
#     ----------
#     value : int, float, string, custom class, ... , optional
#         the value to set the described item at initialization
#         can be whathever python type, built-in or custom
#     valid : function 
#         a function that accepts a type( value ) object and
#         either returns it as is or raises a ValueError exception
#     """
    
#     def __init__ ( self, value=None, valid=None, string=str ) :
#         super( SmartAttribute, self ).__init__()
#         self.value = value
#         self.validate = valid
#         self.string = string
        
#     def __set__ ( self, obj, value ) :
#         """ Validates the input value before setting
#         """
#         self.value = self.validate( value )
        
#     def __get__ ( self, obj, objtype ) :
#         """ Returns the stored value
#         """
#         return self.value
    
#     def __repr__ ( self ) :
#         return f"{type(self).__name__}(value={self.value},valid={self.validate.__name__},string={self.string.__name__})"
    
#     def __str__ ( self ) :
#         return self.string( self.value )

###################################################################################

# class descAttrBase () :
#     """ Base class that allows described instance-attributes.
#     As described instance-attributes are parameters of the
#     derived classes, when one of them is changed it sets the
#     '_paramset' attribute to 'False'
#     (this is where the magic happens)
#     """

#     def __setattr__(self, attr, value):
#         try:
#             # Try invoking the descriptor protocol __set__ "manually"
#             got_attr = super().__getattribute__(attr)
#             got_attr.__set__(self, value)
#             # If class owns a '_paramset' attribute set it to False
#             # when changing one of the described attributes
#             try:
#                 _ = super().__getattribute__('_paramset')
#                 super().__setattr__('_paramset',False)
#             except AttributeError:
#                 pass
#         except AttributeError:
#             # Attribute is not a descriptor, just set it:
#             super().__setattr__(attr, value)

#     def __getattribute__(self, attr):
#         # If the attribute does not exist, super().__getattribute__()
#         # will raise an AttributeError
#         got_attr = super().__getattribute__(attr)
#         try:
#             # Try "manually" invoking the descriptor protocol __get__()
#             return got_attr.__get__(self, type(self))
#         except AttributeError:
#             # Attribute is not a descriptor, just return it:
#             return got_attr

###################################################################################

# class gxyComp ( descAttrBase, ABC ) :
#     """
#     """

#     def __init__ ( self ) :
#         self._paramset = False

#     def __bool__ ( self ) :
#         return self._paramset

#     def __and__ ( self, other ) :
#         return self._paramset.__and__( other )
#     def __or__ ( self, other ) :
#         return self._paramset.__or__( other )
#     # def __iand__ ( self, other ) :
#     #     return self._paramset.__and__( other )
#     # def __ior__ ( self, other ) :
#     #     return self._paramset.__or__( other )

#     @abstractmethod
#     def set_parameters ( self, *args, **kwargs ) :
#         pass

###################################################################################

# def set_nested(d, parent_obj):
#     """
#     Sets a collection of instance attributes stored in a nested dictionary
#     in the class instance parent_obj and in its sub-classes instances
#     recursevely.
    
#     Parameters
#     ----------
#     d : dictionary
#         Dictionary containing the instance attributes to be set.
#         In case of nested dictionaries the function will search the
#         class instance for attributes named as the key identifying 
#         the nested dictionary and eventually call itself recursively
#         to set its attributes.
#     parent_obj : class instance
#         the instance of a generic type where the new instance attributes 
#         will be set.
        
#     Examples
#     --------
#     >>> specs = { 'attr1' : 1, 'attr2' : 2 }
#     >>> class Foo () :
#     ...     pass
#     ... 
#     >>> foo = Foo()
#     >>> set_nested( specs, foo )
#     >>> foo.attr1
#     1
#     >>> foo.attr2
#     2
    
#     In case of nested classes
    
#     >>> nested_specs = { 'attr_bar' : 3, 'foo' : specs }
#     >>> class Bar () :
#     ...     def __init__ ( self ) :
#     ...             self.foo = Foo()
#     ... 
#     >>> bar = Bar()
#     >>> set_nested( nested_specs, bar )
#     >>> bar.attr_bar
#     3
#     >>> bar.foo.attr1, bar.foo.attr2
#     (1, 2)
#     """
#     for k, v in d.items():
#         if isinstance(v, collections.abc.MutableMapping):
#             set_nested( v, getattr( parent_obj, k ) )
#         else:
#             setattr( parent_obj, k, v )
#     return;

###################################################################################

