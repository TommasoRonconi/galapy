# External imports
# from abc import ABC, abstractmethod
from collections.abc import MutableMapping as MM
import numpy

###################################################################################

def trap_int ( xx, yy ) :
    """ Trapezoid integration

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
    return numpy.sum( 0.5 * ( yy[:-1] + yy[1:] ) * ( xx[1:] - xx[:-1] ) )

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
    """ Function finding the indexes of array closer to the
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
    return El**(-gamma+3) * numpy.exp(-El/Ecut) 

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

# def unwrap_keys ( d ) :
#     """ Generator yielding the keys of a nested dictionary

#     Paramters
#     ---------
#     d : dictionary
#       A nested dictionary (whathever object inheriting from 'dict' 
#       or collections.abc.MutableMapping or similar)
    
#     Yields
#     ------
#     : list of keys    

#     Examples
#     --------
#     >>> d = { 'a' : 1, 'b' : { 'c' : 2, 'd' : 3 } }
#     >>> for i in unwrap_keys( d ) : print( i )
#     ['a']
#     ['b', 'c']
#     ['b', 'd']
#     """
#     from collections.abc import MutableMapping
#     for k, v in d.items() :
#         if isinstance( v, MutableMapping ) :
#             for t in unwrap_keys( v ) :
#                 yield [ k ] + t
#         else :
#             yield [ k ]

# def unwrap_values ( d ) :
#     """ Generator yielding the values of a nested dictionary

#     Paramters
#     ---------
#     d : dictionary
#       A nested dictionary (whathever object inheriting from 'dict' 
#       or collections.abc.MutableMapping or similar)
    
#     Yields
#     ------
#     : values

#     Examples
#     --------
#     >>> d = { 'a' : 1, 'b' : { 'c' : 2, 'd' : 3 } }
#     >>> for i in unwrap_values( d ) : print( i )
#     1
#     2
#     3
#     """
#     from collections.abc import MutableMapping
#     for k, v in d.items() :
#         if isinstance( v, MutableMapping ) :
#             for t in unwrap_values( v ) :
#                 yield t
#         else :
#             yield v

# def unwrap_items ( d ) :
#     """ Generator yielding the items of a nested dictionary

#     Paramters
#     ---------
#     d : dictionary
#       A nested dictionary (whathever object inheriting from 'dict' 
#       or collections.abc.MutableMapping or similar)
    
#     Yields
#     ------
#     : items
#       a 2d tuple in which the first element is the list of nested keys
#       and the second is the associated value

#     Examples
#     --------
#     >>> d = { 'a' : 1, 'b' : { 'c' : 2, 'd' : 3 } }
#     >>> for i in unwrap_items( d ) : print( i )
#     (['a'], 1)
#     (['b', 'c'], 2)
#     (['b', 'd'], 3)
#     """
#     from collections.abc import MutableMapping
#     for k, v in d.items() :
#         if isinstance( v, MutableMapping ) :
#             for tk, tv in unwrap_items( v ) :
#                 yield [ k ] + tk, tv
#         else :
#             yield [ k ], v

###################################################################################

