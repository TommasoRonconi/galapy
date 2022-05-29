""" The Synchrotron module implements a generic parameterized synchrotron emission.
"""

# External imports
import numpy

# Internal imports
from .SYN_core import CSYN

syn_tunables = ( 'alpha_syn', 'nu_self_syn' )

def syn_build_params ( **kwargs ) :
    """ Builds the parameters dictionary for the generic synchrotron.
    """
    
    out = {
        'alpha_syn'   : 0.75,
        'nu_self_syn' : 0.2,  # [GHz]
    }
    for k in set( out.keys() ).intersection(kwargs.keys()) :
        out[k] = kwargs[k]
        
    return out

class SYN () :
    """ Class wrapping the C-core implementation of a generic Synchrotron emission type.   
    
    Parameters
    ----------
    ll : float array
      the wavelenght grid where the synchrotron emission is computed
    """

    def __init__ ( self, ll, **kwargs ) :

        self.l = numpy.ascontiguousarray( ll )
        self.core = CSYN( self.l )
        self.params = syn_build_params( **kwargs )
        self.set_parameters()
        
    def set_parameters ( self, **kwargs ) :
        r""" Function for setting the parameters of the model.
        
        Returns
        -------
        : None
        """

        self.params.update( kwargs )
        self.core.set_params( numpy.asarray( [
            self.params[k]
            for k in syn_tunables
        ], dtype=float) )
        return;

    def opt_depth ( self, il = None ) :        

        if il is None :
            il = numpy.arange( self.l.size, dtype = numpy.uint64 )
        il = numpy.asarray( il )
        scalar_input = False
        if il.ndim == 0 :
            il = il[None] # makes il 1D
            scalar_input = True

        ret = numpy.array( [ self.core.opt_depth( _i ) for _i in il ] )
        if scalar_input :
            return ret.item()
        return ret

    def emission ( self, SynNorm, il = None ) :
        
        if il is None :
            il = numpy.arange( len(self.l), dtype = numpy.uint64 )
        il = numpy.asarray( il )
        scalar_input = False
        if il.ndim == 0 :
            il = il[None] # makes il 1D
            scalar_input = True
        
        ret = self.core.emission( il, SynNorm )
        if scalar_input :
            return ret.item()
        return ret
    
