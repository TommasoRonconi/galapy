""" The Synchrotron module implements a generic parameterized synchrotron emission.
"""

# External imports
import numpy

# Internal imports
from galapy.internal.constants import clight
from .SYN_core import CSYN

__all__ = [ 'SYN', 'SNSYN' ]

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
      the wavelength grid where the synchrotron emission is computed
    """

    def __init__ ( self, ll, **kwargs ) :

        self.l = numpy.ascontiguousarray( ll )
        self.AngUnit = clight['A/s'] / self.l**2
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

    def energy ( self, SynNorm, il = None ) :
        """ Returns the normalized total energy radiated at given wavelength.       
        """
        
        if il is None :
            il = numpy.arange( len(self.l), dtype = numpy.uint64 )
        il = numpy.asarray( il )
        scalar_input = False
        if il.ndim == 0 :
            il = il[None] # makes il 1D
            scalar_input = True
        
        ret = self.core.energy( il, SynNorm )
        if scalar_input :
            return ret.item()
        return ret

    def emission ( self, SynNorm, il = None ) :
        
        return self.energy( SynNorm, il ) * self.AngUnit
    
class SNSYN ( SYN ) :

    NormFact = 1.e+30 # [ erg s^-1 Hz^-1 ]
    Lsyn0    = 3.e+28 # [ erg s^-1 Hz^-1 ]
    eta      = 2.     # [ adimensional ] 

    def __init__ ( self, *args, **kwargs ) :

        super().__init__( *args, **kwargs )
        self.params['RCCSN'] = 0.
        self.params.update(kwargs)

    def emission ( self, il = None ) :
        
        Lsyn = super().energy( self.params['RCCSN'] * SNSYN.NormFact, il )
        wn0 = Lsyn > 0.
        Lsyn[wn0] /= 1 + ( SNSYN.Lsyn0 / Lsyn[wn0] )**SNSYN.eta

        if il is None :
            return Lsyn * self.AngUnit
        return Lsyn * self.AngUnit[ il ]
