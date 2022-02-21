""" Galapy module providing the Nebular Emission component
"""

# External imports
import numpy

# Internal imports
from .NFF_core import CNFF

nff_tunables = ( 'Zgas', 'Zi' )

def nff_build_params ( **kwargs ) :
    
    out = {
        'Zgas' : 0.02,
        'Zi'   : 1.,
    }
    for k in set( out.keys() ).intersection( kwargs.keys() ) :
        out[ k ] = kwargs[ k ]

    return out

class NFF () :
    """ Class wrapping the C-core implementation of the Nebular Free-Free emission type.   
    
    Parameters
    ----------
    ll : float array
      the wavelenght grid where the free-free emission is computed
    """
    
    def __init__ ( self, ll, **kwargs ) :

        self.l = numpy.ascontiguousarray( ll )
        self.core = CNFF( self.l )
        self.params = nff_build_params( **kwargs )
        self.set_parameters()

        
    def set_parameters ( self, **kwargs ) :
        r""" Function for setting the parameters of the model.
        
        Returns
        -------
        : None
        """

        self.params.update( kwargs )
        # self.params.update( { k : v
        #                       for k, v in kwargs.items()
        #                       if k in nff_tunables } )
        self.core.set_params( numpy.asarray( [
            self.params[k]
            for k in nff_tunables
        ], dtype=float) )
        return;

    def gff ( self, il = None, Te = None ) :
        
        if il is None :
            il = numpy.arange( self.l.size, dtype = numpy.uint64 )
        il = numpy.asarray( il )
        scalar_input = False
        if il.ndim == 0 :
            il = il[None] # makes il 1D
            scalar_input = True
        if Te is None :
            Te = self.Te()

        ret = numpy.array( [ self.core.gff( _i, Te ) for _i in il ] )
        if scalar_input :
            return ret.item()
        return ret

    def Te ( self, Zgas = None ) :
        
        if Zgas is None :
            Zgas = self.params[ 'Zgas' ]
        Zgas = numpy.asarray( Zgas )
        scalar_input = False
        if Zgas.ndim == 0 :
            Zgas = Zgas[None] # makes il 1D
            scalar_input = True
            
        ret = numpy.array( [ self.core.Te( _Z ) for _Z in Zgas ] )
        if scalar_input :
            return ret.item()
        return ret

    def emission ( self, Q_H, il = None ) :
        
        if il is None :
            il = numpy.arange( len(self.l), dtype = numpy.uint64 )
        il = numpy.asarray( il )
        scalar_input = False
        if il.ndim == 0 :
            il = il[None] # makes il 1D
            scalar_input = True
        
        ret = self.core.emission( il, Q_H )
        if scalar_input :
            return ret.item()
        return ret
        
    
