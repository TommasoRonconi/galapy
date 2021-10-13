# External imports
import numpy

# Internal imports
from .internal.CPySFH import *
#from .internal.utils import gxyComp

class SFH () :
    """ Class defining the model of Star Formation History.
    The possible models to choose are
    

    Parameters
    ----------
    tau_quench : float
      Eventual abrupt quenching time for star formation. 
      Should be expressed in years. Refers to the age of the galaxy.

    model : string
      One among ( 'insitu', 'constant', 'delayedexp', 'lognormal', 'burst' ).
      Default is 'insitu'.
    """

    def __init__ ( self, *args, **kwargs ) :
        self.csfh = CSFH( *args, **kwargs )
        self.__call__.__func__.__doc__ = self.csfh.__call__.__doc__
        self.Mstar.__func__.__doc__    = self.csfh.Mstar.__doc__
        self.Mdust.__func__.__doc__    = self.csfh.Mdust.__doc__
        self.Mgas.__func__.__doc__     = self.csfh.Mgas.__doc__
        self.Zgas.__func__.__doc__     = self.csfh.Zgas.__doc__
        self.Zstar.__func__.__doc__    = self.csfh.Zstar.__doc__
        

    def __call__ ( self, tau ) :
        return numpy.asarray( self.csfh( tau ) )

    def set_parameters ( self, params ) :
        self.csfh.set_params( numpy.asarray( params, dtype = float ) )
        return;

    def Mstar ( self, tau, npoints = 100 ) :
        """ 

        Parameters
        ----------
        tau : array-like or float
        
        """
        if hasattr( tau, "__len__") :
            return numpy.array( [ self.csfh.Mstar( _t, npoints )
                                  for _t in tau ] )
        return self.csfh.Mstar( tau, npoints )

    def Mdust ( self, tau ) :
        pass

    def Mgas ( self, tau ) :
        pass

    def Zgas ( self, tau ) :
        pass

    def Zstar ( self, tau ) :
        pass


        

    
    
    
