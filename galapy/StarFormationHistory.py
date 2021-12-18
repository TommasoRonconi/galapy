# External imports
import numpy

# Internal imports
from .SFH_core import CSFH
#from .internal.utils import gxyComp

sfh_tunables = {
    # In-Situ SF model
    'insitu' : [ 'psi_max', 'tau_star' ],
    
    # Constant SF model
    'constant' : [ 'psi', 'Mdust', 'Zgs' ],
    
    # Delayed-Exponential model
    'delayedexp' : ['psi_norm', 'k_shape', 'tau_star', 'Mdust', 'Zgs' ],
    
    # Log-Normal model
    'lognormal' : ['psi_norm', 'sigma_star', 'tau_star', 'Mdust', 'Zgs' ]
}
""" Dictionary of tunable parameters
"""
    
def sfh_build_params ( tau_quench = 2.e+10, model = 'insitu', **kwargs ) :
    """ Builds the parameter dictionary for a given SFH model.

    Paramters
    ---------

    tau_quench : float
       eventual time of quenching in units of year, 
       defaults to the arbitrary large value of $2 \cdot 10^9$ years
    model : string
       SFH model
    
    Keyword Arguments
    -----------------
    
    The parameterization depends on the chosen model.

    Returns
    -------
    : dict
       Dictionary containing the parameterization of the chosen SFH model.
       All the parameters that have not been passed to the function are defined
       with their default value.
    """

    _models = {
        # In-Situ SF model
        'insitu' : {
            'psi_max' :  100.,
            'tau_star' : 3.e+8,
        },
        # Constant SF model
        'constant' : {
            'psi'   : 1.,
            'Mdust' : 1.e+8,
            'Zgs'   : 0.1,
        },
        # Delayed-Exponential model
        'delayedexp' : {
            'psi_norm' : 1.,
            'k_shape'  :  0.2,
            'tau_star' : 1.e+8,
            'Mdust'    : 1.e+8,
            'Zgs'      : 0.1,
        },
        # Log-Normal model
        'lognormal' : {
            'psi_norm'   :   100.,
            'sigma_star' : 2.,
            'tau_star'   :   3.e+8,
            'Mdust'      : 1.e+8,
            'Zgs'        : 0.1,
        }
    }
    out = { 'tau_quench' : tau_quench, 'model' : model }
    out.update( _models[ model ] )
    for k in set(out.keys()).intersection(kwargs.keys()) :
        out[k] = kwargs[k]

    return out
        

class SFH () :
    """ Class wrapping the C-core implementation of the Star Formation History type.    
    The possible models to choose are:

    #. 'insitu'
    #. 'constant'
    #. 'delayedexp'
    #. 'lognormal'
    #. 'burst'

    Parameters
    ----------
    tau_quench : float
      Eventual abrupt quenching time for star formation. 
      Should be expressed in years. It refers to the age of the galaxy,
      i.e. it has to be intended as the time passed from the formation of the galaxy.

    model : string
      One among ( 'insitu', 'constant', 'delayedexp', 'lognormal', 'burst' ).
      Default is 'insitu'.

    Note
    ----
    Not for SED fitting, use the galaxy class
    """

    def __init__ ( self, tau_quench = 2.e+20, model = 'insitu', **kwargs ) :
        
        self.core = CSFH( tau_quench, model )
        self.params = sfh_build_params( tau_quench, model, **kwargs )
        self.tunable = set( self.params.keys() )
        self.tunable.remove('model')
        self.set_parameters()

        # steal docstrings from C-core:
        self.__call__.__func__.__doc__ = self.core.__call__.__doc__
        self.Mstar.__func__.__doc__    = self.core.Mstar.__doc__
        self.Mdust.__func__.__doc__    = self.core.Mdust.__doc__
        self.Mgas.__func__.__doc__     = self.core.Mgas.__doc__
        self.Zgas.__func__.__doc__     = self.core.Zgas.__doc__
        self.Zstar.__func__.__doc__    = self.core.Zstar.__doc__
        

    def __call__ ( self, tau ) :
        """ Method for returning the Star Formation Rate (SFR) at given age.

        Parameters
        ----------
        tau : float or array-like
          The time in the galaxy SFH in which to compute the SFR

        Returns
        -------
        SFR : float or array-like
          SFR at given time from galaxy formation.
        """
        tau = numpy.asarray( tau )
        scalar = False
        if tau.ndim == 0 :
            tau = tau[None] # makes 'tau' 1D
            scalar = True
        ret = self.core( tau )
        if scalar :
            return ret.item()
        return ret

    def set_parameters ( self, tau_quench = None, **kwargs ) :
        r""" Function for setting the parameters of the model.

        Parameters
        ----------
        tau_quench : float
          Eventual time of abrupt quenching event, stopping 
          star formation forever.
        \**kwargs : 
          see below

        Keyword Arguments
        -----------------
        All the parameters the user wants, will consider only those with
        a valid key for the SFH model chosen.
        
        Returns
        -------
        : None
        """
        if tau_quench is not None :
            self.core.set_tau_quench( tau_quench )
        self.params.update( kwargs )
        self.core.set_params( numpy.asarray( [
            self.params[k]
            for k in sfh_tunables[self.params['model']]
        ], dtype=float) )
        return;

    def Mstar ( self, tau, npoints = 100 ) :
        """ 

        Parameters
        ----------
        tau : array-like or float
        
        """
        tau = numpy.asarray( tau )
        scalar_input = False
        if tau.ndim == 0 :
            tau = tau[None] # makes 'tau' 1D
            scalar_input = True
        ret = numpy.asarray( [ self.core.Mstar( _t, npoints )
                               for _t in tau ],
                             dtype=numpy.float64 )
        if scalar_input :
            return ret.item()
        return ret

    def Mdust ( self, tau ) :
        tau = numpy.asarray( tau )
        scalar_input = False
        if tau.ndim == 0 :
            tau = tau[None] # makes 'tau' 1D
            scalar_input = True
        ret = numpy.asarray( [ self.core.Mdust( _t )
                               for _t in tau ],
                             dtype=numpy.float64 )
        if scalar_input :
            return ret.item()
        return ret

    def Mgas ( self, tau ) :
        tau = numpy.asarray( tau )
        scalar_input = False
        if tau.ndim == 0 :
            tau = tau[None] # makes 'tau' 1D
            scalar_input = True
        ret = numpy.asarray( [ self.core.Mgas( _t )
                               for _t in tau ],
                             dtype=numpy.float64 )
        if scalar_input :
            return ret.item()
        return ret

    def Zgas ( self, tau ) :
        tau = numpy.asarray( tau )
        scalar_input = False
        if tau.ndim == 0 :
            tau = tau[None] # makes 'tau' 1D
            scalar_input = True
        ret = numpy.asarray( [ self.core.Zgas( _t )
                               for _t in tau ],
                             dtype=numpy.float64 )
        if scalar_input :
            return ret.item()
        return ret

    def Zstar ( self, tau ) :
        tau = numpy.asarray( tau )
        scalar_input = False
        if tau.ndim == 0 :
            tau = tau[None] # makes 'tau' 1D
            scalar_input = True
        ret = numpy.asarray( [ self.core.Zstar( _t )
                               for _t in tau ],
                             dtype=numpy.float64 )
        if scalar_input :
            return ret.item()
        return ret
    
