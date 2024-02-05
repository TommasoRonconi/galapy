""" GalaPy module handling Star Formation History (SFH) models.
"""

# External imports
import numpy

# Internal imports
from .SFH_core import CSFH
#from .internal.utils import gxyComp

sfh_tunables = {
    # In-Situ SF model
    'insitu' : ( 'psi_max', 'tau_star' ),
    
    # Constant SF model
    'constant' : ( 'psi', 'Mdust', 'Zgxy' ),
    
    # Delayed-Exponential model
    'delayedexp' : ('psi_norm', 'k_shape', 'tau_star', 'Mdust', 'Zgxy' ),
    
    # Log-Normal model
    'lognormal' : ('psi_norm', 'sigma_star', 'tau_star', 'Mdust', 'Zgxy' ),
    
    # Interpolated model
    'interpolated' : ('Mdust', 'Zgxy' )
}
""" Dictionary of tunable parameters
"""
    
def sfh_build_params ( tau_quench = 2.e+10, model = 'insitu', **kwargs ) :
    """ Builds the parameter dictionary for a given SFH model.

    Parameters
    ----------
    tau_quench : float
       eventual time of quenching in units of year, 
       defaults to the arbitrary large value of $2 \cdot 10^9$ years
    model : string
       SFH model
    **kwargs : dictionary, optional
       The free-parameters of the SFH model chosen (with argument ``model``)
       The parameterization depends on the chosen model.
    
    Keyword Arguments
    -----------------
    psi_max : float
      ``model = 'insitu'``, optional, default = 100
    tau_star : float
      ``model = ['insitu','delayedexp', 'lognormal']``, optional, 
      default = 3.e+8 if ``insitu`` or ``lognormal``, = 1.e+8 if ``delayedexp``
    psi : float
      ``model = 'constant'``, optional, default = 1.0
    psi_norm : float
      ``model = ['delayedexp', 'lognormal']``, optional, 
      default = 1.0 if ``delayedexp``, = 100 if ``lognormal``
    k_shape : float
      ``model = 'delayedexp'``, optional, default = 0.2
    sigma_star : float
      ``model = 'lognormal'``, optional, default = 2.0
    Mdust : float
      ``model = ['constant','delayedexp','lognormal','interpolated']``, optional, default = 1.e+8
    Zgxy : float
      ``model = ['constant','delayedexp','lognormal','interpolated']``, optional, default = 0.01

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
            'psi_max'  : 100.,
            'tau_star' : 3.e+8,
        },
        # Constant SF model
        'constant' : {
            'psi'   : 1.,
            'Mdust' : 1.e+8,
            'Zgxy'  : 0.01,
        },
        # Delayed-Exponential model
        'delayedexp' : {
            'psi_norm' : 1.,
            'k_shape'  : 0.2,
            'tau_star' : 1.e+8,
            'Mdust'    : 1.e+8,
            'Zgxy'     : 0.01,
        },
        # Log-Normal model
        'lognormal' : {
            'psi_norm'   : 100.,
            'sigma_star' : 2.,
            'tau_star'   : 3.e+8,
            'Mdust'      : 1.e+8,
            'Zgxy'       : 0.01,
        },
        # Interpolated model
        'interpolated' : {
            'Mdust'      : 1.e+8,
            'Zgxy'       : 0.01,
        }
    }
    out = { 'tau_quench' : tau_quench, 'model' : model }
    out.update( _models[ model ] )
    for k in set(out.keys()).intersection(kwargs.keys()) :
        out[k] = kwargs[k]

    return out
        

class SFH () :
    r""" Class wrapping the C-core implementation of the Star Formation History type.   
 
    The possible models to choose are:

    #. ``'insitu'``
       
       .. math ::
          
          \psi(\tau) = \psi_\text{max} \left[ \exp(-x) - \exp(-s\gamma x) \right]
       
       with :math:`x = \tau / (s\tau_\star)`, :math:`s=3`, 
       :math:`\gamma = 1 - \mathcal{R} + 3 \psi_\text{max}^{-0.3}`, 
       :math:`\mathcal{R} = 0.45` the instantaneous recycling factor, 
       :math:`\tau_\star =`:code:`tau_star` and :math:`\psi_\text{max} =`:code:`psi_max`  
    
    #. ``'constant'``
       
       .. math ::
          
          \psi(\tau) = \psi_0
       
       with :math:`\psi_0 =`:code:`psi`
    
    #. ``'delayedexp'``
       
       .. math ::
          
          \psi(\tau) = \psi_\text{norm} \tau^{\kappa}\, \exp{(-\tau/\tau_\star)}
       
       with :math:`\psi_\text{norm} =`:code:`psi_norm`, 
       :math:`\kappa =`:code:`k_shape` and :math:`\tau_\star =`:code:`tau_star`
    
    #. ``'lognormal'``
       
       .. math ::
          
          \psi(\tau) = \psi_\text{norm} \dfrac{1}{\tau}\, \dfrac{1}{\sqrt{2\pi\sigma_\star^2}}\, 
          \exp\left[-\dfrac{\ln^2(\tau/\tau_\star)}{2\,\sigma_\star^2}\right]
       
       with :math:`\psi_\text{norm} =`:code:`psi_norm`, :math:`\sigma_\star =`:code:`sigma_star` 
       and :math:`\tau_\star =`:code:`tau_star`
    
    #. ``'interpolated'``: from a pre-computed grid of SFRs computed at given times.
       If this model is chosen, two positional arguments are mandatory, 
       i.e. the two grids ``args = (tau,sfr)`` with ``len(tau)==len(sfr)``, see below.

    Parameters
    ----------
    *args :
      positional arguments are optional, if passed they should be the two-iterables 
      corresponding to gridded values of time and SFR to interpolate over: ``args = (tau, sfr)``
      with ``len(tau)==len(sfr)``. Note that in this case the 'model' argument will be over-written
      and fixed to 'interpolated'
    tau_quench : float
      Eventual abrupt quenching time for star formation. 
      Should be expressed in years. It refers to the age of the galaxy,
      i.e. it has to be intended as the time passed from the formation of the galaxy.

    model : string
      One among ( 'insitu', 'constant', 'delayedexp', 'lognormal', 'interpolated', 'burst' ).
      Default is 'insitu'.

    **kwargs : dictionary, optional
       The free-parameters of the SFH model chosen (with argument ``model``)
       The parameterization depends on the chosen model.
    
    Keyword Arguments
    -----------------
    psi_max : float
      ``model = 'insitu'``, optional, default = 100
    tau_star : float
      ``model = ['insitu','delayedexp', 'lognormal']``, optional, 
      default = 3.e+8 if ``insitu`` or ``lognormal``, = 1.e+8 if ``delayedexp``
    psi : float
      ``model = 'constant'``, optional, default = 1.0
    psi_norm : float
      ``model = ['delayedexp', 'lognormal']``, optional, 
      default = 1.0 if ``delayedexp``, = 100 if ``lognormal``
    k_shape : float
      ``model = 'delayedexp'``, optional, default = 0.2
    sigma_star : float
      ``model = 'lognormal'``, optional, default = 2.0
    Mdust : float
      ``model = ['constant','delayedexp','lognormal','interpolated']``, optional, default = 1.e+8
    Zgxy : float
      ``model = ['constant','delayedexp','lognormal','interpolated']``, optional, default = 0.01

    Note
    ----
    Not for SED fitting, use the galaxy class
    """

    def __init__ ( self, *args, tau_quench = 2.e+20, model = 'insitu', **kwargs ) :

        if len(args) == 2 :
            try :
                tau = numpy.ascontiguousarray(args[0])
                sfr = numpy.ascontiguousarray(args[1])
                model = 'interpolated'
            except :
                raise;
        if 'tau' in kwargs and 'sfr' in kwargs :
            try :
                tau = numpy.ascontiguousarray(kwargs['tau'])
                sfr = numpy.ascontiguousarray(kwargs['sfr'])
                model = 'interpolated'
            except :
                raise;            
            
        self.core = CSFH( tau_quench, model )
        self.params = sfh_build_params( tau_quench, model, **kwargs )
        self.tunable = set( self.params.keys() )
        self.tunable.remove('model')
        if model == 'interpolated' :
            self.core.set_interpolator(tau, sfr)
        self.set_parameters()

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
        **kwargs : dictionary, optional
          The free-parameters of the SFH model chosen (with argument ``model``)
          The parameterization depends on the chosen model.
    
        Keyword Arguments
        -----------------
        psi_max : float
          ``model = 'insitu'``, optional, default = 100
        tau_star : float
          ``model = ['insitu','delayedexp', 'lognormal']``, optional, 
          default = 3.e+8 if ``insitu`` or ``lognormal``, = 1.e+8 if ``delayedexp``
        psi : float
          ``model = 'constant'``, optional, default = 1.0
        psi_norm : float
          ``model = ['delayedexp', 'lognormal']``, optional, 
          default = 1.0 if ``delayedexp``, = 100 if ``lognormal``
        k_shape : float
          ``model = 'delayedexp'``, optional, default = 0.2
        sigma_star : float
          ``model = 'lognormal'``, optional, default = 2.0
        Mdust : float
          ``model = ['constant','delayedexp','lognormal','interpolated']``, optional, default = 1.e+8
        Zgxy : float
          ``model = ['constant','delayedexp','lognormal','interpolated']``, optional, default = 0.01
        
        Returns
        -------
        : None
        """
        if tau_quench is not None :
            self.core.set_tau_quench( tau_quench )
            self.params['tau_quench'] = tau_quench
        self.params.update( kwargs )
        self.core.set_params( numpy.asarray( [
            self.params[k]
            for k in sfh_tunables[self.params['model']]
        ], dtype=float) )
        return;

    def Mstar ( self, tau, npoints = 100 ) :
        """ Computes the stellar mass at a given age of the galaxy.
        
        It approximates the integral:

        .. math::
        
          M_\\ast(\\tau') = \\int_0^{\\tau'}\\text{d}\\tau 
                            \\bigl[1 - \\mathcal{R}_\\text{IMF}(\\tau)\\bigr]\\psi(\\tau)

        Parameters
        ----------
        tau : float or array-like of floats
          galaxy age in years.
        npoints : int
          thinness for approximated integral computation (default is 100)

        Returns
        -------
        : float or array-like of floats
          the stellar mass of the galaxy at time :math:`\\tau`        
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
        """ Returns the dust mass at a given age of the galaxy.
  
        For empirical models of star formation (i.e. :code:`const`, 
        :code:`delayedexp`, :code:`lognorm`) this is a free parameter.

        For the :code:`insitu` model, the dust mass is given by 

        .. math::
           
           M_\\text{dust} = M_\\text{gas}(\\tau)D(\\tau)

        where 

        .. math::
           
           M_\\text{gas}=\\psi(\\tau)\\tau_\\ast
 
        and where :math:`D(\\tau)` is the gas mass ratio 
        (for an analytic expression of this quantity 
        see Pantoni et al. 2019 and Lapi et al. 2020).

        Parameters
        ----------
        tau : float or array-like of floats
          age of the galaxy in years.

        Returns
        -------
        : float or array-like of floats
          Dust content in solar masses (:math:`M_\\odot`) at give time :math:`\\tau`.
        """
        
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
        """ Returns the gas mass at a given age of the galaxy.

        For empirical models of star formation (i.e. :code:`const`, 
        :code:`delayedexp`, :code:`lognorm`) 
        this is given by 

        .. math::
           
           M_\\text{gas} = M_\\text{dust}/D
        
        where :math:`D \\sim 0.01 (Z_\\text{gas}/Z_\\odot)^{-0.85}` is the 
        dust-to-gas mass ratio, derived from observations.

        For the :code:`insitu` model, the gas mass is given by

        .. math::
           
           M_\\text{gas}=\\psi(\\tau)\\tau_\\ast

        Parameters
        ----------
        tau : float or array-like of floats
          age of the galaxy in years.

        Returns
        -------
        : float or array-like of floats 
          Gas content in solar masses (:math:`M_\\odot`) at give time :math:`\\tau`.
        """
        
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
        """ Returns the gas metallicity at a given age of the galaxy.
        
        For empirical models of star formation (i.e. :code:`const`, 
        :code:`delayedexp`, :code:`lognorm`) this is a free parameter with 
        :math:`Z_\\text{gas} = Z_\\ast`.

        For the :code:`insitu` model, it is instead given by

        .. math::

          Z_\\text{gas}=\\dfrac{s y_Z}{s\\gamma-1}
                        \\biggl[
                           1 - \\dfrac{(s\\gamma-1)x}{e^{(s\\gamma-1)x}-1}
                        \\biggr]
    
        where :math:`x\\equiv\\tau/s\\tau_\\ast` and :math:`y_Z\\approx0.04` 
        is the metal production yield (including recycling) for a Chabrier IMF.

        Parameters
        ----------
        tau : float or array-like of floats
          age of the galaxy in years.
    
        Returns
        -------
        : float or array-like of floats
          Gas absolute metallicity at give time :math:`\\tau`.
        """
    
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
        """Returns the stellar metallicity at a given age of the galaxy.

        For empirical models of star formation (i.e. :code:`const`, 
        :code:`delayedexp`, :code:`lognorm`) this is a free parameter with 
        :math:`Z_\\ast = Z_\\text{gas}`.

        For the :code:`insitu` model, it is instead given by

        .. math::

          Z_\\ast=\\dfrac{y_Z}{\\gamma-1}
                  \\biggl[1 - \\dfrac{s\\gamma}{(s\\gamma-1}
                  \\dfrac{e^{-x}-e^{-s\\gamma x}[1 + (s\\gamma -1)x]}
                  {s\\gamma -1 + e^{-s\\gamma x}- s\\gamma e^{-x}}\\biggr]

        where :math:`x\\equiv\\tau/s\\tau_\\ast` and :math:`y_Z\\approx0.04` 
        is the metal production yield (including recycling) for a Chabrier IMF.

        Parameters
        ----------
        tau : float or array-like of floats
          age of the galaxy in years.

        Returns
        -------
        : float or array-like of floats
          Stellar absolute metallicity at give time :math:`\\tau`.
        """
        
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
    
