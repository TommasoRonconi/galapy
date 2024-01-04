""" The Inter-Stellar Medium module implements and wraps the absorbing and re-radiating media of the Galaxy.
"""

# External Imports
import numpy

# Internal imports
from .ISM_core import CMC, CDD, total_attenuation

ism_tunables = {
    # Molecular Clouds
    'mc' : [ 'f_MC', 'norm_MC', 'N_MC', 'R_MC', 'Zgas', 'tau_esc', 'Mgas', 'dMClow', 'dMCupp' ],
    
    # Diffuse Dust
    'dd' : [ 'f_MC', 'norm_DD', 'Mdust', 'Rdust', 'f_PAH', 'dDDlow', 'dDDupp' ],
}
""" Dictionary with the tunable parameters of the ISM phases

#. :code:`mc` : **Molecular Clouds**
#. :code:`dd` : **Diffuse Dust**

"""

def ism_build_params ( phase, **kwargs ) :
    """ Builds the parameters dictionary for given phase.
    """
    if phase not in ism_tunables.keys() :
        av_phases = ''
        for k in ism_tunables.keys() :
            av_phases += f'"{k}" '
        raise ValueError( f'Phase "{phase}" chosen not available. '
                          f'Available phases are: {av_phases}')
    if phase == 'mc' :
        out = {
            'f_MC'    : 0.5,
            'norm_MC' : 1.e+02,
            'N_MC'    : 1.e+03,
            'R_MC'    : 10.,
            'Zgas'    : 0.01,
            'tau_esc' : 1.e+07,
            'Mgas'    : 1.e+09,
            'dMClow'  : 1.3,
            'dMCupp'  : 1.6,
        }
        
    if phase == 'dd' :
        out = {
            'f_MC'    : 0.5,
            'norm_DD' : 1.e+00,
            'Mdust'   : 1.e+07,
            'Rdust'   : 1.e+03,
            'f_PAH'   : 0.2,
            'dDDlow'  : 0.7,
            'dDDupp'  : 2.0,
        }
        
    for k in set( out.keys() ).intersection(kwargs.keys()) :
        out[k] = kwargs[k]
    return out
    
class ismPhase ():
    """ ISM phase base class.
    """

    def __init__ ( self, phase, builder, T = None, **kwargs ) :
        
        self.phase  = phase
        self.params = ism_build_params( self.phase, **kwargs )
        self.core   = builder()
        self.set_parameters( **kwargs )
        self.T = None
        if T :
            self.set_temperature( T )
                
        # steal docstrings from C-core:
        self.set_temperature.__func__.__doc__ = self.core.set_temperature.__doc__
        self.temperature.__func__.__doc__     = self.core.temperature.__doc__
        self.emission.__func__.__doc__        = self.core.emission.__doc__
        self.attenuation.__func__.__doc__     = self.core.attenuation.__doc__
        self.extinction.__func__.__doc__      = self.core.extinction.__doc__
        self.A_V.__func__.__doc__             = self.core.A_V.__doc__

    def set_parameters ( self, **kwargs ) :
        """ Function for setting the free parameters of the class.

        It is implemented to take any number of keyword arguments and
        automathically ignores arguments that are not free-parameters of the class.
        """
        self.params.update( kwargs )
        self.core.set_params( numpy.asarray( [
            self.params[k]
            for k in ism_tunables[self.phase]
        ], dtype = float ) )
    
    def set_temperature ( self, T ) :
        """ Manually set the temperature of the ISM phase
        
        Parameters
        ----------
        T : float
          Temperature in Kelvin :math:`[K]`
        """
        self.T = T
        self.core.set_temperature( self.T )
        return;
    
    def set_slopes ( self, lower, upper ) :
        r""" Manually set the slopes of the ISM phase extinction
        
        Parameters
        ----------
        lower : float
          slope of extinction at wavelengths <= 100 :math:`\mu m`
        lower : float
          slope of extinction at wavelengths > 100 :math:`\mu m`
          
        """
        self.params[ism_tunables[self.phase][-2]] = lower 
        self.params[ism_tunables[self.phase][-1]] = upper
        self.core.set_slopes( lower, upper )
        return;

    def temperature ( self, Etot ) :
        r""" Computes the temperature of the ISM at given total energy.

        With the assumption that the total energy density :math:`E_\text{abs}^\text{phase}`
        absorbed by the ISM phase is re-radiated as a grey body 
        with luminosity :math:`L_\lambda[\tau, T]`,
        the ISM temperature is computed by solving the integral
        
        .. math:: 
          
          \int_0^\infty \text{d}\lambda L_\lambda^\text{phase}[\tau|T_\text{phase}(\tau)] =
          E_\text{abs}^\text{phase}(\tau)
        
        Parameters
        ----------
        Etot : scalar float
          Total energy absorbed by the ISM phase.

        Returns
        -------
        T : scalar float
          Temperature in Kelvin :math:`[K]` of the given ISM phase.
        """
        self.T = self.core.temperature( Etot )
        return self.T

    def emission ( self, wavelength, T = None ) :
        r""" Computes the ISM emission at given wavelength.
        
        We assume the ISM radiates as a gray body with emission spectrum
        
        .. math::
        
          L_\lambda^\text{phase}(\tau) = 
          \mathcal{N}_\text{phase}\; \bigl[1 
          - 10^{-0.4\,A_\lambda^\text{phase}(\tau)}\bigr]\;
          B_\lambda(T_\text{phase})
        
        where :math:`\mathcal{N}_\text{phase}` is a normalization depending
        on the model parameters, 
        the factor :math:`1 - 10^{-0.4\,A_\lambda^\text{phase}(\tau)}` 
        is the optical depth of the ISM phase and :math:`B_\lambda(T_\text{phase})`
        is the temperature of the medium.          

        Parameters
        ----------
        wavelength : array or scala float
          wavelength in angstroms :math:`[\mathring{A}]`

        Keyword Arguments
        -----------------
        T : float
          temperature in Kelvin :math:`[K]` of the ISM phase 
          (optional, default = :code:`None`)

        Returns
        -------
        L_ISM : float of array-like of floats
          Luminosity at given wavelength 
        """
        if T :
            self.set_temperature( T )
        return self.core.emission( wavelength )

    def attenuation ( self, wavelength ) :
        """ Computes the ISM attenuation at given wavelength.

        Parameters
        ----------
        wavelength : array or scalar float
          wavelength in angstroms :math:`[\mathring{A}]`
        
        Returns
        -------
         : array or scalar float
        """
        return self.core.attenuation( wavelength )

    def extinction ( self, wavelength ) :
        """ Computes the ISM extinction at given wavelength.

        Parameters
        ----------
        wavelength : array or scalar float
          wavelength in angstroms :math:`[\mathring{A}]`
        
        Returns
        -------
         : array or scalar float
        """
        return self.core.extinction( wavelength )

    def A_V ( self ) :
        """ Returns the extinction value in the visible band.
        """
        return self.core.A_V()

class MC ( ismPhase ) :
    """Class implementing the Molecular-Cloud medium.
    
    Parameters
    ----------
    T : float
       (optional, default = None) average temperature of Molecular Clouds 
    kwargs : dictionary
       (optional) can contain any parameter value (key-value pairs)
    """

    def __init__ ( self, T = None, **kwargs ) :
        super().__init__( 'mc', CMC, T=T, **kwargs )

    def eta ( self, tt ) :
        r"""Implements the function regulating the time-dependent
        evaporation of Molecular clouds:

        .. math::
            
            \eta(\tau) 
            = \begin{cases}
                1  & \tau\leq \tau_\text{esc}\\
                2-\frac{\tau}{\tau_\text{esc}} & \tau_\text{esc}<\tau\leq 2\,\tau_\text{esc}\\
                0 & \tau>2\,\tau_\text{esc}
            \end{cases}

        where :math:`\tau_\text{esc}` is a free parameter of the MC class and :math:`\tau` is the
        age at which the hosting galaxy is being observed.

        Parameters
        ----------
        tt : array or scalar float
            age of the hosting galaxy (in years)
        
        Returns
        -------
        : array or scalar float
           value of the :math:`\eta(\tau)` parameter at the input age(s)
        """
        return self.core.eta( tt )

    def time_attenuation ( self, wavelength, tt ) :
        r"""Computes the time-dependent attenuation due to Molecular Clouds:
        
        .. math::
            
            \mathcal{A}_\text{MC}(\lambda) = 1-\eta(\tau)\left[ 1 - 10^{-0.4\, A_\text{MC}(\lambda)} \right]
        
        where :math:`\eta(\tau)` is computed with :code:`MC.eta()`.
        
        Parameters
        ----------
        wavelength : array or scalar float
            wavelength in angstroms :math:`[\mathring{A}]`
        tt : array or scalar float
            age of the hosting galaxy (in years)
        
        Returns
        -------
        : float or iterable
            Output value of the time-dependent attenutation from MCs. The output shape depends
            on the inputs and is = :code:`(len(wavelength), len(tt))`
        """
        return (
            1 - ( 1 - self.attenuation( wavelength ) )[:,numpy.newaxis]
            * self.eta( tt )[numpy.newaxis,:]
        )
        
class DD ( ismPhase ) :
    """Class implementing the Diffuse-Dust medium.
    
    Parameters
    ----------
    T : float
       (optional, default = None) average temperature of Diffuse Dust 
    kwargs : dictionary
       (optional) can contain any parameter value (key-value pairs)
    """

    def __init__ ( self, T = None, **kwargs ) :
        super().__init__( 'dd', CDD, T=T, **kwargs )

class ISM () :
    """Class implementing the Inter-Stellar Medium.
    Wraps up the combined effect of Molecular Clouds and Diffuse Dust.
    
    Parameters
    ----------
    TMC : float
       (optional, default = None) average temperature of Molecular Clouds 
    TDD : float
       (optional, default = None) average temperature of Diffuse Dust 
    kwargs : dictionary
       (optional) can contain any parameter value (key-value pairs) of both
       MCs and DD
    """

    def __init__ ( self, TMC = None, TDD = None, **kwargs ) :
        self.mc = MC( TMC,
                      **{ k : v
                          for k, v in kwargs.items()
                          if k in ism_tunables[ 'mc' ] } )
        self.dd = DD( TDD,
                      **{ k : v
                          for k, v in kwargs.items()
                          if k in ism_tunables[ 'dd' ] } )

    def set_parameters ( self, **kwargs ) :
        """Change the value of the free parameters.
        
        Keyword Arguments
        -----------------
        f_MC    : float
          Fraction of the total dust mass in MCs 
        norm_MC : float
          Normalization of the attenuation law of MCs
        N_MC    : float
          Number of MCs within the hosting galaxy
        R_MC    : float
          Average radius of a MC
        Zgas    : float
          Average metallicity of the gas in the hosting galaxy
        tau_esc : float
          Average escape time of SSPs from MCs
        Mgas    : float
          Total gas mass
        dMClow  : float
          spectral index of the MC extinction at short wavelengths
        dMCupp  : float
          spectral index of the MC extinction at long wavelengths
        norm_DD : float
          Normalization of the attenuation law of DD
        Mdust   : float
          Total dust mass
        Rdust   : float
          Radius of the DD region
        f_PAH   : float
          fraction of the total emission from DD contributed by PAH
        dDDlow  : float
          spectral index of the DD extinction at short wavelengths
        dDDupp  : float
          spectral index of the DD extinction at long wavelengths
        """
        self.mc.set_parameters( **{ k : v
                                    for k, v in kwargs.items()
                                    if k in ism_tunables[ 'mc' ] } )
        self.dd.set_parameters( **{ k : v
                                    for k, v in kwargs.items()
                                    if k in ism_tunables[ 'dd' ] } )
        return;

    def total_attenuation ( self, wavelength, tt ) :
        """Function computing the total attenuation due to MCs and DD.
        
        Note that the output is flattened and the total size depends on the
        size of the input arrays.
        Both output arrays (e.g. :code:`arr`) can be reshaped to their original dimensions by calling
        
        .. code::
            
            arr.reshape( wavelength.size, tt.size )
        
        for the case when both :code:`wavelength` and :code:`tt` are arrays.
        If one of them is not, flattening the output does not have any effect.
        
        Parameters
        ----------
        wavelength : array or scalar float
          wavelength in angstroms :math:`[\mathring{A}]`
        tt : array or scalar float
          age of the hosting galaxy (in years)
        
        Returns
        -------
        : 1d array
          Total time attenuation due to MCs for each wavelength and age in the input grid.
          (uses function :code:`MC.time_attenuation` and flattens the output array.
        : 1d array
          Total time attenuation due to the combined effect of MCs first and, subsequently, DD.
          (the output array is flattened)
        """
        attMC = self.mc.time_attenuation( wavelength, tt )
        attDD = attMC * self.dd.attenuation( wavelength )[:,numpy.newaxis]
        return ( numpy.ascontiguousarray( attMC.ravel() ),
                 numpy.ascontiguousarray( attDD.ravel() ) )

