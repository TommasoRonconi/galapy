""" The Inter-Stellar Medium module implements and wraps the absorbing media of the Galaxy.
"""

# External Imports
import numpy

# Internal imports
from .ISM_core import CMC, CDD, total_attenuation

ism_tunables = {
    # Molecular Clouds
    'mc' : [ 'f_MC', 'norm_MC', 'N_MC', 'R_MC', 'Zgas', 'tau_esc', 'Mgas' ],
    
    # Diffuse Dust
    'dd' : [ 'f_MC', 'norm_DD', 'Mdust', 'Rdust', 'f_PAH' ],
}
""" Dictionary with the tunable parameters of the ISM phases

#. :code:`mc` : **Molecular Clouds**
#. :code:`dd` : **Diffuse Dust**

"""

def ism_build_params ( phase, **kwargs ) :
    """ Builds the parameters dictionary for given phase.
    """
    if phase == 'mc' :
        out = {
            'f_MC'    : 0.5,
            'norm_MC' : 1.e+02,
            'N_MC'    : 1.e+03,
            'R_MC'    : 10.,
            'Zgas'    : 0.5,
            'tau_esc' : 1.e+07,
            'Mgas'    : 1.e+09,
        }
        for k in set( out.keys() ).intersection(kwargs.keys()) :
            out[k] = kwargs[k]
        return out
    if phase == 'dd' :
        return {
            'f_MC'    : 0.5,
            'norm_DD' : 1.e+00,
            'Mdust'   : 1.e+07,
            'Rdust'   : 1.e+03,
            'f_PAH'   : 0.2
        }
        for k in set( out.keys() ).intersection(kwargs.keys()) :
            out[k] = kwargs[k]
        return out
    av_phases = ''
    for k in ism_tunables.keys() :
        av_phases += f'"{k}" '
    raise ValueError( f'Phase "{phase}" chosen not available. '
                      f'Available phases are: {av_phases}')
    
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
        """ Function for ...
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

    def emission ( self, wavelenght, T = None ) :
        r""" Computes the ISM emission at given wavelenght.
        
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
        wavelenght : float of array-like of floats
          wavelenght in angstroms :math:`[\mathring{A}]`

        Keyword Arguments
        -----------------
        T : float
          temperature in Kelvin :math:`[K]` of the ISM phase 
          (optional, default = :code:`None`)

        Returns
        -------
        L_ISM : float of array-like of floats
          Luminosity at given wavelenght 
        """
        if T :
            self.set_temperature( T )
        return self.core.emission( wavelenght )

    def attenuation ( self, wavelenght ) :
        """ Computes the ISM attenuation at given wavelenght.

        Parameters
        ----------
        wavelenght : array or scalar float
        
        Returns
        -------
         : array or scalar float
        """
        return self.core.attenuation( wavelenght )

    def extinction ( self, wavelenght ) :
        """ Computes the ISM extinction at given wavelenght.

        Parameters
        ----------
        wavelenght : array or scalar float
        
        Returns
        -------
         : array or scalar float
        """
        return self.core.extinction( wavelenght )

    def A_V ( self ) :
        """ Returns the extinction value in the visible band-
        """
        return self.core.A_V()

class MC ( ismPhase ) :

    def __init__ ( self, T = None, **kwargs ) :
        super().__init__( 'mc', CMC, T )

    def eta ( self, tt ) :
        return self.core.eta( tt )

    def time_attenuation ( self, wavelenght, tt ) :
        return (
            1 - ( 1 - self.attenuation( wavelenght ) )[:,numpy.newaxis]
            * self.eta( tt )[numpy.newaxis,:]
        )
        
class DD ( ismPhase ) :

    def __init__ ( self, T = None, **kwargs ) :
        super().__init__( 'dd', CDD, T )

class ISM () :

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
        self.mc.set_parameters( **{ k : v
                                    for k, v in kwargs.items()
                                    if k in ism_tunables[ 'mc' ] } )
        self.dd.set_parameters( **{ k : v
                                    for k, v in kwargs.items()
                                    if k in ism_tunables[ 'dd' ] } )
        return;

    def total_attenuation ( self, wavelenght, tt ) :
        attMC = self.mc.time_attenuation( wavelenght, tt )
        attDD = attMC * self.dd.attenuation( wavelenght )[:,numpy.newaxis]
        return ( numpy.ascontiguousarray( attMC.ravel() ),
                 numpy.ascontiguousarray( attDD.ravel() ) )

