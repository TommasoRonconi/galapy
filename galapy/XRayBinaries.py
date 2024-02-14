"""The XRayBinaries module implements the contribution to the X-band emission from stellar binary systems.

This module provides classes for modeling X-ray emission from both High Mass Stars (HMXRB) and Low Mass Stars (LMXRB). 
"""

########################################################################################

# External imports
from abc import ABC, abstractmethod
import numpy

# Internal imports
from galapy.internal.utils import trap_int, powerlaw_exp_cutoff, poly_N
from galapy.internal.constants import Ang_to_keV, Lsun, sunL
from galapy.internal.interp import lin_interp

########################################################################################

xrb_tunables = {
    # High Mass Stars
    'hm' : [ 'psi', 'Zstar' ],
    
    # Low Mass Stars
    'lm' : [ 'Mstar', 'age' ],
}
""" Dictionary with the tunable parameters of the X-ray contribution from binary star systems

#. :code:`hm` : **High Mass Star**
#. :code:`lm` : **Low Mass Stars**

"""

def xrb_build_params ( startype, **kwargs ) :
    """Builds the parameters dictionary for a given star type.
    
    Parameters
    ----------
    startype : str
        Type of star ('hm' for High Mass Star, 'lm' for Low Mass Star).
    **kwargs : dict
        Additional keyword arguments for specifying parameters.
        The :py:data:`galapy.XRayBinaries.xrb_tunables` contains the 
        available free-parameters for the 2 different types of Xray Binaries:
        high and low mass stars.

    Returns
    -------
    dict
        Dictionary containing the parameters for the specified star type.
    """
    if startype not in xrb_tunables.keys() :
        av_types = ''
        for k in xrb_tunables.keys() :
            av_types += f'"{k}" '
        raise ValueError( f'Star type "{startype}" chosen not available. '
                          f'Available phases are: {av_types}')
    
    if startype == 'hm' :
        out = {
            'psi'   : 1.,
            'Zstar' : 0.02,
        }
        
    if startype == 'lm' :
        out = {
            'Mstar' : 1.e+10,
            'age'   : 1.e+8,
        }

    for k in set( out.keys() ).intersection(kwargs.keys()) :
        out[k] = kwargs[k]
    return out

########################################################################################

class _XRBase ( ABC ) :
    r"""Abstract base class for X-ray emission models.

    Implements an input X-ray luminosity on a power-law spectrum with exponential cut-off

    .. math::
    
       L_\Gamma^X(\lambda)\propto E^{-\Gamma+3}(\lambda)\, e^{-E(\lambda)/E_\text{cut}}

    with :math:`E(\lambda)=h_P\, \nu_\lambda=h_P\, c/\lambda` the energy of a photon 
    with wavelength :math:`\lambda`, :math:`E_\text{cut}` the characteristic energy
    of the exponential cut-off and :math:`\Gamma` the photon index.
        
    Parameters
    ----------
    startype : str
        Type of star ('hm' for High Mass Star, 'lm' for Low Mass Star).
    gamma : float
        Power-law index for the X-ray emission.
    Ecut : float
        Cutoff energy for the X-ray emission.
    lmin : float, optional
        Minimum wavelength, default is 1.
    lmax : float, optional
        Maximum wavelength, default is 1.e+10.
    **kwargs : dict
        Additional keyword arguments for specifying parameters.
    """

    @property
    @abstractmethod
    def coeff ( self ) :
        """Coefficients for normalization."""
        pass

    def __init__ ( self, startype, gamma, Ecut,
                   lmin = 1., lmax = 1.e+10, **kwargs ) :
        """Initialize the X-ray base model.
        
        .. warning::
           Being an abstract class, this will give an error if the
           abstract methods are not defined.
        """
        self.startype = startype
        self.params   = xrb_build_params( self.startype, **kwargs )
        self.lmin     = lmin
        self.lmax     = lmax
        self.NormStarType = None

        ###########################################################
        # Compute X-template

        # generate wavelength grid
        ll = numpy.logspace( numpy.log10( self.lmin ),
                             numpy.log10( self.lmax ),
                             256 )
        
        # convert wavelength to energy
        El = Ang_to_keV( ll )
        
        # find interval for hard-X normalization
        wE = ( 2. <= El ) & ( El <= 10. )

        # compute emission law normalization
        Lnorm = -1. / trap_int( El[wE], powerlaw_exp_cutoff( El[wE], gamma=gamma, Ecut=Ecut ) )

        # Compute the normalized power-law:
        # - high energy cut-off at 300 keV
        # - low energy cut-off at 50 Angstrom (~0.25 keV)
        ret = numpy.zeros_like( ll )
        wL = ( ll <= 5.e+1 ) 
        ret[wL] = powerlaw_exp_cutoff( El[wL], gamma=gamma, Ecut=Ecut ) * Lnorm

        # store normalized interpolator object for X-ray emission
        self.f_norm_X = lin_interp( numpy.ascontiguousarray( ll ),
                                    numpy.ascontiguousarray( ret ) )

    @abstractmethod
    def set_parameters ( self ) :
        """Set parameters for the X-ray emission model."""
        pass
    
    def _set_norm ( self, xx, coeff, fact ) :
        """Set normalization for the X-ray emission model."""
        logNorm = poly_N( xx, coeff )
        self.NormStarType = fact * 10.**logNorm * sunL
        return;

    def emission ( self, ll ) :
        r"""Compute X-ray emission for the given wavelengths.

        Given as :math:`\mathcal{N}_\text{type} \cdot L_\Gamma^X(\lambda)`
        where the normalisation factor depends on the stellar type.
        
        Parameters
        ----------
        ll : array-like
            Wavelengths at which to compute the X-ray emission.

        Returns
        -------
        array-like
            X-ray emission values corresponding to the input wavelengths.
        """
        return self.NormStarType * self.f_norm_X( ll )

    
########################################################################################

class HMXRB ( _XRBase ) :
    r"""Class for handling X-ray emission from High Mass Stars.

    Distributes the X-ray luminosity on a power-law spectrum with exponential cut-off

    .. math::
    
       L_\Gamma^X(\lambda)\propto E^{-\Gamma+3}(\lambda)\, e^{-E(\lambda)/E_\text{cut}}

    with :math:`E(\lambda)=h_P\, \nu_\lambda=h_P\, c/\lambda` the energy of a photon 
    with wavelength :math:`\lambda`, :math:`E_\text{cut}` the characteristic energy
    of the exponential cut-off and :math:`\Gamma = 2` 
    (`Fabbiano, 2006 <https://doi.org/10.1146/annurev.astro.44.051905.092519>`_) 
    the photon index.

    The normalisation factor of the exponential cut-off emission is given by

    .. math::
    
       \log(L_\text{HMXB}/\text{erg}\,s^{-1}) \approx\\ \log({\dot M_\star/M_\odot}\text{yr}^{-1}) + 
                            40.28-62.12\, Z_\star\\+569.44\,Z_\star^2-1883.80\,Z_\star^3+1968.33\, Z_\star^4

    where :math:`{\dot M_\star} \equiv \psi(\tau)` is the star formation rate and 
    :math:`Z_\star` is the stellar absolute metallicity.
    
    The emission is therefore computed as
    
    .. math::
       
       L(\lambda, Z_\star) = L_{\Gamma=2}^X(\lambda) \cdot L_\text{HMXB}
        
    Parameters
    ----------
    lmin : float
        Minimum wavelength for the X-ray emission.
    lmax : float
        Maximum wavelength for the X-ray emission.
    **kwargs : dict
        Additional keyword arguments for specifying parameters.
    
    Keyword Arguments
    -----------------
    psi : float
       star formation rate, default is 1.
    Zstar : float
       stellar average absolute metallicity, default is 0.02
    """

    coeff = [ 1968.33, -1883.80, 569.44, -62.12, 40.28 ]

    def __init__ ( self, lmin, lmax, **kwargs ) :
        """Initialize HMXRB object."""
        super().__init__( startype = 'hm',
                          gamma    = 2.0,
                          Ecut     = 100.,
                          lmin = lmin, lmax = lmax,
                          **kwargs )
        self.set_parameters()

    def set_parameters ( self, **kwargs ) :
        """Set parameters for HMXRB model.
    
        Keyword Arguments
        -----------------
        psi : float
           star formation rate.
        Zstar : float
           stellar average absolute metallicity.
        """

        self.params.update( kwargs )
        super()._set_norm( self.params['Zstar'],
                           type( self ).coeff,
                           self.params['psi'] )
        return;

    
########################################################################################

class LMXRB ( _XRBase ) :
    r"""Class for handling X-ray emission from Low Mass Stars.

    Distributes the X-ray luminosity on a power-law spectrum with exponential cut-off

    .. math::
    
       L_\Gamma^X(\lambda)\propto E^{-\Gamma+3}(\lambda)\, e^{-E(\lambda)/E_\text{cut}}

    with :math:`E(\lambda)=h_P\, \nu_\lambda=h_P\, c/\lambda` the energy of a photon 
    with wavelength :math:`\lambda`, :math:`E_\text{cut}` the characteristic energy
    of the exponential cut-off and :math:`\Gamma = 1.6` 
    (`Fabbiano, 2006 <https://doi.org/10.1146/annurev.astro.44.051905.092519>`_) 
    the photon index.

    The normalisation factor of the exponential cut-off emission is given by

    .. math::
    
       \log(L_\text{LMXB}/\text{erg}\,s^{-1}) \approx\\ \log(M_\star/M_\odot) + 40.276-1.503\, 
                                \theta-0.423\,\theta^2+0.425\,\theta^3+0.136\, \theta^4

    where :math:`M_\star` is the stellar mass and :math:`\theta \equiv \log(\tau/\text{Gyr})`
    is the logarithm of the stellar population age.
    
    Assuming :math:`\Gamma = 1.6` 
    (`Fabbiano, 2006 <https://doi.org/10.1146/annurev.astro.44.051905.092519>`_) 
    the emission is computed as
    
    .. math::
       
       L(\lambda, Z_\star) = L_{\Gamma=1.6}^X(\lambda) \cdot L_\text{LMXB}
        
    Parameters
    ----------
    lmin : float
        Minimum wavelength for the X-ray emission.
    lmax : float
        Maximum wavelength for the X-ray emission.
    **kwargs : dict
        Additional keyword arguments for specifying parameters.
    
    Keyword Arguments
    -----------------
    age : float
        age of the stellar populations, default is 1.e+8
    Mstar : float
        total stellar mass of the stellar populations, default is 1.e+10
    """

    coeff = [ 0.136, 0.425, -0.423, -1.503, 40.276 ]

    def __init__ ( self, lmin, lmax, **kwargs ) :
        """Initialize LMXRB object."""
        super().__init__( startype = 'lm',
                          gamma    = 1.6,
                          Ecut     = 100.,
                          lmin = lmin, lmax = lmax,
                          **kwargs )
        self.set_parameters()

    def set_parameters ( self, **kwargs ) :
        """Set parameters for LMXRB model.
    
        Keyword Arguments
        -----------------
        age : float
            age of the stellar populations
        Mstar : float
            total stellar mass of the stellar populations
        """
        
        self.params.update( kwargs )
        super()._set_norm( numpy.log10( self.params['age'] * 1.e-9 ),
                           type( self ).coeff,
                           self.params['Mstar'] * 1.e-10 )
        return;

########################################################################################

class XRB () :
    r"""Class for handling X-ray emission from both High and Low Mass Stars.

    The total emission is given as
    
    .. math::
       
       L_\text{XRB}(\lambda, Z_\star, \tau) = L_\text{HMXB}(Z_\star) L_{\Gamma=2}^X(\lambda) + 
                                              L_\text{LMXB}(\tau) L_{\Gamma=1.6}^X(\lambda)
    
    where the first term gives the contribution from High Mass binaries and
    the second term that from Low Mass binaries.
        
    Parameters
    ----------
    lmin : float
        Minimum wavelength for the X-ray emission.
    lmax : float
        Maximum wavelength for the X-ray emission.
    **kwargs : dict
        Additional keyword arguments for specifying parameters.
    
    Keyword Arguments
    -----------------
    psi : float
       star formation rate, default is 1.
    Zstar : float
       stellar average absolute metallicity, default is 0.02
    age : float
        age of the stellar populations, default is 1.e+8
    Mstar : float
        total stellar mass of the stellar populations, default is 1.e+10
    """

    def __init__ ( self, lmin, lmax, **kwargs ) :

        self.params = {}
        self.hm = HMXRB(
            lmin, lmax,
            **{ k : v
                for k, v in kwargs.items()
                if k in xrb_tunables[ 'hm' ] }
        )
        self.params.update( self.hm.params )
        
        self.lm = LMXRB(
            lmin, lmax,
            **{ k : v
                for k, v in kwargs.items()
                if k in xrb_tunables[ 'lm' ] }
        )
        self.params.update( self.lm.params )

    def set_parameters ( self, **kwargs ) :
        """Set parameters of the model.
    
        Keyword Arguments
        -----------------
        psi : float
           star formation rate.
        Zstar : float
           stellar average absolute metallicity.
        age : float
            age of the stellar populations
        Mstar : float
            total stellar mass of the stellar populations
        """
        
        self.params.update( kwargs )
        self.hm.set_parameters( **{ k : v
                                    for k, v in kwargs.items()
                                    if k in xrb_tunables[ 'hm' ] } )
        self.lm.set_parameters( **{ k : v
                                    for k, v in kwargs.items()
                                    if k in xrb_tunables[ 'lm' ] } )
        return;

    def emission ( self, ll ) :
        """Compute X-ray emission due to X-Ray Binaries for the given wavelengths.
        
        Parameters
        ----------
        ll : array-like
            Wavelengths at which to compute the X-ray emission.

        Returns
        -------
        array-like
            X-ray emission values corresponding to the input wavelengths.
        """
        return (
            self.hm.emission( ll ) +
            self.lm.emission( ll )
        )
    
########################################################################################

