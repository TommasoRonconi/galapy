""" The Cosmology module defines a Cosmology class whose just a wrapper for 
interpolation on grid of Luminosity Distances and Ages
"""

########################################################################################

# External imports
import numpy
from collections.abc import MutableMapping as MM

# Internal imports
import galapy.internal.globs as GP_GBL
from galapy.internal.constants import Lsun, clight, Mpc_to_cm
from galapy.internal.interp import lin_interp
from galapy.internal.data import DataFile

########################################################################################

class CSM () :
    r""" Class for cosmological computations. It is built using pre-computed redshift-dependent
    quantities useful for transforming energies radiated per unit wavelength into fluxes from
    objects at a given distance.
    
    Parameters
    ----------
    cosmo : string or dictionary
      If a string is passed, it should name one of the pre-computed cosmologies available in the
      database. Available cosmologies are :code:`('WMAP7', 'WMAP9', 'Planck15', 'Planck18')`.
      If a dictionary is passed, the class expects to find 3 key-value couples:
      (i) key = 'redshift', value = an array of redshift values; 
      (ii) key = 'luminosity_distance', value = an array of luminosity distances corresponding 
      to the redshift values of the first key-value couple;
      (iii) key = 'age', value = an array of ages of the Universe corresponding to the redshift 
      values of the first key-value couple.
      All these arrays should have the same length.
    
    Attributes
    ----------
    DL : :py:class:`galapy.internal.interp.lin_interp`
       interpolator interface to luminosity distance
    UA : :py:class:`galapy.internal.interp.lin_interp`
       interpolator interface to age of the universe

    Examples
    --------
    Interpolator attributes can be called as functions that take as input a redshift value 
    (or an array of redshift values):

    >>> from galapy.Cosmology import CSM
    >>> csm = CSM('Planck18')
    >>> csm.DL(1.0)
    6791.26894

    The returned value is in MegaParsecs for the luminosity distance and in years for the age.
    """
    
    def __init__ ( self, cosmo ) :
        
        if isinstance( cosmo, str ) :
            zz, DL, UA = numpy.loadtxt( DataFile( f'{cosmo:s}.txt',
                                                  GP_GBL.COSMO_DIR ).get_file(),
                                        unpack = True )
        elif isinstance( cosmo, MM ) :
            zz, DL, UA = cosmo['redshift'], cosmo['luminosity_distance'], cosmo['age']
        else :
            raise ValueError( 'Argument `cosmo` should be either a string '
                              'or a formatted dictionary.' )
        zz = numpy.ascontiguousarray(zz)
        DL = numpy.ascontiguousarray(DL)
        TF = numpy.ascontiguousarray(
            1.e+26 * (1 + zz) * Lsun /
            ( 4 * numpy.pi * clight['A/s'] * DL**2 * Mpc_to_cm**2 )
        )
        UA = numpy.ascontiguousarray(UA)
        self.DL  = lin_interp( zz, DL )
        self.age = lin_interp( zz, UA )
        self._TF = lin_interp( zz, TF )

    def to_flux ( self, redshift, restframe_wavelength, luminosity )  :
        r"""Converts a restframe luminosity to the flux received at 
        a given redshift.
        
        .. math::
        
          S_{\lambda_O} = \lambda_R^2 L_\text{tot}(\lambda_R) (1+z) / (4 \pi c D_L^2 )
        
        Parameters
        ----------
        redshift : float
           redshift of the source
        restframe_wavelength : array-like 
           rest-frame wavelength grid
        luminosity : array-like
           rest-frame luminosity grid
        
        Returns
        -------
        : array-like
          a flux in milliJansky
        """
        return luminosity * restframe_wavelength**2 * self._TF( redshift ) 

########################################################################################

