""" The Cosmology module defines a Cosmology class whose just a wrapper for interpolation on grid of Luminosity Distances and Ages
"""

# External imports
import numpy

# Internal imports
import galapy.internal.globs as GP_GBL
from galapy.internal.constants import Lsun, clight, Mpc_to_cm
from galapy.internal.interp import lin_interp
from galapy.internal.data import DataFile


class CSM () :
    
    def __init__ ( self, cosmo ) :
        """"""
        
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

    def to_flux ( self, redshift, restframe_wavelenght, luminosity )  :
        """ Converts a restframe luminosity to the flux received at 
        a given redshift.
        
        .. math::
        
          S_{\lambda_O} = \lambda_R^2 * L_tot(\lambda_R) * (1+z)/(4 * \pi * c * D_L^2 )
        
        Parameters
        ----------
        redshift : float
        restframe_wavelenght : array-like 
        luminosity : array-like
        
        Returns
        -------
        : array-like
          a flux in milliJansky
        """
        return luminosity * restframe_wavelenght**2 * self._TF( redshift ) 

        


