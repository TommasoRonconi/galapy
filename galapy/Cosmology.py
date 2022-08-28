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
        self.DL = lin_interp( zz, DL )
        self.age = lin_interp( zz, UA )


