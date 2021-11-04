# External imports
import numpy
import os 

# Internal imports
from .CSP_core import loadSSP, CCSP
import galapy.internal.globs as GP_GBL

SSP_LIB = {
    'bc03.basel.chab.extend' : os.path.join( os.path.dirname( GP_GBL.__file__ ),
                                              GP_GBL.bc03_basel_chab_zeros ),
    'bc03.stelib.chab.extend' : os.path.join( os.path.dirname( GP_GBL.__file__ ),
                                              GP_GBL.bc03_stelib_chab_zeros ),
    'bc03.stelib.chab.extrap' : os.path.join( os.path.dirname( GP_GBL.__file__ ),
                                              GP_GBL.bc03_stelib_chab_extrap ),
    }
""" Dictionary providing available SSP-libraries.

Available SSP formats:
----------------------
* bc03.basel.chab.extend :
* bc03.stelib.chab.(extend/extrap) :
"""

def print_ssp_libs () :
    list_libs = ''
    for k in SSP_LIB.keys() :
        list_libs += f'* {k};\n'
    print( 'Available SSP formats\n'
           '---------------------\n'
           + list_libs )
    return;

class CSP () :

    def __init__ ( self, ssp_lib = 'bc03.basel.chab.extend' ) :
        
        self.ssp_lib = None
        if ssp_lib in set(SSP_LIB.keys()) :
            self.ssp_lib = ssp_lib
        else :
            raise ValueError( f'SSP library "{ssp_lib}" not available; '
                              'you can see the list of available libraries by running '
                              'galapy.CompositeStellarPopulation.print_ssp_libs()' )
        self.l, self.t, self.Z, self.L = loadSSP(SSP_LIB[self.ssp_lib])
        self.core = CCSP( self.l, self.t, self.Z, self.L )
        self._timetuple = None

        # steal docstrings from C-core:
        self.SSP.__func__.__doc__      = self.core.SSP.__doc__
        self.emission.__func__.__doc__ = self.core.emission.__doc__

    def set_parameters ( self, age, sfh ) :
        from .StarFormationHistory import SFH
        from .SFH_core import CSFH
        
        if isinstance( sfh, CSFH ) :
            self._timetuple = sfh.time_grid( age, self.t, self.Z )
        elif isinstance( sfh, SFH ) :
            self._timetuple = sfh.core.time_grid( age, self.t, self.Z )
        else :
            raise ValueError( 'Argument sfh must be an instance of either '
                              'galapy.StarFormationHistory.SFH or '
                              'galapy.SFH_core.' )
        
        self.core.set_params( *self._timetuple )
        return;
            
    def SSP ( self, il, it, iz ) :
        return self.core.SSP( il, it, iz )

    def emission ( self, age, sfh, il = None, ftau = None ) :

        self.set_parameters( age, sfh )
        if not il :
            il = numpy.arange( len(self.l), dtype = numpy.uint64 )
        if ftau :
            return self.core.emission( il, ftau )
            
        return self.core.emission( il )
        
    
