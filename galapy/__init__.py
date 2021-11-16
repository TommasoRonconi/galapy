""" This package simulates galactic spectra
"""

import galapy.StarFormationHistory
import galapy.CompositeStellarPopulation
import galapy.InterStellarMedium

__all__ = [ 'StarFormationHistory', 'CompositeStellarPopulation', 'InterStellarMedium' ]

# import galapy.internal.globs as GP_GBL
# import os 

# SSP_LIB = {
#     'bc03.basel.chab.extend' : os.path.join( os.path.dirname( GP_GBL.__file__ ),
#                                               GP_GBL.bc03_basel_chab_zeros ),
#     'bc03.stelib.chab.extend' : os.path.join( os.path.dirname( GP_GBL.__file__ ),
#                                               GP_GBL.bc03_stelib_chab_zeros ),
#     'bc03.stelib.chab.extrap' : os.path.join( os.path.dirname( GP_GBL.__file__ ),
#                                               GP_GBL.bc03_stelib_chab_extrap ),
#     }
# """ Dictionary providing available SSP-libraries.

# Available SSP formats:
# ----------------------
# * bc03.basel.chab.extend :
# * bc03.stelib.chab.(extend/extrap) :
# """

