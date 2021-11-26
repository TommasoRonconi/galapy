import os

################################################################
#                      DATA DIRECTORY:                         #
# contains binary and ASCII files to be used along the library #
################################################################

DATA_DIR = 'data'

################################################################
# Bruzual and Charlot (2003) SSP libraries

# directory
BC03     = os.path.join( DATA_DIR, 'BC03' )

# BaSeL (Low Resolution) SSP + Chabrier IMF + null extension
bc03_basel_chab_zeros  = os.path.join(
    BC03, 'basel_chab_zero_extended.dat'
)

# Stelib (High Resolution) SSP + Chabrier IMF + null extension
bc03_stelib_chab_zeros  = os.path.join(
    BC03, 'stelib_chab_zero_extended.dat'
)

# Stelib (High Resolution) SSP + Chabrier IMF + extrapolation
bc03_stelib_chab_extrap = os.path.join(
    BC03, 'stelib_chab_extrap.dat'
)

################################################################
# Photometric Filters

FILT = os.path.join( DATA_DIR, 'filters' )

################################################################
