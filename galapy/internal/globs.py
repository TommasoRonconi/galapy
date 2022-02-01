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

# BaSeL (Low Resolution) SSP + Chabrier IMF + null ext. + thinner lambda-grid
bc03_basel_chab_zeros_refined  = os.path.join(
    BC03, 'basel_chab_zero_refined.dat'
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

FILT_DIR = os.path.join( DATA_DIR, 'filters' )

################################################################
# Luminosity Distances

DL_DIR = os.path.join( DATA_DIR, 'DL' )

################################################################
# AGN Templates (Fritz et al., 2006)

AGN_FILE = os.path.join( DATA_DIR, 'AGN_templates',
                         'ct{0:d}al{1:.1f}be{2:.2f}ta{3:.1f}rm{4:d}ia{5:.3f}' )

################################################################
