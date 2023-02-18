import os
# from galapy.configuration import rcParams

################################################################
#                      DATA DIRECTORY:                         #
# contains binary and ASCII files to be used along the library #
################################################################

DATABASE = 'https://api.github.com/repos/TommasoRonconi/galapy_data/releases/'
DATA_VERSION = '0.0.1'
DATA_URL = 'https://github.com/TommasoRonconi/galapy_data/blob/main/{:s}?raw=true'
DATA_DIR = 'data'

################################################################
# Bruzual and Charlot (2003) SSP libraries

# directory
BC03 = ( DATA_DIR, 'SSP', 'BC03' )

# BaSeL (Low Resolution) SSP + Chabrier IMF + null extension
bc03_basel_chab_zeros = ( 'basel_chab_zero_extended.dat', BC03  )

# BaSeL (Low Resolution) SSP + Chabrier IMF + null ext. + thinner lambda-grid
bc03_basel_chab_zeros_refined = ( 'basel_chab_zero_refined.dat', BC03 )

# Stelib (High Resolution) SSP + Chabrier IMF + null extension
bc03_stelib_chab_zeros = ( 'stelib_chab_zero_extended.dat', BC03 )

# Stelib (High Resolution) SSP + Chabrier IMF + extrapolation
bc03_stelib_chab_extrap = ( 'stelib_chab_extrap.dat', BC03 )

################################################################
# Bressan SSP libraries

# directory
Br22 = ( DATA_DIR, 'SSP', 'Br22' )

# Custom SISSA SSPs including non-thermal SN emission
br22_NT = ( 'Bressan_NT_formatted.dat', Br22  )

# Custom SISSA SSPs including non-thermal SN emission + Nebular Free-free and lines
br22_NTL = ( 'Bressan_NTL_formatted.dat', Br22  )

# Custom SISSA SSPs including non-thermal SN emission + Nebular Free-free and lines
# + thinner lambda-grid
br22_NT_refined = ( 'Bressan_NT_formatted_refined.dat', Br22  )

# Custom SISSA SSPs including non-thermal SN emission + Nebular Free-free and lines
# + thinner lambda-grid
br22_NTL_refined = ( 'Bressan_NTL_formatted_refined.dat', Br22  )

################################################################
# Photometric Filters

FILT_DIR = ( DATA_DIR, 'filters' )

################################################################
# Cosmology

COSMO_DIR = ( DATA_DIR, 'Cosmologies' )

################################################################
# AGN Templates (Fritz et al., 2006)

AGN_DIR = ( DATA_DIR, 'AGN_templates' )
AGN_FILE = 'ct{0:d}al{1:.1f}be{2:.2f}ta{3:.1f}rm{4:d}ia{5:.3f}'

################################################################
# IGM Fitting values (Inoue et al., 2014)

IGM_DIR = ( DATA_DIR, 'IGM' )
IGM_Inoue14 = ( 'Inoue+2014.txt', IGM_DIR )

################################################################
