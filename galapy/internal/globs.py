import os

################################################################
#                      DATA DIRECTORY:                         #
# contains binary and ASCII files to be used along the library #
################################################################

DATABASE = 'https://api.github.com/repos/TommasoRonconi/galapy_database/releases'
DATA_VERSION = '0.1.0'
DATA_URL = 'https://raw.githubusercontent.com/TommasoRonconi/galapy_database/main/{:s}'
DATA_DIR = 'galapy_database'

################################################################
# Simple Stellar Population libraries

SSP_DIR = ( DATA_DIR, 'SSP' )

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
