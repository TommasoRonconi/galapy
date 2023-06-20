""" Definitions of physical and mathematical constants
"""

# Solar luminosity [ erg / s ]
Lsun = 3.828e+33
sunL = 1. / Lsun

# Speed of light in vacuum
clight = {
    'cm/s' : 2.99792458e+10,
    'm/s'  : 2.99792458e+8,
    'A/s'  : 2.99792458e+18,
}

# Mpc to cm
Mpc_to_cm = 3.086e+24

# Planck constant
hP = {
    'eV/Hz' : 4.1357e-15,
    'erg*s' : 6.6262e-27,
}

# Angstrom to keV
def Ang_to_keV ( wavelength ) :
    return 1.e-3 * hP['eV/Hz'] * clight['A/s'] / wavelength

