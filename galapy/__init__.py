""" This package simulates galactic spectra
"""

from galapy.configuration import rcParams

import galapy.StarFormationHistory
import galapy.CompositeStellarPopulation
import galapy.InterStellarMedium
import galapy.ActiveGalacticNucleus
import galapy.XRayBinaries
import galapy.NebularFreeFree
import galapy.Synchrotron
import galapy.PhotometricSystem
import galapy.Galaxy

__all__ = [
    'Galaxy',
    'StarFormationHistory',
    'CompositeStellarPopulation',
    'InterStellarMedium',
    'ActiveGalacticNucleus',
    'NebularFreeFree',
    'Synchrotron',
    'XRayBinaries',
    'PhotometricSystem'
]

