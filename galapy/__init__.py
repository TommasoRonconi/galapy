"""GalaPy models the Spectral Energy Distribution (SED) of galaxies from the X-ray to the Radio band.
 
It is an API as well as the emission itself, models the evolution of 
galaxy components and dust attenuation/reradiation. 
GalaPy provides functions, classes and terminal commands for Bayesian inference of 
galaxy properties from panchromatic photometric data, as well as for the analysis 
and simulation of galaxies spectral properties.
"""

from galapy.configuration import rcParams

from galapy import StarFormationHistory
from galapy import CompositeStellarPopulation
from galapy import InterStellarMedium
from galapy import NebularFreeFree
from galapy import Synchrotron
from galapy import XRayBinaries
from galapy import ActiveGalacticNucleus
from galapy import Cosmology
from galapy import Galaxy
from galapy import InterGalacticMedium
from galapy import PhotometricSystem
from galapy import Noise
from galapy import Handlers
from galapy import sampling, internal, io, analysis

__all__ = [
    'StarFormationHistory', 'CompositeStellarPopulation', 'InterStellarMedium',
    'NebularFreeFree', 'Synchrotron', 'XRayBinaries', 'ActiveGalacticNucleus',
    'Cosmology', 'Galaxy', 'InterGalacticMedium', 'PhotometricSystem',
    'Noise', 'Handlers',
    'sampling', 'internal', 'io', 'analysis',
]
__version__ = 'v0.5.1'
