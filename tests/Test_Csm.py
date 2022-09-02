
#------------------------------------------------------------------------------#
#                      Test of galapy.Cosmology module.
#------------------------------------------------------------------------------#

"""

Testing the class CSM whose just a wrapper for interpolation on grid of 
Luminosity Distances and Ages.

    Parameters:
    ----------

    Cosmo: string or formatted dictionary.
    Cosmological parameters used as a function of redshift.
    Formatted as a table: Redshift, Luminosity Distance, Age of the universe.


"""


# External imports
import numpy as np
import pytest

# Internal imports
import galapy
from galapy import Cosmology as gpcsm


@pytest.fixture
def csm () :

    """
    CSM class is initialize with a pytest fixture to avoid repetitions.
    Required parameters: cosmology 
  

    """
    
    csm =  gpcsm.CSM( 'Planck18' )
    return csm

#------------------------------------------------------------------------------#

def test_csm_init ( csm ):

    """
    Check the correct initialiation of the class.
    
    """

    assert isinstance( csm, galapy.Cosmology.CSM )

#------------------------------------------------------------------------------#


# def test_csm_error ( ):

#     cosmo={}
#     with pytest.raises( ValueError, match =  'Argument `cosmo` should be either a string '
#                               'or a formatted dictionary.'):
#       csm=gpcsm.CSM(cosmo)
    
