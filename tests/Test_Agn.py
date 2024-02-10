
#------------------------------------------------------------------------------#
#                  Test of galapy.ActiveGalacticNucleus module.
#------------------------------------------------------------------------------#

"""

Testing the Class "AGN" implementing the AGN contribution by loading templates 
from Fritz et al., 2006"

The AGN is divided into 3 components:

  -  Accretion disk around the central SMBH
  -  Scattered emission by the surrounding dusty torus
  -  Thermal dust emission associated to the dusty torus

The templates are computed accounting for the variation of 6 structural 
parameters:


  - cl : the covering angle of the torus 
  - al, be : parameters of the dust density distribution in spherical 
             coordinates.
  - rm : the ratio between the maximum to minimum radii of the dusty torus;
  - ta : the optical depth at 9.7 micrometers;
  - ia : the angle between the rotation axis and the line-of-sight

"""


# External imports
import numpy as np
import pytest

# Internal imports
import galapy
from galapy import ActiveGalacticNucleus as gpagn
from galapy.internal.constants import Ang_to_keV

def check_2leveldict_approxequal ( A, B ) :
    return np.all([
        va == pytest.approx( vb )
        for va, vb in zip( A.values(), B.values() )
    ])
        

def test_agn_wrong_key () :

    """
    Invalid key should report an error.

    """
    
    key = 'wrong key' 
    with pytest.raises( AttributeError, match = f"Parameter '{key}' is not a valid template parameter."):
        gpagn.find_template_par( key, 1. )

#------------------------------------------------------------------------------#

def test_agn_nearest_params () :
    
    """
    If a key is not present, the nearest one will be chosen. 

    """

    assert gpagn.find_template_par( 'ct', 45. ) == 40

#------------------------------------------------------------------------------#

def test_agn_build_params () :

    """
    Check the module parameters after setting some different from default ones.

    """

    assert check_2leveldict_approxequal(
        gpagn.agn_build_params( 1.e-3 ),
        { 'fAGN': 0.001,
          'template': { 'ct': 40,
                        'al': 0.0,
                        'be': -0.5,
                        'ta': 6.0,
                        'rm': 60,
                        'ia': 0.001 } }
    )

    assert check_2leveldict_approxequal(
        gpagn.agn_build_params( 1.e-3, ct = 60 ),
        { 'fAGN': 0.001,
          'template': { 'ct': 60,
                        'al': 0.0,
                        'be': -0.5,
                        'ta': 6.0,
                        'rm': 60,
                        'ia': 0.001 } }
    )

#------------------------------------------------------------------------------#

@pytest.fixture
def agn () :

    """
    AGN class is initialize with a pytest fixture to avoid repetitions.
    Required parameters: lmin, lmax
  

    """
    
    agn =  gpagn.AGN( 1.e+0, 1.e+10, Xray = True )
    return agn
 
#------------------------------------------------------------------------------#

def test_agn_init ( agn ) :

    """
    Check the correct initialiation of the class.
    
    """

    assert isinstance( agn, galapy.ActiveGalacticNucleus.AGN )

#------------------------------------------------------------------------------#

def test_agn_compute_X_template ():

    """
    lmin must be smaller than 2 keV.
    
    """

    with pytest.raises( RuntimeError, match =
                        "Cannot build the X-ray spectrum for "
                                "a wavelength grid starting at lambda > "
                                "6 Angstrom ~ 2 keV! "
                                "Set a smaller `lmin` value."):

        agn=gpagn.AGN( lmin = 1.e+2, lmax = 1.e+10 )

#------------------------------------------------------------------------------#

def test_agn_set_parameters (agn) :

    """
    Set the class parameters.
    
    """

    agn.set_parameters( fAGN = 0.1, ia = 65.1 )
    assert check_2leveldict_approxequal(
        agn.params,
        { 'fAGN': 0.1,
          'template': { 'ct': 40,
                        'al': 0.0,
                        'be': -0.5,
                        'ta': 6.0,
                        'rm': 60,
                        'ia': 60.1 } }
    )

#------------------------------------------------------------------------------#

def test_agn_emission (agn):
   
    """
    Check the AGN emission at different wavelengths.
    
    """
 
    ll= np.logspace(0, 8, 5 )
    assert np.all ( agn.emission( ll, 1.e+12 ) ==
                    pytest.approx( [ 7.56985544e+09,
                                     1.18091098e-14,
                                     1.84215088e-05,
                                     5.14417599e+01,
                                     0. ] ) )

    
