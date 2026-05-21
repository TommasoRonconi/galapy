
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
from galapy.AGN_core.Fritz2006 import Fritz2006, find_template_par, agn_build_params
from galapy.AGN_core.Panchromatic import Panchromatic
from galapy.internal.constants import Ang_to_keV

# Shared lgrid used throughout
_LGRID = np.logspace(0, 10, 500)

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
        find_template_par( key, 1. )

#------------------------------------------------------------------------------#

def test_agn_nearest_params () :

    """
    If a key is not present, the nearest one will be chosen.

    """

    assert find_template_par( 'ct', 45. ) == 40

#------------------------------------------------------------------------------#

def test_agn_build_params () :

    """
    Check the module parameters after setting some different from default ones.

    """

    assert check_2leveldict_approxequal(
        agn_build_params( 1.e-3 ),
        { 'fAGN': 0.001,
          'template': { 'ct': 40,
                        'al': 0.0,
                        'be': -0.5,
                        'ta': 6.0,
                        'rm': 60,
                        'ia': 0.001 } }
    )

    assert check_2leveldict_approxequal(
        agn_build_params( 1.e-3, ct = 60 ),
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
    """Fritz2006 AGN instance shared across tests."""
    return gpagn.AGN( _LGRID, model = 'Fritz2006', do_Xray = True )
 
#------------------------------------------------------------------------------#

def test_agn_init ( agn ) :

    """
    Check the correct initialiation of the class.
    
    """

    assert isinstance( agn, galapy.ActiveGalacticNucleus.AGN )

#------------------------------------------------------------------------------#

def test_agn_compute_X_template ():

    """
    lgrid must start below 6 Angstrom (~2 keV) when do_Xray=True.

    """

    with pytest.raises( RuntimeError, match =
                        "Cannot build the X-ray spectrum for "
                        "a wavelength grid starting at lambda > "
                        "6 Angstrom ~ 2 keV! "
                        "Set a smaller `lmin` value." ):
        gpagn.AGN( np.logspace( 2, 10, 500 ), model = 'Fritz2006', do_Xray = True )

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
    Check Fritz2006 AGN emission values at reference wavelengths.
    Last point (1e8 Ang) is outside the template domain and must be zero.

    """

    ll = np.logspace( 0, 8, 5 )     # [1, 1e2, 1e4, 1e6, 1e8] Ang
    em = agn.emission( ll, 1.e+12 )
    assert em.shape == ll.shape
    assert not np.any( np.isnan( em ) )
    assert np.all( em >= 0. )
    expected = np.array( [ 7.56985544e+09, 1.18091098e-14,
                            1.84215088e-05, 5.14417599e+01,
                            0.00000000e+00 ] )
    np.testing.assert_allclose( em, expected, rtol = 1e-5 )


def test_agn_invalid_model () :

    """
    Requesting an unknown AGN model should raise RuntimeError.

    """

    with pytest.raises( RuntimeError, match = 'not available' ) :
        gpagn.AGN( _LGRID, model = 'NonExistent' )


# ===========================================================================
# Panchromatic model
# ===========================================================================

@pytest.fixture
def agn_pan () :
    """Panchromatic AGN instance shared across tests."""
    return gpagn.AGN( _LGRID, model = 'Panchromatic' )


def test_agn_panchromatic_init ( agn_pan ) :

    """Panchromatic wrapper builds correctly."""

    assert isinstance( agn_pan, galapy.ActiveGalacticNucleus.AGN )
    assert isinstance( agn_pan.core, Panchromatic )
    assert agn_pan.params['model'] == 'Panchromatic'


def test_agn_panchromatic_params ( agn_pan ) :

    """Panchromatic params dict exposes all expected physical parameters."""

    expected = { 'Lbol', 'theta_view', 'delta', 'TH', 'Delta', 'EBV', 'RL' }
    assert expected <= set( agn_pan.params.keys() )


def test_agn_panchromatic_set_parameters ( agn_pan ) :

    """set_parameters updates the internal state without error."""

    agn_pan.set_parameters( Lbol = 1.e+45, theta_view = 0.5 )
    assert agn_pan.params['Lbol'] == pytest.approx( 1.e+45 )
    assert agn_pan.params['theta_view'] == pytest.approx( 0.5 )


def test_agn_panchromatic_emission ( agn_pan ) :

    """Panchromatic emission matches reference values over the full wavelength range."""

    ll = np.logspace( 0, 8, 20 )
    em = agn_pan.emission( ll )
    assert em.shape == ll.shape
    assert not np.any( np.isnan( em ) )
    assert np.all( em >= 0. )
    expected = np.array( [ 1.40697062e+09, 4.90034180e+08, 1.56345563e+08,
                            4.74456243e+07, 1.14968625e+07, 5.94710472e+06,
                            1.13201935e+07, 1.85147815e+07, 7.55356380e+06,
                            1.24909433e+06, 2.21994768e+05, 5.66239762e+04,
                            6.63297207e+04, 3.14606799e+04, 4.11599468e+03,
                            2.60765330e+03, 8.06368988e+02, 1.93090266e+02,
                            4.37826172e+01, 9.47538245e+00 ] )
    np.testing.assert_allclose( em, expected, rtol = 1e-5 )

    
