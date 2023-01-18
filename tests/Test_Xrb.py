
#------------------------------------------------------------------------------#
#                      Test of galapy.XRayBinaries module.
#------------------------------------------------------------------------------#

"""

Testing the Class XRB which implements the contribution to the X-band emission 
from stellar binary systems. 
  
    Parameters
    ----------
    - HM (HIGH MASS STARS):
            
    Psi: Star Formation History;
    Zstar: Star Metallicity.

    ----------
    - LM (LOW MASS STARS):
            
    Age: Stellar Age;
    Mstar: Stellar Mass.


"""



# External imports
import numpy as np
import pytest

# Internal imports
import galapy
from galapy import XRayBinaries as gpxrb


def test_xrb_build_params_wrong_startype () :

    """
    Startype accepted are only "hm" or "lm".

    """
    
    startype='wrong_startype'
    av_types='"hm" "lm" '
    with pytest.raises( ValueError, match =  f'Star type "{startype}" chosen not available. '
                          f'Available phases are: {av_types}'):

        gpxrb.xrb_build_params(startype)

#------------------------------------------------------------------------------#
        
def test_xrb_build_params_hm () :

    """
    Build module parameters for high mass stars using default values or 
    different ones.

    """
    assert gpxrb.xrb_build_params('hm') == {'psi': 1.,
                                            'Zstar': 0.02 }

    assert gpxrb.xrb_build_params('hm', psi = 1.5, Zstar = 0.05 ) == {'psi': 1.5,  'Zstar': 0.05 }

#------------------------------------------------------------------------------#

def test_xrb_build_params_lm () :

    """
    Build module parameters for low mass stars using default values or 
    different ones.

    """
    assert gpxrb.xrb_build_params('lm') == { 'Mstar' : 1.e+10
                                             ,'age'   : 1.e+8 }

    assert gpxrb.xrb_build_params('lm', Mstar = 1.e+8, age = 1.e+7 ) == { 'Mstar' : 1.e+8,'age'   : 1.e+7 }

#------------------------------------------------------------------------------#

@pytest.fixture
def xrb () :
    
    """
    XRB class is initialize with a pytest fixture to avoid repetitions.
    Required parameters: wavelength grid
  
    """
    
    xrb =  gpxrb.XRB( 1., 1.e+10 )
    return xrb

#------------------------------------------------------------------------------#

def test_xrb_init ( xrb ) :
   
    """
    Check the correct initialiation of the class.
    
    """
 
    assert isinstance( xrb, galapy.XRayBinaries.XRB )

#------------------------------------------------------------------------------#
 
def test_xrb_hm_set_parameters ( xrb ) :

    """
    Set the class parameters for high mass stars.
    
    """

    xrb.hm.set_parameters( psi = 1.5 )
    assert xrb.hm.params == {'psi': 1.5, 'Zstar': 0.02 }

#------------------------------------------------------------------------------#
 
def test_xrb_lm_set_parameters ( xrb ):

    """
    Set the class parameters for low mass stars.
    
    """

    xrb.lm.set_parameters( Mstar = 1.e+8 )
    assert xrb.lm.params ==  { 'Mstar' : 1.e+8, 'age' : 1.e+8 }

#------------------------------------------------------------------------------#

def test_xrb_emission ( xrb ):
   
    """
    Check the X-ray emission from binary systems at different wavelengths.
    
    """
    
    ll = np.logspace( 0, 3, 5 )
    assert np.all( xrb.emission( ll ) == pytest.approx( [ 11087592.11500118,
                                                          1108286.39951714,
                                                          102964.04627411,
                                                          0.,
                                                          0. ] ) )

#------------------------------------------------------------------------------#
    
def test_xrb_emission_sum ( xrb ) : 
   
    """
    Check that the total X-ray emission from binary systems is the sum of 
    the emission from low and high mass startype
    
    """
    
    hm_em = xrb.hm.emission( 1.e-1 )
    lm_em = xrb.lm.emission( 1.e-1 )
    em = hm_em + lm_em
    assert xrb.emission( 1.e-1 ) == pytest.approx( em )
