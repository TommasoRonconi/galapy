
#------------------------------------------------------------------------------#
#                  Test of galapy.PhotometricSystem module.
#------------------------------------------------------------------------------#

"""
Testing the class PMS that builds the Photometric System. 
    
    Parameters
    ----------
    *args : strings
       As many strings as necessary. Each string has to name
       one of the available filters already present in the database.
    **kwargs : 'named' dictionaries
       Keyword arguments can be used for providing user-defined 
       transmissions. 
       Custom transmissions are passed as properly formatted dictionaries:
       keyword = { 'wavelengths' : array-like,
                   'photons' : array-like }
       The two arrays provided, as the keys suggest, must define the
       wavelength grid and the corresponding transmission in photon units
       (and thus they must have the same size).
       Note that the chosen keyword will be used as unique identifier of
       the custom transmission.

"""


# External imports
import numpy as np
import pytest

# Internal imports
import galapy
from galapy import PhotometricSystem as gppms


@pytest.fixture
def pms ():
         
    """
    PMS class is initialize with a pytest fixture to avoid repetitions.
  
    """
    
    pms =  gppms.PMS( )
    return pms

#------------------------------------------------------------------------------#

def test_pms_init ( pms ) :
   
    """
    Check the correct initialiation of the class.
    
    """

    assert isinstance( pms, galapy.PhotometricSystem.PMS )

#------------------------------------------------------------------------------#

# def test_pms_error () :

#      """
#      A ValueError should be reported if the filter name is not present in the 
#      database.

#      """
     
#      arg = 'wrong_name'
#      with pytest.raises( ValueError, match = f'Argument "{arg}" provided does not name a known filter.' ):
#         pms = gppms.PMS( arg )

#------------------------------------------------------------------------------#


def test_pms_get_fluxes ():
     
     """
     Check the value of the average-flux in a certain band using a filter 
     included in the database.

     """

     pms = gppms.PMS( 'GOODS.b' )
     ll = np.logspace( 3.5, 3.7, 100 )
     fl = np.logspace( 0., 5. , 100 )
     assert pms.get_fluxes( ll, fl ) == pytest.approx( [ 6314.613587755724 ] )
