
#------------------------------------------------------------------------------#
#                    Test of galapy.NebularFreeFree module.
#------------------------------------------------------------------------------#

"""

Testing the Galapy module providing the Nebular Emission component. (see Murphy
et al. 2011; Bressan et al. 2002; Mancuso et al. 2017)

    Parameters
    ----------

    Zgas: Gas metallicity

    Zi: Gas atomic number ( = 1 for a pure H plasma )

"""


# External imports
import numpy as np
import pytest

# Internal imports
import galapy
from galapy import NebularFreeFree as gpnff
from galapy import CompositeStellarPopulation as gpcsp




@pytest.fixture
def nff () :
    
    """
    NFF class is initialize with a pytest fixture to avoid repetitions.
    Required parameters: wavelength grid
  
    """
    
    ll = np.logspace( 0., 10., 5 )
    nff =  gpnff.NFF( ll )
    return nff
 
#------------------------------------------------------------------------------#

def test_nff_init ( nff ) :
   
    """
    Check the correct initialiation of the class.
    
    """
 
    assert isinstance( nff, galapy.NebularFreeFree.NFF ) 
    
#------------------------------------------------------------------------------#

def test_nff_build_params () :

    """
    Check the module parameters after initializing the class with parameter 
    values different from default ones.

    """
    ll = ( 1.e+0, 1.e+10 )
    nff = gpnff.NFF( ll, Zgas = 0.01 )
    assert nff.params == { 'Zgas' : 0.01, 'Zi' : 1.0 }

#------------------------------------------------------------------------------#

def test_nff_set_params ( nff ) :
    
    """
    Set the class parameters.
    
    """

    nff.set_parameters( Zgas = 0.01 )
    assert nff.params == { 'Zgas' : 0.01, 'Zi' : 1.0 }

#------------------------------------------------------------------------------#

def test_nff_Te ( nff ) :

    """
    Check the electron temperature for a certain gas metallicity

    """

    Tem = nff.Te( Zgas = 0.1 )
    assert Tem == 3502.185479394832

#------------------------------------------------------------------------------#

def test_nff_gff ( nff ) :

    """
    Check the value of the gaunt factor for a certain wavelengt grid and at a 
    certain electron temperature.

    """

    gff = nff.gff()
    assert np.all( gff == pytest.approx( [ 1.00039697,
                                           1.00944316,
                                           1.20434471,
                                           2.85905773,
                                           5.87106787 ] ) )

#------------------------------------------------------------------------------#

def test_nff_emission ( nff ) :
   
    """
    Check the NFF emission at different wavelengths for a set intrinsic 
    photoionization rate Q_H.
    
    """
    
    assert np.all( nff.emission( 1.e+20 ) == pytest.approx( [ 0.,
                                                              1.76571524e-19,
                                                              5.00431080e+01,
                                                              1.42909655e-03,
                                                              2.93636138e-08 ]
    ) )

