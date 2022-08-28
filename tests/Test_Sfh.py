
#------------------------------------------------------------------------------#
#                  Test of galapy.StarFormationHistory module.
#------------------------------------------------------------------------------#

"""

Testing the Class "SFH" wrapping the C-core implementation of the 
Star Formation History type .    

    Parameters
    ----------
    tau_quench : float
      Eventual abrupt quenching time for star formation ( yr ) . 
      It is expressed in terms of the time passed from the formation 
      of the galaxy .
      Default is 2.e+20 . 

    model : string
      One among ( 'insitu', 'constant', 'delayedexp', 'lognormal', 'burst' ) .
      Default is 'insitu'.

For the empirical SFH models (i.e. constant, delayed-exponential and log-normal)
two further parameters can be tuned :

    ----------
    Mdust : float
      Total dust content of the galaxy ( M_sun ) .
      Default is 1.e+8 .

    Zgs : float
      Average gas metallicity of the components .
      The metallicity of the stellar and gas components is the same . 
      Default is 1% .

"""

# External imports
import numpy as np
import pytest

# Internal imports
import galapy
from galapy import StarFormationHistory as gpsfh


#Chosen value for tau_quench: 800000000.0 yr

tau_quench=8.e+8




def test_sfh_init_models () :
    
    """
    Only models indicate above are valid choices for initialize the class SFH .
 
    """
    
    with pytest.raises( TypeError, match = "SFH model not valid. Valid models are: 'insitu', 'constant', 'delayedexp', 'lognormal', 'burst'" ):
        sfh = gpsfh.SFH( model = '' )
        
#------------------------------------------------------------------------------#


################################################################################
############################## SFH MODEL = INSITU ##############################
################################################################################


"""
    In-situ model additional parameters
    ----------
    psi_max : float
      Normalization of the SFR ( M_sun / yr ) .
      Default is 100.0 .

    tau_star : float
      Characteristic time ( yr ) .
      Default is 3.e+8 .

"""

@pytest.fixture
def sfh_insitu ():
    
     sfh_insitu = gpsfh.SFH( tau_quench = tau_quench, model = 'insitu' )
     return sfh_insitu
 
def test_sfh_insitu_init ( sfh_insitu ) :
    
    """
    Test initialization of the class SFH .
 
    """
    
    assert isinstance( sfh_insitu, galapy.StarFormationHistory.SFH )

#------------------------------------------------------------------------------#

#@pytest.mark.parametrize("models",['insitu','constant','delayedexp','lognormal'])
def test_sfh_insitu_parameters ( ):

    """
    Test default parameters for 'insitu' model .

    """
    sfh_insitu = gpsfh.SFH( model = 'insitu' )
    assert sfh_insitu.params == { 'tau_quench': 2.e+20, 'model': 'insitu',
                                  'psi_max': 100.0, 'tau_star': 3.e+8 }

#------------------------------------------------------------------------------#

def test_sfh_insitu_call ( sfh_insitu ) :
    
    """ 
    Test calling the class at specific time ( array-like or float type ) .

    """

    assert sfh_insitu( 1.e+8 ) == pytest.approx( 24.726516804071363 )
    tau = np.logspace( 6, 10, 4 )
    assert np.all( sfh_insitu( tau ) == pytest.approx( [ 0.32252989,
                                                         6.57127224,
                                                         46.399231, 0. ] ) )
#------------------------------------------------------------------------------#

def test_sfh_insitu_psimax ( sfh_insitu ) :
    
    """ 
    Normalization factor equal to zero will return a null sfh . 
    
    """
    
    sfh_insitu.set_parameters( psi_max = 0. )
    tau = np.logspace( 6, 10, 50 )
    psi_insitu = sfh_insitu( tau )
    assert all ( [ a == 0. for a in psi_insitu ] )

#------------------------------------------------------------------------------#

def test_sfh_insitu_quench ( sfh_insitu ) :
    
    """
    After tau_quench, sfh must be equal to 0 .

    """
    tau = np.logspace( 9, 11, 100 )
    psi_insitu = sfh_insitu( tau )
    assert all ( [ a == 0. for a in psi_insitu ] )
    
#------------------------------------------------------------------------------#

def test_sfh_insitu_Mstar ( sfh_insitu ) :
    
    """
    Check the value of the stellar mass at quench time and at present time .

    """
    assert sfh_insitu.Mstar( 8.e+8 ) == pytest.approx( 20694798192.00076 ) 
    assert sfh_insitu.Mstar( 1.e+9 ) == pytest.approx( 19645102756.25652 ) 

#------------------------------------------------------------------------------#

def test_sfh_insitu_Mdust ( sfh_insitu ) :
    
    """
    Check the value of the total dust content at quench and at present time .

    """
    assert sfh_insitu.Mdust( 8.e+8 ) == pytest.approx( 140171330.69560972 ) 
    assert sfh_insitu.Mdust( 1.e+9 ) == 0.
    
#------------------------------------------------------------------------------#

def test_sfh_insitu_Mgas ( sfh_insitu ) :
    
    """
    Check the value of total gas content at quench and at present time .
    
    """
    assert sfh_insitu.Mgas( 8.e+8 ) == pytest.approx( 5702802608.720548 ) 
    assert sfh_insitu.Mgas( 1.e+9 ) == 0.
    
#------------------------------------------------------------------------------#

def test_sfh_insitu_Zstar ( sfh_insitu ) :
    
    """
    Check the value of the stellar metallicity at present time ( ~0.02 ) .

    """
    assert sfh_insitu.Zstar( 1.e+9 ) == pytest.approx( 0.023774376958190433 )
    
#------------------------------------------------------------------------------#

def test_sfh_insitu_Zgas ( sfh_insitu ) :
    
    """
    Check the value of the gas metallicity at present time ( ~0.034 ) .
    """
    assert sfh_insitu.Zgas( 1.e+9 ) == pytest.approx( 0.03575910183762231 ) 

#------------------------------------------------------------------------------#

#Chosen values for Mdust and Zgs for empirical SFH models:

md = 1.e+7
zz = 0.02

################################################################################
############################## SFH MODEL = CONSTANT ############################
################################################################################


"""
    Constant model additional parameters
    ----------
    psi : float
      SFR constant value ( M_sun / yr ) .
      Default is 1.0 .

"""
@pytest.fixture
def sfh_const ():
    
     sfh_const = gpsfh.SFH( tau_quench = tau_quench, model = 'constant',
                            Mdust = md, Zgs = zz )
     return sfh_const
 

def test_sfh_const_init ( sfh_const ) :
    
    """
    Test initialization of the class SFH.
 
    """
    assert isinstance( sfh_const, galapy.StarFormationHistory.SFH )

#------------------------------------------------------------------------------#

def test_sfh_const_parameters( ):

    """
    Test default parameters for 'constant' model .

    """
    
    sfh_const = gpsfh.SFH( model = 'constant' )
    assert sfh_const.params == { 'tau_quench': 2e+20, 'model': 'constant',
                                 'psi': 1.0, 'Mdust': 100000000.0, 'Zgs': 0.1 }

#------------------------------------------------------------------------------#

def test_sfh_const_call ( sfh_const ) :
    
    """ 
    Test calling the class at specific time ( array-like or float type ) .

    """
    assert sfh_const( 1.e+8 ) == 1.
    tau = np.logspace( 6, 10, 4 )
    assert np.all( sfh_const( tau ) == pytest.approx( [ 1., 1., 1., 0. ] ) )

#------------------------------------------------------------------------------#

def test_sfh_const_quench ( sfh_const ) :
    
    """
    After tau_quench, sfh must be equal to 0 .

    """
    
    tau = np.logspace( 9, 11, 100 )
    psi_const = sfh_const( tau )
    assert all ( [ a == 0. for a in psi_const ] )
    
#------------------------------------------------------------------------------#

def test_sfh_const_Mstar ( sfh_const ) :
    
    """
    Check the value of the stellar mass at quench and at present time .

    """
    assert sfh_const.Mstar( 8.e+8 ) == pytest.approx( 535695178.04297197 ) 
    assert sfh_const.Mstar( 1.e+9 ) == pytest.approx( 509519015.14239407 ) 

#------------------------------------------------------------------------------#

def test_sfh_const_Mdust ( sfh_const ) :
    
    """
    Check the value of total dust content at quench and at present time
    ( must be equal to md ) .

    """
    assert sfh_const.Mdust( 8.e+8 ) == md 
    assert sfh_const.Mdust( 1.e+9 ) == md
    
#------------------------------------------------------------------------------#

def test_sfh_const_Mgas ( sfh_const ) :
    
    """
    Check the value of total gas content at present time .
    
    """
    assert sfh_const.Mgas( 1.e+9 ) == pytest.approx( 1405518018.9671905 )
    
#------------------------------------------------------------------------------#

def test_sfh_const_Zstar ( sfh_const ) :
    
    """
    Check the value of the stellar metallicity at present time 
    ( must be equal to zz ) .

    """
    assert sfh_const.Zstar( 1.e+9 ) == zz
    
#------------------------------------------------------------------------------#

def test_sfh_const_Zgas ( sfh_const ) :
    
    """
    Check the value of the Gas Metallicity at present time 
    ( must be equal to zz ) .
    """
    assert sfh_const.Zgas( 1.e+9 ) == zz 

#------------------------------------------------------------------------------#

################################################################################
######################### SFH MODEL = DELAYED EXPONENTIAL ######################
################################################################################

"""
    Delayed Exponential model additional parameters
    ----------
    psi_norm : float
      Normalization of the SFR ( M_sun / yr ) .
      Default is 1. .
    k : float
      Shape parameter ( adimensional ) .
      Default is 0.2 . 
    tau_star : float
      Characteristic time ( yr ) .
      Default is 3.e+8 .

"""
@pytest.fixture
def sfh_dexp ():
    
     sfh_dexp = gpsfh.SFH( tau_quench = tau_quench, model = 'delayedexp',
                            Mdust = md, Zgs = zz )
     return sfh_dexp
 

def test_sfh_dexp_init ( sfh_dexp ) :
    
    """
    Test initialization of the class SFH .
 
    """
    assert isinstance( sfh_dexp, galapy.StarFormationHistory.SFH )

#------------------------------------------------------------------------------#

def test_sfh_dexp_parameters () :

    """
    Test default parameters for 'delayedexp' model .

    """
    
    sfh_dexp = gpsfh.SFH( model = 'delayedexp' )
    assert sfh_dexp.params == { 'tau_quench': 2e+20, 'model': 'delayedexp',
                                'psi_norm': 1.0, 'k_shape': 0.2,
                                'tau_star': 100000000.0, 'Mdust': 100000000.0,
                                'Zgs': 0.1 }

#------------------------------------------------------------------------------#

def test_sfh_dexp_call ( sfh_dexp ) :
    
    """ 
    Test calling the class at specific time ( array-like or float type ) .

    """
    assert sfh_dexp( 1.e+8 ) == 14.645544342956468
    tau = np.logspace( 6, 10, 4 )
    assert np.all( sfh_dexp( tau ) == pytest.approx( [ 15.69123242, 23.61025931,
                                                       0.52181543, 0. ] ) )

#------------------------------------------------------------------------------#

def test_sfh_dexp_quench ( sfh_dexp ) :
    
    """
    After tau_quench, sfh must be equal to 0 .

    """
    tau = np.logspace( 9, 11, 100 )
    psi_dexp = sfh_dexp( tau )
    assert all ( [ a == 0. for a in psi_dexp ] )
    
#------------------------------------------------------------------------------#

def test_sfh_dexp_Mstar ( sfh_dexp ) :
    
    """
    Check the value of the stellar mass at quench time and at present time .

    """
    assert sfh_dexp.Mstar( 8.e+8 ) == pytest.approx( 2302675603.296462 ) 
    assert sfh_dexp.Mstar( 1.e+9 ) == pytest.approx( 2255427119.0862346 ) 

#------------------------------------------------------------------------------#

def test_sfh_dexp_Mdust ( sfh_dexp ) :
    
    """
    Check the value of total dust content at quench and at present time .
    ( must be equal to md ) .

    """
    assert sfh_dexp.Mdust( 8.e+8 ) == md 
    assert sfh_dexp.Mdust( 1.e+9 ) == md
    
#------------------------------------------------------------------------------#

def test_sfh_dexp_Mgas ( sfh_dexp ) :
    
    """
    Check the value of total gas content at present time .
    
    """
    assert sfh_dexp.Mgas( 1.e+9 ) == pytest.approx( 1405518018.9671905 )
    
#------------------------------------------------------------------------------#

def test_sfh_dexp_Zstar ( sfh_dexp ) :
    
    """
    Check the value of the stellar metallicity at present time 
    ( must be equal to zz ) .

    """
    assert sfh_dexp.Zstar( 1.e+9 ) == zz
    
#------------------------------------------------------------------------------#

def test_sfh_dexp_Zgas ( sfh_dexp ) :
    
    """
    Check the value of the gas metallicity at present time 
    ( must be equal to zz ) .
    """
    assert sfh_dexp.Zgas( 1.e+9 ) == zz 

#------------------------------------------------------------------------------#


################################################################################
############################# SFH MODEL = LOG-NORMAL ###########################
################################################################################

"""
    Log-normal model additional parameters
    ----------
    psi_norm : float
      Normalization of the SFR ( M_sun / yr ) .
      Default is 100. .
    sigma_star : float
      Time span of the evolution ( adimensional ) .
      Default is 2. . 
    tau_star : float
      Characteristic time ( yr ) .
      Default is 3.e+8 .

"""
@pytest.fixture
def sfh_lnorm ():
    
     sfh_lnorm = gpsfh.SFH( tau_quench = tau_quench, model = 'lognormal',
                            Mdust = md, Zgs = zz )
     return sfh_lnorm
 

def test_sfh_lnorm_init ( sfh_lnorm ) :
    
    """
    Test initialization of the class SFH .
 
    """
    assert isinstance( sfh_lnorm, galapy.StarFormationHistory.SFH )

#------------------------------------------------------------------------------#

def test_sfh_lnorm_parameters () :

    """
    Test default parameters for 'lognormal' model .

    """
    
    sfh_lnorm = gpsfh.SFH( model = 'lognormal' )
    assert sfh_lnorm.params == { 'tau_quench': 2e+20, 'model': 'lognormal',
                                 'psi_norm': 100.0, 'sigma_star': 2.0,
                                 'tau_star': 300000000.0, 'Mdust': 100000000.0,
                                 'Zgs': 0.1 }


#------------------------------------------------------------------------------#

def test_sfh_lnorm_call ( sfh_lnorm ) :
    
    """ 
    Test calling the class at specific time ( array-like or float type ) .

    """
    assert sfh_lnorm( 1.e+8 ) == 17.153733592792385
    tau = np.logspace( 6, 10, 4 )
    assert np.all( sfh_lnorm( tau ) == pytest.approx( [ 0.34179049, 8.38175995,
                                                        19.47777367, 0. ] ) )

#------------------------------------------------------------------------------#

def test_sfh_lnorm_quench ( sfh_lnorm ) :
    
    """
    After tau_quench, sfh must be equal to 0 .

    """
    tau = np.logspace( 9, 11, 100 )
    psi_lnorm = sfh_lnorm( tau )
    assert all ( [ a == 0. for a in psi_lnorm ] )
    
#------------------------------------------------------------------------------#

def test_sfh_lnorm_Mstar ( sfh_lnorm ) :
    
    """
    Check the value of the stellar mass at quench and at present time .

    """
    assert sfh_lnorm.Mstar( 8.e+8 ) == pytest.approx( 9751556099.681376 ) 
    assert sfh_lnorm.Mstar( 1.e+9 ) == pytest.approx( 9273433696.49997 ) 

#------------------------------------------------------------------------------#

def test_sfh_lnorm_Mdust ( sfh_lnorm ) :
    
    """
    Check the value of total dust content at quench time and at present time 
    ( must be equal to md ) .

    """
    assert sfh_lnorm.Mdust( 8.e+8 ) == md 
    assert sfh_lnorm.Mdust( 1.e+9 ) == md
    
#------------------------------------------------------------------------------#

def test_sfh_lnorm_Mgas ( sfh_lnorm ) :
    
    """
    Check the value of total gas content at present time .
    
    """
    assert sfh_lnorm.Mgas( 1.e+9 ) == pytest.approx( 1405518018.9671905 )
    
#------------------------------------------------------------------------------#

def test_sfh_lnorm_Zstar ( sfh_lnorm ) :
    
    """
    Check the value of the stellar metallicity at present time 
    ( must be equal to zz ) .

    """
    assert sfh_lnorm.Zstar( 1.e+9 ) == zz
    
#------------------------------------------------------------------------------#

def test_sfh_lnorm_Zgas ( sfh_lnorm ) :
    
    """
    Check the value of the Gas Metallicity at present time 
    ( must be equal to zz ) .
    """
    assert sfh_lnorm.Zgas( 1.e+9 ) == zz 

#------------------------------------------------------------------------------#
