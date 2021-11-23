
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


def test_sfh_insitu_init () :
    
    """
    Test initialization of the class SFH .
 
    """
    
    sfh_insitu = gpsfh.SFH( model = 'insitu' )
    assert isinstance( sfh_insitu, galapy.StarFormationHistory.SFH )

#------------------------------------------------------------------------------#

def test_sfh_insitu_parameters () :

    """
    Test default parameters for 'insitu' model .

    """
    
    sfh_insitu = gpsfh.SFH( model = 'insitu' )
    assert sfh_insitu.params == { 'tau_quench': 2.e+20, 'model': 'insitu',
                                  'psi_max': 100.0, 'tau_star': 3.e+8 }

#------------------------------------------------------------------------------#

def test_sfh_insitu_call () :
    
    """ 
    Test calling the class at specific time ( array-like or float type ) .

    """
    
    sfh_insitu = gpsfh.SFH( tau_quench = tau_quench, model = 'insitu' )
    assert sfh_insitu( 1.e+8 ) == pytest.approx( 28.37004662 )
    tau = np.logspace( 6, 10, 4 )
    assert np.all( sfh_insitu( tau ) == pytest.approx( [ 0.4510077, 8.19005838,
                                                         48.25602023, 0. ] ) )
#------------------------------------------------------------------------------#

def test_sfh_insitu_psimax () :
    
    """ 
    Normalization factor equal to zero will return a null sfh . 
    
    """
    
    sfh_insitu = gpsfh.SFH( tau_quench = tau_quench, model = 'insitu',
                           psi_max = 0. )
    tau = np.logspace( 6, 10, 50 )
    psi_insitu = sfh_insitu( tau )
    assert all ( [ a == 0. for a in psi_insitu ] )

#------------------------------------------------------------------------------#

def test_sfh_insitu_quench () :
    
    """
    After tau_quench, sfh must be equal to 0 .

    """
    
    sfh_insitu = gpsfh.SFH( tau_quench = tau_quench, model = 'insitu' )
    tau = np.logspace( 9, 11, 100 )
    psi_insitu = sfh_insitu( tau )
    assert all ( [ a == 0. for a in psi_insitu ] )
    
#------------------------------------------------------------------------------#

def test_sfh_insitu_Mstar () :
    
    """
    Check the value of the stellar mass at quench time and at present time .

    """
    
    sfh_insitu = gpsfh.SFH( tau_quench = tau_quench, model = 'insitu' )
    assert sfh_insitu.Mstar( 8.e+8 ) == pytest.approx( 21817709705.887363 ) 
    assert sfh_insitu.Mstar( 1.e+9 ) == pytest.approx( 20730890486.269768 ) 

#------------------------------------------------------------------------------#

def test_sfh_insitu_Mdust () :
    
    """
    Check the value of the total dust content at quench and at present time .

    """
    sfh_insitu = gpsfh.SFH( tau_quench = tau_quench, model = 'insitu' )
    assert sfh_insitu.Mdust( 8.e+8 ) == pytest.approx( 136530596.03205943 ) 
    assert sfh_insitu.Mdust( 1.e+9 ) == 0.
    
#------------------------------------------------------------------------------#

def test_sfh_insitu_Mgas () :
    
    """
    Check the value of total gas content at quench and at present time .
    
    """
    
    sfh_insitu = gpsfh.SFH( tau_quench = tau_quench, model = 'insitu' )
    assert sfh_insitu.Mgas( 8.e+8 ) == pytest.approx( 5781720297.528179 ) 
    assert sfh_insitu.Mgas( 1.e+9 ) == 0.
    
#------------------------------------------------------------------------------#

def test_sfh_insitu_Zstar () :
    
    """
    Check the value of the stellar metallicity at present time ( ~0.02 ) .

    """
    sfh_insitu = gpsfh.SFH( tau_quench = tau_quench, model = 'insitu' )
    assert sfh_insitu.Zstar( 1.e+9 ) == pytest.approx( 0.023106097623840163 )
    
#------------------------------------------------------------------------------#

def test_sfh_insitu_Zgas () :
    
    """
    Check the value of the gas metallicity at present time ( ~0.034 ) .
    """
    sfh_insitu = gpsfh.SFH( tau_quench = tau_quench, model = 'insitu' )
    assert sfh_insitu.Zgas( 1.e+9 ) == pytest.approx( 0.03440966614222255 ) 

#------------------------------------------------------------------------------#

#Chosen values for Mdust and Zgs:

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

def test_sfh_const_init () :
    
    """
    Test initialization of the class SFH.
 
    """
    sfh_const = gpsfh.SFH( model ='constant' )
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

def test_sfh_const_call () :
    
    """ 
    Test calling the class at specific time ( array-like or float type ) .

    """
    
    sfh_const = gpsfh.SFH( tau_quench = tau_quench, model = 'constant' )
    sfh_const.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_const( 1.e+8 ) == 1.
    tau = np.logspace( 6, 10, 4 )
    assert np.all( sfh_const( tau ) == pytest.approx( [ 1., 1., 1., 0. ] ) )

#------------------------------------------------------------------------------#

def test_sfh_const_quench () :
    
    """
    After tau_quench, sfh must be equal to 0 .

    """
    
    sfh_const = gpsfh.SFH( tau_quench = tau_quench, model = 'constant' )
    sfh_const.set_parameters( Mdust = md, Zgs = zz )
    tau = np.logspace( 9, 11, 100 )
    psi_const = sfh_const( tau )
    assert all ( [ a == 0. for a in psi_const ] )
    
#------------------------------------------------------------------------------#

def test_sfh_const_Mstar () :
    
    """
    Check the value of the stellar mass at quench and at present time .

    """
    
    sfh_const = gpsfh.SFH( tau_quench = tau_quench, model = 'constant' )
    sfh_const.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_const.Mstar( 8.e+8 ) == pytest.approx( 535695178.04297197 ) 
    assert sfh_const.Mstar( 1.e+9 ) == pytest.approx( 509519015.14239407 ) 

#------------------------------------------------------------------------------#

def test_sfh_const_Mdust () :
    
    """
    Check the value of total dust content at quench and at present time
    ( must be equal to md ) .

    """
    sfh_const = gpsfh.SFH( tau_quench = tau_quench, model = 'constant' )
    sfh_const.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_const.Mdust( 8.e+8 ) == md 
    assert sfh_const.Mdust( 1.e+9 ) == md
    
#------------------------------------------------------------------------------#

def test_sfh_const_Mgas () :
    
    """
    Check the value of total gas content at present time .
    
    """
    
    sfh_const = gpsfh.SFH( tau_quench = tau_quench, model = 'constant' )
    sfh_const.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_const.Mgas( 1.e+9 ) == pytest.approx( 1405518018.9671905 )
    
#------------------------------------------------------------------------------#

def test_sfh_const_Zstar () :
    
    """
    Check the value of the stellar metallicity at present time 
    ( must be equal to zz ) .

    """
    sfh_const = gpsfh.SFH( tau_quench = tau_quench, model = 'constant' )
    sfh_const.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_const.Zstar( 1.e+9 ) == zz
    
#------------------------------------------------------------------------------#

def test_sfh_const_Zgas () :
    
    """
    Check the value of the Gas Metallicity at present time 
    ( must be equal to zz ) .
    """
    sfh_const = gpsfh.SFH( tau_quench = tau_quench, model = 'constant' )
    sfh_const.set_parameters( Mdust = md, Zgs = zz )
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

def test_sfh_dexp_init () :
    
    """
    Test initialization of the class SFH .
 
    """
    sfh_dexp = gpsfh.SFH( model ='delayedexp' )
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

def test_sfh_dexp_call () :
    
    """ 
    Test calling the class at specific time ( array-like or float type ) .

    """
    
    sfh_dexp = gpsfh.SFH( tau_quench = tau_quench, model = 'delayedexp' )
    sfh_dexp.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_dexp( 1.e+8 ) == 14.645544342956468
    tau = np.logspace( 6, 10, 4 )
    assert np.all( sfh_dexp( tau ) == pytest.approx( [ 15.69123242, 23.61025931,
                                                       0.52181543, 0. ] ) )

#------------------------------------------------------------------------------#

def test_sfh_dexp_quench () :
    
    """
    After tau_quench, sfh must be equal to 0 .

    """
    
    sfh_dexp = gpsfh.SFH( tau_quench = tau_quench, model = 'delayedexp' )
    sfh_dexp.set_parameters( Mdust = md, Zgs = zz )
    tau = np.logspace( 9, 11, 100 )
    psi_dexp = sfh_dexp( tau )
    assert all ( [ a == 0. for a in psi_dexp ] )
    
#------------------------------------------------------------------------------#

def test_sfh_dexp_Mstar () :
    
    """
    Check the value of the stellar mass at quench time and at present time .

    """
    
    sfh_dexp = gpsfh.SFH( tau_quench = tau_quench, model = 'delayedexp' )
    sfh_dexp.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_dexp.Mstar( 8.e+8 ) == pytest.approx( 2302675603.296462 ) 
    assert sfh_dexp.Mstar( 1.e+9 ) == pytest.approx( 2255427119.0862346 ) 

#------------------------------------------------------------------------------#

def test_sfh_dexp_Mdust () :
    
    """
    Check the value of total dust content at quench and at present time .
    ( must be equal to md ) .

    """
    sfh_dexp = gpsfh.SFH( tau_quench = tau_quench, model = 'delayedexp' )
    sfh_dexp.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_dexp.Mdust( 8.e+8 ) == md 
    assert sfh_dexp.Mdust( 1.e+9 ) == md
    
#------------------------------------------------------------------------------#

def test_sfh_dexp_Mgas () :
    
    """
    Check the value of total gas content at present time .
    
    """
    
    sfh_dexp = gpsfh.SFH( tau_quench = tau_quench, model = 'delayedexp' )
    sfh_dexp.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_dexp.Mgas( 1.e+9 ) == pytest.approx( 1405518018.9671905 )
    
#------------------------------------------------------------------------------#

def test_sfh_dexp_Zstar () :
    
    """
    Check the value of the stellar metallicity at present time 
    ( must be equal to zz ) .

    """
    sfh_dexp = gpsfh.SFH( tau_quench = tau_quench, model = 'delayedexp' )
    sfh_dexp.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_dexp.Zstar( 1.e+9 ) == zz
    
#------------------------------------------------------------------------------#

def test_sfh_dexp_Zgas () :
    
    """
    Check the value of the gas metallicity at present time 
    ( must be equal to zz ) .
    """
    sfh_dexp = gpsfh.SFH( tau_quench = tau_quench, model = 'delayedexp' )
    sfh_dexp.set_parameters( Mdust = md, Zgs = zz )
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

def test_sfh_lnorm_init () :
    
    """
    Test initialization of the class SFH .
 
    """
    sfh_lnorm = gpsfh.SFH( model ='lognormal' )
    assert isinstance( sfh_lnorm, galapy.StarFormationHistory.SFH )

#------------------------------------------------------------------------------#

def test_sfh_lnorm_parameters () :

    """
    Test default parameters for 'lognormal' model .

    """
    
    sfh_lnorm = gpsfh.SFH( model = 'lognormal' )
    print(sfh_lnorm.params)
    assert sfh_lnorm.params == { 'tau_quench': 2e+20, 'model': 'lognormal',
                                 'psi_norm': 100.0, 'sigma_star': 2.0,
                                 'tau_star': 300000000.0, 'Mdust': 100000000.0,
                                 'Zgs': 0.1 }


#------------------------------------------------------------------------------#

def test_sfh_lnorm_call () :
    
    """ 
    Test calling the class at specific time ( array-like or float type ) .

    """
    
    sfh_lnorm = gpsfh.SFH( tau_quench = tau_quench, model = 'lognormal' )
    sfh_lnorm.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_lnorm( 1.e+8 ) == 17.153733592792385
    tau = np.logspace( 6, 10, 4 )
    assert np.all( sfh_lnorm( tau ) == pytest.approx( [ 0.34179049, 8.38175995,
                                                        19.47777367, 0. ] ) )

#------------------------------------------------------------------------------#

def test_sfh_lnorm_quench () :
    
    """
    After tau_quench, sfh must be equal to 0 .

    """
    
    sfh_lnorm = gpsfh.SFH( tau_quench = tau_quench, model = 'lognormal' )
    sfh_lnorm.set_parameters( Mdust = md, Zgs = zz )
    tau = np.logspace( 9, 11, 100 )
    psi_lnorm = sfh_lnorm( tau )
    assert all ( [ a == 0. for a in psi_lnorm ] )
    
#------------------------------------------------------------------------------#

def test_sfh_lnorm_Mstar () :
    
    """
    Check the value of the stellar mass at quench and at present time .

    """
    
    sfh_lnorm = gpsfh.SFH( tau_quench = tau_quench, model = 'lognormal' )
    sfh_lnorm.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_lnorm.Mstar( 8.e+8 ) == pytest.approx( 9751556099.681376 ) 
    assert sfh_lnorm.Mstar( 1.e+9 ) == pytest.approx( 9273433696.49997 ) 

#------------------------------------------------------------------------------#

def test_sfh_lnorm_Mdust () :
    
    """
    Check the value of total dust content at quench time and at present time 
    ( must be equal to md ) .

    """
    sfh_lnorm = gpsfh.SFH( tau_quench = tau_quench, model = 'lognormal' )
    sfh_lnorm.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_lnorm.Mdust( 8.e+8 ) == md 
    assert sfh_lnorm.Mdust( 1.e+9 ) == md
    
#------------------------------------------------------------------------------#

def test_sfh_lnorm_Mgas () :
    
    """
    Check the value of total gas content at present time .
    
    """
    
    sfh_lnorm = gpsfh.SFH( tau_quench = tau_quench, model = 'lognormal' )
    sfh_lnorm.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_lnorm.Mgas( 1.e+9 ) == pytest.approx( 1405518018.9671905 )
    
#------------------------------------------------------------------------------#

def test_sfh_lnorm_Zstar () :
    
    """
    Check the value of the stellar metallicity at present time 
    ( must be equal to zz ) .

    """
    sfh_lnorm = gpsfh.SFH( tau_quench = tau_quench, model = 'lognormal' )
    sfh_lnorm.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_lnorm.Zstar( 1.e+9 ) == zz
    
#------------------------------------------------------------------------------#

def test_sfh_lnorm_Zgas () :
    
    """
    Check the value of the Gas Metallicity at present time 
    ( must be equal to zz ) .
    """
    sfh_lnorm = gpsfh.SFH( tau_quench = tau_quench, model = 'lognormal' )
    sfh_lnorm.set_parameters( Mdust = md, Zgs = zz )
    assert sfh_lnorm.Zgas( 1.e+9 ) == zz 

#------------------------------------------------------------------------------#
