########################################################################################

# External imports
from collections.abc import MutableMapping
import os

# import all from SADFit library
# import sadfit
# from .internal.utils import SmartAttribute, unwrap_values, unwrap_keys, unwrap_items

########################################################################################
"""
validate_* functions check that the input variable is valid.
If true they return it as it is
or transformed in a convenient form.
If False raise a ValueError
"""
def validate_time ( t ) :
    return t

def validate_string ( s ) :
    if isinstance( s, str ) :
        return s
    else:
        raise ValueError( f'"{s}" of type {type(s)} should be a string' )

def validate_path_exists ( s ) :
    """If s is a path, return s, else raise RuntimeError"""
    
    if isinstance( s, str ) and os.path.exists(s):
        return s
    
    else:
        raise ValueError( f'"{s}" of type {type(s)} should be an existing path but is not' )

def validate_model ( m ) :

    valid = [ 'insitu', 'constant', 'delayedexp', 'lognormal', 'burst' ]
    if m in valid :
        return m
    else :
        raise ValueError( f'Value "{m}" should be a valid model but is not.\n'
                          f'Valid models are {valid}.' )

def validate_float ( f ) :
    if isinstance( f, float ) :
        return f
    else :
        raise ValueError( f'"{f}" should be a float but has type={type(f)}' )

def validate_float_positive ( f ) :
    if isinstance( f, float ) and not f < 0. :
        return f
    else :
        raise ValueError( f'"{f}" of type {type(f)} should be a positive float but is not' )

def validate_time ( t ) :
    return validate_float_positive( t )

def validate_fraction ( f ) :
    if isinstance( f, float ) and not f < 0. and not f > 1. :
        return f
    else :
        raise ValueError( f'"{f}" of type {type(f)} should be a float in [0,1] but is not' )

########################################################################################
"""
str_* functions define how a given parameter is transformed into string.
When referring to physical quantities, 
they provide human-readable access to the assumed units.
"""

def str_time ( t ) :
    return '{:.5e} [ yr ]'.format( t )

def str_model ( m ) :
    sfh_models = {
        1 : 'In-Situ',
        2 : 'Constant',
        3 : 'Delayed-Exponential',
        4 : 'Log-Normal'
    }
    return sfh_models[ m ]

def str_mass ( m ) :
    return '{:.5e} [ Msol ]'.format( m )

def str_sfr ( s ) :
    return '{:.5e} [ Msol * yr^-1 ]'.format( s )

def str_lenght ( l ) :
    return '{:.5e} [ pc ]'.format( l )

########################################################################################
"""
Dictionary containing default values of parameters
"""
    
# sfh_model_choice = {
#     'insitu' : {
#         'psi_max' :  SmartAttribute( value = 100.,
#                                      valid = validate_float_positive,
#                                      string = str_sfr ),
#         'tau_star' : SmartAttribute( value = 3.e+8,
#                                      valid = validate_time,
#                                      string = str_time )
#     },
#     # Constant SF model
#     'constant' : {
#         'psi' : SmartAttribute( value = 1.,
#                                 valid = validate_float_positive,
#                                 string = str_sfr )
#     },
#     # Delayed-Exponential model
#     'delayedexp' : {
#         'psi_norm' : SmartAttribute( value = 1.,
#                                      valid = validate_float_positive,
#                                      string = str_sfr ),
#         'k_shape' :  SmartAttribute( value = 0.2,
#                                      valid = validate_float ),
#         'tau_star' : SmartAttribute( value = 1.e+8,
#                                      valid = validate_time,
#                                      string = str_time )
#     },
#     # Log-Normal model
#     'lognormal' : {
#         'psi_norm' :   SmartAttribute( value = 100.,
#                                        valid = validate_float_positive,
#                                        string = str_sfr ),
#         'sigma_star' : SmartAttribute( value = 2.,
#                                        valid = validate_float ),
#         'tau_star' :   SmartAttribute( value = 3.e+8,
#                                        valid = validate_time,
#                                        string = str_time )
#     }
# }

# def _build_param_dict ( **kwargs ) :
#     out = {
#         # Age of the galaxy
#         'age' : SmartAttribute( value = 1.e+8,
#                                 valid = validate_time,
#                                 string = str_time ),

#         # Simple-Stellar-Population library path
#         'csp' : {
#             'ssp_library' : SmartAttribute( value = '/path/to/nowhere',
#                                             valid = validate_path_exists ),
#         },

#         # Star Formation History parameters
#         'sfh' : {
#             'tau_quench' : SmartAttribute( value = 1.e+11,
#                                            valid = validate_time,
#                                            string = str_time ),
#             'model' :      SmartAttribute( value = 'insitu',
#                                            valid = validate_string ),
#         },

#         # Inter-Stellar Medium parameters
#         'ism' : {
#             'f_MC' : SmartAttribute( value = 0.5,
#                                      valid = validate_fraction ),
#             # Diffuse-Dust phase
#             'dd' : {
#                 'norm' :  SmartAttribute( value = 1.,
#                                           valid = validate_float_positive ),
#                 'Mdust' : SmartAttribute( value = 1.e+7,
#                                           valid = validate_float_positive,
#                                           string = str_mass ),
#                 'Rdust' : SmartAttribute( value = 1000.,
#                                           valid = validate_float_positive,
#                                           string = str_lenght  ),
#                 'f_PAH' : SmartAttribute( value = 0.,
#                                           valid = validate_fraction ),
#             },
#             # Molecular-Cloud phase
#             'mc' : {
#                 'norm' :    SmartAttribute( value = 100.,
#                                             valid = validate_float_positive ),
#                 'M_MC' :    SmartAttribute( value = 1.e+6,
#                                             valid = validate_float_positive,
#                                             string = str_mass  ),
#                 'R_MC' :    SmartAttribute( value = 16,
#                                             valid = validate_float_positive,
#                                             string = str_lenght  ),
#                 'Zgas' :    SmartAttribute( value = 0.5,
#                                             valid = validate_float ),
#                 'tau_esc' : SmartAttribute( value = 1.e+8,
#                                             valid = validate_time,
#                                             string = str_time  ),
#                 'Mgas' :    SmartAttribute( value = 1.e+8,
#                                             valid = validate_float_positive,
#                                             string = str_mass  ),
#             },
#         },

#         # # Active Galactic Nucleus parameters
#         # 'agn' : {
#         #     'ct' : SmartAttribute( value = ,
#         #                            valid = ,
#         #                            string = ),
#         #     'al' : SmartAttribute( value = ,
#         #                            valid = ,
#         #                            string = ),
#         #     'be' : SmartAttribute( value = ,
#         #                            valid = ,
#         #                            string = ),
#         #     'ta' : SmartAttribute( value = ,
#         #                            valid = ,
#         #                            string = ),
#         #     'rm' : SmartAttribute( value = ,
#         #                            valid = ,
#         #                            string = ),
#         #     'ia' : SmartAttribute( value = ,
#         #                            valid = ,
#         #                            string = ),

#         # },
#     }
    
#     out['sfh'].update( sfh_model_choice[ out['sfh']['model'].value ] )
    
#     return out
