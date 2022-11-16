import numpy
from scipy.special import erf

#############################################################################################

def simple_uplim ( data, error, model ) :
    data = numpy.asarray( data )
    error = numpy.asarray( error )
    model = numpy.asarray( model )
    
    if numpy.any( model > data ) :
        return -numpy.inf
    return 0.

#############################################################################################

# For upper limits (using eq. A10 of Sawicki2012):
def chi2_uplim ( data, error, model ) :
    
    # get argument of logarithm
    value = ( numpy.sqrt( 0.5 * numpy.pi ) *
              error *
              ( 1 + erf( numpy.sqrt( 0.5 ) *
                         ( data - model ) /
                         error ) ) )
    
    # check if there are values less than/equal to zero
    if numpy.any( ~( value > 0. ) ) :
        return numpy.inf

    # compute logarithm
    chi = numpy.log( value )

    return -2 * numpy.sum( chi )

#############################################################################################

def chi2 ( data, error, model ) :
    d, e, m = numpy.asarray(data), numpy.asarray(error), numpy.asarray(model)
    chi = ( d - m ) / e
    return numpy.sum( chi * chi )

#############################################################################################

f_uplims = { 'simple' : simple_uplim, 'chi2' : chi2, 'S12' : chi2_uplim }

#############################################################################################

def gaussian_loglikelihood ( data, error, model, uplims, method_uplims = 'simple' ) :
    """
    Uplim methods are:
    - 'simple' : 
    - 'chi2' :
    - 'S12' :
    """

    try :
        uplims = numpy.asarray( uplims )
    except KeyError as ke :
        raise RuntimeError( f"the method '{method_uplim}' is not valid, valid methods are {list(f_uplim.keys())}" )
    c2 = chi2( data[~uplims], error[~uplims], model[~uplims] )
    u2 = f_uplims[method_uplims]( data[uplims], error[uplims], model[uplims] )
    if numpy.isnan( c2 ) :
        return -numpy.inf
    if numpy.isnan( u2 ) :
        return -numpy.inf
        
    return -0.5 * ( c2 + u2 )

#############################################################################################

def transform_to_prior_unit_cube ( values, prior_limits ) :
    """ Maps the values of the parameters in the unit-cube prior, considered as uniform.
    
    Parameters
    ----------
    values : scalar or array-like
      value(s) of the parameters 
    prior_limits : array-like
      array with the lower/upper bounds of the uniform prior
      associated to each parameter listed in 'values'. 
      Shape should be Nx2 where N == len(values).

    Returns
    -------
    : scalar or array-like
    values mapped within the limits of the prior (all values in the [0.,1.) interval)
    """
    
    x = numpy.asarray( values ) 

    pmin, pmax = numpy.asarray( prior_limits ).T
    delta = pmax-pmin
    x *= delta
    x += pmin
        
    return x    

#############################################################################################
