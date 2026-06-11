######################################################################################
# External imports

import numpy

######################################################################################
# Internal imports

from galapy.internal.utils import (
    quantile_weighted,
    get_credible_interval,
    find_nearest,
)

######################################################################################
# General posterior estimators
# These are stateless functions that operate on raw (values, weights) arrays,
# generalising the class methods of galapy.sampling.Results.GalaxyResults.

def weighted_quantile ( values, quantiles, weights = None ) :
    """ Weighted quantile of an array.

    Parameters
    ----------
    values : array-like, shape (N,) or (N, M)
    quantiles : scalar or array-like, values in [0, 1]
    weights : array-like, shape (N,), optional

    Returns
    -------
    numpy.ndarray
    """
    return quantile_weighted( values, quantiles, weights = weights )

def weighted_mean ( values, weights = None ) :
    """ Weighted mean along axis 0.

    Parameters
    ----------
    values : array-like, shape (N,) or (N, M)
    weights : array-like, shape (N,), optional

    Returns
    -------
    scalar or numpy.ndarray
    """
    if weights is None :
        return numpy.mean( values, axis = 0 )
    w = numpy.asarray( weights )
    mask = w > 0
    return numpy.average( numpy.asarray( values )[mask],
                          weights = w[mask], axis = 0 )

def weighted_std ( values, weights = None ) :
    """ Weighted standard deviation along axis 0.

    Parameters
    ----------
    values : array-like, shape (N,) or (N, M)
    weights : array-like, shape (N,), optional

    Returns
    -------
    scalar or numpy.ndarray
    """
    mean = weighted_mean( values, weights )
    v    = numpy.asarray( values )
    if weights is None :
        return numpy.std( v, axis = 0 )
    w    = numpy.asarray( weights )
    mask = w > 0
    return numpy.sqrt(
        numpy.average( ( v[mask] - mean ) ** 2,
                       weights = w[mask], axis = 0 )
    )

def credible_interval ( values, logl, percent = 0.68,
                        centre = 'bestfit', weights = None ) :
    """ Credible interval for a 1-D marginal posterior.

    Parameters
    ----------
    values : array-like, shape (N,)
        Samples of the quantity of interest.
    logl : array-like, shape (N,)
        Log-likelihood values corresponding to each sample.
    percent : float, optional
        Probability mass enclosed by the interval (default 0.68).
    centre : str, optional
        'bestfit' (default), 'mean', or 'median'.
    weights : array-like, shape (N,), optional
        Sample weights. Uniform if None.

    Returns
    -------
    low, upp : float
        Distances from the centre to the lower and upper bounds.
        Returns (-inf, upp) for an upper limit and (low, +inf) for a
        lower limit when the distribution is one-sided.
    """
    v = numpy.asarray( values )
    if centre == 'bestfit' :
        idcentre = numpy.asarray( logl ).argmax()
    elif centre == 'mean' :
        idcentre = find_nearest( v, weighted_mean( v, weights ) )
    elif centre == 'median' :
        idcentre = find_nearest( v, weighted_quantile( v, 0.5, weights ) )
    else :
        raise ValueError(
            "``centre`` must be one of ('bestfit', 'mean', 'median')"
        )
    low, upp = get_credible_interval( v, idcentre, percent, weights )
    if low is None :
        return -numpy.inf, upp
    if upp is None :
        return low, +numpy.inf
    return v[idcentre] - low, upp - v[idcentre]

######################################################################################
# Bayesian model comparison

def bayes_factor ( logz1, logz2 ) :
    """ Log Bayes factor B_12 = exp(logz1 - logz2).

    A positive value favours model 1; a negative value favours model 2.

    Parameters
    ----------
    logz1, logz2 : float
        Natural log-evidences of the two models.

    Returns
    -------
    float
        ln(B_12) = logz1 - logz2.
    """
    return logz1 - logz2

def jeffreys_scale ( log_bf ) :
    """ Classify a log Bayes factor using the Jeffreys (1961) scale.

    Parameters
    ----------
    log_bf : float
        ln(B_12) as returned by :func:`bayes_factor`.

    Returns
    -------
    str
        Qualitative strength of evidence in favour of the model with the
        higher evidence (positive log_bf → model 1; negative → model 2).
    """
    abf = abs( log_bf )
    if abf < 1.0 :
        label = 'Inconclusive'
    elif abf < 2.5 :
        label = 'Moderate'
    elif abf < 5.0 :
        label = 'Strong'
    else :
        label = 'Decisive'
    favoured = 1 if log_bf >= 0 else 2
    return f'{label} evidence in favour of model {favoured}'
