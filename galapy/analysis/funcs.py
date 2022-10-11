######################################################################################
# External imports

import numpy

######################################################################################
# Internal imports

from galapy.sampling.Results import Results

######################################################################################

def get_parameters_summary_statistics ( res, stat_type = 'quantiles', quantile = None ) :
    """ Function returning a dictionary with summary statistics of the samples in 
    the Results instance of some sampling run
    
    Parameters
    ----------
    res : galapy.sampling.Results.Results
        An instance of type ``Results``
    stat_type : str
        The desired summary statistics. Available statistics are
        'mean_and_std' and 'quantiles'
    quantile : scalar or sequence
        The quantiles requested (used only is ``stat_type='quantiles'``.
        Default is (0.16,0.5,0.84), i.e. median and 68% bound.
    
    Returns
    -------
    : ndarray
        - if ``stat_type='quantiles'`` each row is an array with size = ``len(quantile)``
        - if ``stat_type='mean_and_std'`` each row is an array of size = 2 where the
          first element is the mean and the second element is the standard deviation.
    """
    
    if stat_type not in ('mean_and_std', 'quantiles') :
        raise ValueError(
            "Allowed values for argument ``stat_type`` are 'mean&std' and 'quantiles'"
        )
        
    if stat_type == 'quantiles' :
        if quantile is None or len(quantile)!=3:
            print( 'Falling back to default quantiles: (0.16,0.5,0.84)')
            quantile = [0.16, 0.5, 0.84]
        summary = res.get_quantile('samples', quantile).T
        summary[:,0] = abs(summary[:,0]-summary[:,1])
        summary[:,2] = abs(summary[:,2]-summary[:,1])
        
    if stat_type == 'mean_and_std' :
        summary = numpy.array([res.get_mean('samples'), 
                               res.get_std('samples')]).T
    return summary

######################################################################################

def get_parameters_summary_strings ( res, digits = 2,
                                     stat_type = 'quantiles', quantile = None ) :
    """ Function returning a list of formatted strings with the chosen summary 
    statistics of the samples in the Results instance of some sampling run
    
    Parameters
    ----------
    res : galapy.sampling.Results.Results
        An instance of type ``Results``
    digits : integer
        Number of significant digits to show.
        Each statistic will be represented in strings formatted as ``{:.{digits}f}``
        Default is 2
    stat_type : str
        The desired summary statistics. Available statistics are
        'mean_and_std' and 'quantiles'
    quantile : scalar or sequence
        The quantiles requested (used only is ``stat_type='quantiles'``.
        Default is (0.16,0.5,0.84), i.e. median and 68% bound.
    
    Returns
    -------
    : ndarray
        An array with Ndim TeX-formatted summary statistics, where Ndim is the
        size of the free-parameters space of the sampling run.
    """
    
    summary = get_parameters_summary_statistics( res, stat_type, quantile )
    fstring = None
    if stat_type == 'quantiles' :
        fstring = f'{{1:.{digits}f}}_{{{{-{{0:.{digits}f}}}}}}^{{{{+{{2:.{digits}f}}}}}}'
    if stat_type == 'mean_and_std' :
        fstring = f'{{0:.{digits}f}} \\pm {{1:.{digits}f}}'
    return numpy.array([fstring.format(*s) for s in summary])

######################################################################################
