""" Plotting utilities for the analysis of the results and to generate easily plots 
from GalaPy structures and functions.
"""
######################################################################################
# External imports
import os
import warnings
import fnmatch   

# As the sane of mind do
import numpy

# For loading stored files
# import pickle

# For plotting
import matplotlib.pyplot as plt

# Matplotlib setup of the GalaPy default layout
plt.rcParams['figure.figsize'] = (7,5)
plt.rcParams['font.size'] = 14
plt.rcParams['lines.linewidth'] = 2.
#plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.style.use('seaborn-deep')
clr = plt.rcParams['axes.prop_cycle'].by_key()['color']

######################################################################################
# Internal imports

import galapy as gp
from galapy.GalaxyParameters import gxy_params_defaults as gpdefault
from galapy.sampling.Results import Results
from galapy.analysis.funcs import get_parameters_summary_strings, get_parameters_label_strings

######################################################################################
# Basic functions

# a function for formatting the axes' ticks:
def format_axes_ticks(fig, labelsize=14):
    for i, ax in enumerate(fig.axes):
        ax.tick_params(labelsize=labelsize)

######################################################################################
# Default values of the plots settings

_plot_specs_default = {
    'ax_kw' : {
        'xscale' : 'log', 
        'yscale' : 'log',
        'xlim' : (1.e+0,1.e+10),
        'ylim' : (1.e-10, 1.e+10),
        'xlabel' : '$\\lambda\\ [\\,\\AA\\,]$',
        'ylabel' : '$S_\\lambda\\ [\\mathrm{mJy}]$',
    },
    'getdist_settings' : { 
        'mult_bias_correction_order' : 1,
        'smooth_scale_2D':4, 
        'smooth_scale_1D':4, 
        'fine_bins': 128, 'fine_bins_2D' : 128,
        'contours' : [0.68, 0.95, 0.99],
    },
    'titles_kw' : { 
        'digits' : 2, 
        'stat_type' : 'quantiles', 
        'quantile' : None 
    },
    'marker_kw' : { 'ls' : '--', 'lw' : 1.5, 'color' : 'black' },
    'triangle_kw' : { 
        'contour_colors' : ['gray'],
        'line_args' : [{'lw':2, 'color':'gray'}],
    },
}

def show_default_dict ( name = None ) :
    """ Shows the default dictionary for the kwargs dictionary with name 'name'.
    """
    if name is not None :
        return _plot_specs_default[ name ]
    return _plot_specs_default

######################################################################################
# Functions for plotting SEDs

def sed_plot_layout ( redshift, frame, ax = None, **kwargs ) :
    """ Produces/modifies a matplotlib.axes.Axes instance with the grid a labels for
    plotting SED fluxes.
    
    Parameters
    ----------
    redshift : scalar
      the redshift of the observed frame
    frame : str
      one among ('obs', 'rest', 'both'), whether to scale the x-axis
      to the observed frame, the rest-frame or both.
      In the latter case a secondary x-axis will be added on top of the grid
    ax : matplotlib.axes.Axes
      (Optional) 
    kwargs : dictionary
      (Optional) the eventual keyword arguments to pass to the ```set``` function
      of ```ax```.
    
    Returns
    -------
    : matplotlib.axes.Axes
      the generated axes
    """
    
    # If not passed, attach to the latest axis or generate new one
    if ax is None :
        ax = plt.gca()
    if hasattr(ax, 'sed_layout_set' ) :
        if ax.sed_layout_set :
            return(ax)
    
    # Add secondary axis REST-FRAME WAVELENGHT
    if frame == 'both' :
        def red(x):
            return x * (1+redshift)
        def blue(x):
            return x / (1+redshift)
        secax = ax.secondary_xaxis('top', functions=(blue, red))
        _ = secax.set_xlabel('$\\lambda_\\mathrm{rest}\\ [\\,\\AA\\,]$')
        frame = 'obs'
        
    # Set axis parameters
    default_ax_kw = _plot_specs_default['ax_kw']
    default_ax_kw['xlabel'] = f'$\\lambda_\\mathrm{{{frame}}}\\ [\\,\\AA\\,]$'
    default_ax_kw.update(kwargs)
    _ = ax.set(**default_ax_kw)

    # Set background colors
    ax.set_facecolor('white')
    
    # Set an attribute to check the method has already been called
    try :
        _ = setattr(ax, 'sed_layout_set', True)
    except :
        raise Error('Something went wrong.')
    
    return(ax)

def plot_sed_obs ( ll, ff, ee, lo, redshift = None, frame = 'rest', ax = None, ax_kwargs = {} ) :

    ax = sed_plot_layout(redshift, frame, ax, **ax_kwargs)
    
    Pdata = ax.errorbar(ll, ff, ee, uplims = lo, 
                        label='data', 
                        ls='none', marker='o', 
                        markerfacecolor='white',
                        markeredgewidth=2.,
                        markersize=8)
    
    return Pdata

def plot_sed_components ( ll, components,
                          redshift = None, frame = 'rest', ax = None, ax_kwargs = {} ) :
    ax = sed_plot_layout(redshift, frame, ax, **ax_kwargs)
    Pcomp = []
    for k,L in components.items() :
        Pcomp.append(ax.plot(ll, L, '-', lw = 1.5, label=k)[0])
    return Pcomp

def plot_sed_flux ( ll, flux,
                    color = 'black', label = 'total', 
                    redshift = None, frame = 'rest', ax = None, ax_kwargs = {} ) :
    ax = sed_plot_layout(redshift, frame, ax, **ax_kwargs)
    return ax.plot( ll, flux, ls='-', color=color, lw=1.5, label=label)[0]

def plot_sed_1sigma2sigma ( ll, flux, err, center = False, color = 'gray', center_label = 'mean', 
                            redshift = None, frame = 'rest', ax = None, ax_kwargs = {} ) :
    ax = sed_plot_layout( redshift, frame, ax, **ax_kwargs )
    handles = []
    if center :
        handles.append( ax.plot( ll, flux,
                                 ls='-', color=color,
                                 lw=1.5, label=center_label )[0] )
    handles.append(
        ax.fill_between(ll, flux-err, flux+err,
                        label='$1$-$\\sigma$',
                        edgecolor='none', color=color, alpha=0.5)
    )
    handles.append(
        ax.fill_between(ll, flux-2*err, flux+2*err,
                        label='$2$-$\\sigma$',
                        edgecolor='none', color=color, alpha=0.2)
    )
    return handles

def plot_sed_flux_res ( res,
                        model = None,
                        observation = None,
                        plot_observation = False,
                        plot_components = False, 
                        plot_contours = True,
                        frame = 'both', 
                        ax = None, 
                        ax_kwargs = {} ) :
    """ Plots the flux
    
    Parameters
    ----------
    res : Results instance
        A ``Results`` instance from a sampling run.
    model : GXY instance
        (Optional) A galaxy model or any class derived from ``galapy.Galaxy.GXY``.
    observation : Observation instance
        (Optional) an ``Observation`` instance. Note that if ``plot_observation=False``
        it will be ignored anyways.
    plot_observation : bool
        Whether to plot the observational data used in the sampling (dafault is ``False``).
    plot_components : bool
        Whether to plot the different components contributing to the best-fit SED 
        (dafault is ``False``).
    plot_contours : bool
        Whether to plot the 1- and 2-sigma contours around the mean SED (dafault is ``False``).
    frame : string
        One among {``'rest'``, ``'obs'``, ``'both'``}, choose the frame, observed- or rest-frame,
        to plot the wavelenght-axis (default is ``'both'``)
    ax : matplotlib.axes.Axes
        (Optional) an instance of matplotlib axes
    ax_kwargs : dict
        (Optional) keyword arguments to pass to the function ``sed_plot_layout``
    
    Returns
    -------
    : matplotlib.axes.Axes
        An instance of ``matplotlib.axes.Axes`` with the plotted fluxes
    """
    
    if not isinstance( res, Results ) :
        raise AttributeError( 
            'Attribute "res" should be an instance of type ``Results``'
        )
        
    if model is None :
        model = res.get_model()
    else :
        if not isinstance( model, gp.Galaxy.GXY ) :
            raise AttributeError(
                'Attribute "model" should be an instance of type ``GXY``'
            )
    if plot_observation :
        if observation is None :
            observation = res.get_observation()
        else :
            if not isinstance( observation,
                               gp.sampling.Observation.Observation ) :
                raise AttributeError(
                    'Attribute "observation" should be an instance of type ``Observation``'
                )

    # Set the wavelenght (lambda) axis:
    # If 'redshift' is a free parameter the local variable
    # will be set to the best-fit value, otherwise the
    # spectroscopic redshift will be extracted from the model.
    redshift = None
    lx = None
    obs_scale = frame in set(['obs', 'both'])
    if 'redshift' in res.get_sampling_params().par_free :
        redshift = res.get_bestfit( 'params' )['redshift']
        lx = model.wl()
        if obs_scale :
            lx *= ( 1 + redshift )
    else :
        redshift = model.redshift
        lx = model.wl( obs = obs_scale )
    
    # Set the layout
    ax = sed_plot_layout( redshift, frame, ax = ax, **ax_kwargs )
    
    ##################################################################
    # Here the plots start
    
    legend_primary   = []
    legend_secondary = []
    
    # Plot the observational dataset
    if plot_observation :
        ll, ff, ee, lo = ( observation.pms.lpiv,
                           observation.fluxes,
                           observation.errors,
                           observation.uplims )
        ff[lo] = 3*ee[lo]
        legend_primary.append(
            plot_sed_obs( ll, ff, ee, lo,
                          redshift = redshift,
                          frame = frame,
                          ax=ax )
        )
            
    # Plot the different components
    flux_bestf = None
    if plot_components :
        bestfit_p = res.get_bestfit( 'params' )
        model.set_parameters(**bestfit_p)
        flux_bestf = model.get_SED()
        components = model.components_to_flux()
        legend_secondary += plot_sed_components( lx, components,
                                                 redshift = redshift,
                                                 frame = frame,
                                                 ax=ax )
        
    if flux_bestf is None :
        flux_bestf = res.get_bestfit( 'SED' )
    legend_primary.append( plot_sed_flux( lx, flux_bestf, label = 'best-fit',
                                          redshift = redshift,
                                          frame = frame,
                                          ax=ax ) )
    
    if plot_contours :
        flux_mean = res.get_mean( 'SED' )
        flux_err  = res.get_std( 'SED' )
        legend_primary += plot_sed_1sigma2sigma( lx, flux_mean, flux_err,
                                                 redshift = redshift,
                                                 frame = frame,
                                                 ax=ax )
    
    ##################################################################
    # Finalise and return
    
    # LEGEND:
    lpri = ax.legend(handles = legend_primary, loc = 'upper right', fontsize=12)
    ax.add_artist(lpri)
    if len(legend_secondary)>0 :
        lsec = ax.legend(handles = legend_secondary, 
                         loc = 'upper left', ncol = 2, fontsize=12)
        ax.add_artist(lsec)
    
    return(ax)

######################################################################################

def plot_corner_res ( res, handler = None, which_params = None, getdist_settings = None,
                      param_limits = 'auto', plot_titles = True, mark = 'bestfit',
                      titles_kw = {}, triangle_kw = {}, marker_kw = {} ) :
    """ Plot the triangle plot with the 2D posteriors of a sampling run and 1D marginals on the diagonal.
    
    Parameters
    ----------
    res :  Results instance
        A ``Results`` instance from a sampling run.
    handler : GXYParameters instance
        (Optional) The ``GXYParameters`` corresponding to the given sampling run
    which_params : str or sequence of str
        Either a single string or a sequence of strings. Name of the parameters to show on the 
        triangle plot. Also accepts wildcards (e.g. ``which_params = 'sfh*'`` will show all the
        parameters that contain the sub-string ``sfh``).
    getdist_settings : dict
        a dictionary of custom analysis settings to pass to the ``getdist.MCSamples`` settings argument 
        (for further informations see the documentation of ``getdist``)
        (to see defaults call ``galapy.analysis.plots.show_default_dict('getdist_settings')``)
    param_limits : str or sequence or dict
        If a string is passed it has to be one among
        - 'auto' : set the axes limits automatically
        - 'prior' : set the axes to the limits of the prior
        Otherwise it can be a 
        - 2D sequence (list or tuple) with dimensions (ndim, 2), where
          ndim is the number of free parameters of the sampling run. Each of the
          ndim couples will be assigned to the ordered list of free-parameters.
        - a dictionary where the keys are the names of the free parameters of the
          sampling run and the values are sequences of lenght = 2 with the limits
          for the corresponding free parameter.
    plot_titles : bool
        Whether to plot titles above the diagonal marginal posteriors with summary
        statistics computed from the corresponding posteriors (default summary is 
        the median with 68% uncertainty, to change this behaviour modify the 
        argument ``titles_kw``)
    mark : str
        What position in the posteriors space to highlight with lines, 
        default is the 'bestfit', other available positions are 
        'mean' and 'median'.
    titles_kw : dict
        dictionary to modify the titles above the diagonal marginal posteriors.
        It will be passed to function ``galapy.analysis.funcs.get_parameters_summary_strings()``
        (to see defaults call ``galapy.analysis.plots.show_default_dict('titles_kw')``)
    triangle_kw : dict
        Dictionary to modify the aspect of the triangle plot, these are the
        keyword arguments passed to ``triangle_plot()`` function of ``getdist``
        (to see defaults call ``galapy.analysis.plots.show_default_dict('triangle_kw')``)
    marker_kw : dict
        Dictionary to modify the aspect of the markers (color, linewidth, ...)
        (to see defaults call ``galapy.analysis.plots.show_default_dict('marker_kw')``)
    
    Returns
    -------
    : matplotlib.pyplot.Figure instance
    : matplotlib.axes.Axes instance
    """
    from getdist import plots, MCSamples
    import getdist

    ############################################################################
    # Set default kwarg-dictionaries
    
    default_getdist_settings = _plot_specs_default['getdist_settings']
    default_titles_kw = _plot_specs_default['titles_kw']
    default_marker_kw = _plot_specs_default['marker_kw']
    default_triangle_kw = _plot_specs_default['triangle_kw']
    
    ############################################################################
    # Check inputs
    
    if not isinstance( res, Results ) :
        raise AttributeError( 
            'Attribute "res" should be an instance of type ``Results``'
        )
    
    if handler is None :
        handler = res.get_sampling_params()
    else :
        if not isinstance( handler, gp.GalaxyParameters.GXYParameters ) :
            raise AttributeError(
                'Attribute "handler" should be an instance of type ``GXYParameters``'
            )
    
    if which_params is None :
        which_params = handler.par_free
    else :
        if isinstance(which_params, str) :
            which_params = fnmatch.filter(handler.par_free, which_params)
        elif isinstance(which_params, (list,tuple)) :
            which_params = [item for sublist in 
                            [fnmatch.filter(handler.par_free, s) 
                             for s in which_params] 
                            for item in sublist]
        else :
            raise TypeError(
                "Argument ``which_params`` must be a string, a list of strings or tuple of strings"
            )
    
    if getdist_settings is None :
        getdist_settings = {}
    default_getdist_settings.update(getdist_settings)
    
    if isinstance(param_limits, str) :
        if param_limits not in ('auto', 'prior') :
            raise ValueError(
                "Allowed string-values for argument ``param_limits`` are 'auto' and 'prior'"
            )
    elif isinstance(param_limits, MM) :
        pass
    elif isinstance(param_limits, (list, tuple)) :
        param_limits = dict(zip(handler.par_free, param_limits))
    else :
        raise AttributeError(
            "Argument ``param_limits`` must be a string or a dictionary"
        )

    markers = {}
    if mark is not None :
        if mark == 'bestfit' :
            markers = dict(zip(handler.par_free, res.get_bestfit('samples')))
        elif mark == 'mean' :
            markers = dict(zip(handler.par_free, res.get_mean('samples')))
        elif mark == 'median' :
            markers = dict(zip(handler.par_free, res.get_quantile('samples')))
        else :
            warnings.warn( "Argument ``mark`` should be one among 'bestfit', 'median', 'mean'." )
    
    # Update keyword arguments dictionaries
    default_triangle_kw.update(triangle_kw)
    default_titles_kw.update(titles_kw)
    default_marker_kw.update(marker_kw)
    
    ############################################################################
    # Prepare samples for plotting
    
    # Prepare dictionary with Prior-limits
    ranges = dict(zip(handler.par_free, handler.par_prior))
    
    # Prepare list with labels
    labels = get_parameters_label_strings( handler )
    
    # Extract the sampler's name from results
    sampler = None
    if res.sampler == 'dynesty' :
        sampler = 'nested'
    if res.sampler == 'emcee' :
        sampler = 'mcmc'
    
    # Build Samples-object
    mcsamples = MCSamples(
        samples = res.samples,
        loglikes = -res.logl,
        weights = res.weights,
        ranges = ranges,
        settings = default_getdist_settings,
        names = handler.par_free, 
        sampler = sampler,
        labels = list(labels.values())
    )
    
    # Prepare dictionary with parameters limits
    if param_limits == 'auto' : param_limits = {}
    if param_limits == 'prior' : param_limits = ranges 
     
    ############################################################################
    # Start plotting
    
    # Make triangle plot
    g = plots.get_subplot_plotter()
    g.settings.axis_tick_max_labels = 3
    g.settings.num_plot_contours = len(default_getdist_settings)
    g.triangle_plot(
        [mcsamples], params=which_params, 
        filled=True, title_limit=0, 
        param_limits = param_limits,
        markers = markers,
        marker_args = default_marker_kw,
        **default_triangle_kw
    )

    # Extract matplotlib objects
    fig, axes = g.fig, g.subplots
    
    # Set marker-args also on marginals:
    if len(markers) > 0 :
        for i in range(len(axes)) :
            axes[i,i].get_lines()[1].set(**default_marker_kw)
    
    # Plot formatted titles if requested
    if plot_titles :
        title_dict = dict(zip(
            handler.par_free,
            get_parameters_summary_strings( res, **default_titles_kw )
        ))
        for i, (ax, key) in enumerate(zip(axes, which_params)) :
            ax[i].set_title( ' '.join(('$', labels[key], '=', title_dict[key], '$')),
                             fontsize = g.settings.axes_fontsize )

    # Axes-ticks formatting
    format_axes_ticks( fig, labelsize = g.settings.axes_fontsize )

    ############################################################################
    # done and dusted.
    
    return fig, axes

######################################################################################
