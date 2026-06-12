""" Plotting utilities for the analysis of the results and to generate easily plots 
from GalaPy structures and functions.
"""
######################################################################################
# External imports
import os
import warnings
import fnmatch
import difflib
from collections.abc import MutableMapping as MM

# As the sane of mind do
import numpy

# For plotting
import matplotlib.pyplot as plt
from matplotlib.style import available as mpl_sty_av
from matplotlib import get_backend
_mpl_backend = get_backend()

# Matplotlib setup of the GalaPy default layout
plt.rcParams['figure.figsize'] = (7,5)
plt.rcParams['font.size'] = 14
plt.rcParams['lines.linewidth'] = 2.
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"

# this should select seaborn-v0_8-deep
# in all matplotlib installations
# more recent than version 3.6
plt.style.use(difflib.get_close_matches('seaborn-v0_8-deep', mpl_sty_av)[0]) 
clr = plt.rcParams['axes.prop_cycle'].by_key()['color']

######################################################################################
# Internal imports

import galapy as gp
from galapy.internal.utils import filter_strings, shorten_string, quantile_weighted
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

_specs_default = {
    'sed_ax_kw' : {
        'xscale' : 'log', 
        'yscale' : 'log',
        'xlim'   : (1.e+0,1.e+10),
        'ylim'   : (1.e-10, 1.e+10),
        'xlabel' : '$\\lambda\\ [\\,\\AA\\,]$',
        'ylabel' : '$S_\\nu\\ [\\mathrm{mJy}]$',
    },
    'sed_legend_kw' : {
        'l1' : { 'loc' : 'upper right',
                 'fontsize' : 12 },
        'l2' : { 'loc' : 'upper left',
                 'ncol' : 2,
                 'fontsize' : 12 },
    },
    'pms_ax_kw' : {
        'xscale' : 'log',
        'yscale' : 'log',
        'ylim'   : (1.e-2,1.e+1),
        'ylabel' : 'Transmission'
    },
    'residuals_ax_kw' : {
        'xlim' : (1., 1.e+10),
        'ylim' : (-1.0, +1.0),
        'xscale' : 'log',
        'yscale' : 'linear',
        'ylabel' : 'residuals',
    },
    'textbox_kw' : {
        'loc' : 'upper right',
        'borderpad' : 0.2, 
        'frameon' : False, 
        'prop' : {
            'fontsize' : 12,
            'bbox' : {
                'boxstyle' : 'round',
                'facecolor' : 'white',
                'alpha' : 0.7
            },
        },
    },
    'getdist_settings' : { 
        'mult_bias_correction_order' : 1,
        'smooth_scale_2D' : 4, 
        'smooth_scale_1D' : 4, 
        'fine_bins': 128,
        'fine_bins_2D' : 128,
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

_derived_quantity_meta = {
    'Mstar' : { 'label' : r'M_\star',          'log' : True  },
    'Mdust' : { 'label' : r'M_\mathrm{dust}',  'log' : True  },
    'Mgas'  : { 'label' : r'M_\mathrm{gas}',   'log' : True  },
    'Zstar' : { 'label' : r'Z_\star',           'log' : False },
    'Zgas'  : { 'label' : r'Z_\mathrm{gas}',   'log' : False },
    'SFR'   : { 'label' : r'\psi',              'log' : True  },
    'TMC'   : { 'label' : r'T_\mathrm{MC}',    'log' : False },
    'TDD'   : { 'label' : r'T_\mathrm{DD}',    'log' : False },
}

def show_default_dict ( name = None ) :
    """ Shows the default dictionary for the kwargs dictionary with name 'name'.
    """
    if name is not None :
        return _specs_default[ name ]
    return _specs_default

######################################################################################
# Functions for plotting SEDs

def sed_layout ( redshift, frame, ax = None, **kwargs ) :
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
      (Optional) an instance of matplotlib axes
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
    default_ax_kw = dict( _specs_default['sed_ax_kw'] )
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

def _errorbar_with_uplims( xx, yy, ee, lo, ax = None, **kwargs ) :

    if ax is None :
        ax = plt.gca()
    kwargs = dict(kwargs)
    marker = 'o'
    if 'marker' in kwargs.keys() :
        marker = kwargs.pop('marker')
    
    Pdata = ax.errorbar(
        xx[~lo], yy[~lo], ee[~lo],
        label='data',
        marker=marker,
        markerfacecolor='white',
        markeredgewidth=2.0, 
        markersize=8.0,
        elinewidth=2.0, 
        linestyle='none',
        zorder=3,
        **kwargs
    )
    color = None
    if 'c' in kwargs.keys() :
        color = kwargs.pop( 'c' )
    elif 'color' in kwargs.keys() :
        color = kwargs.pop( 'color' )
    else :
        color = Pdata.lines[0].get_color()
    Puplims = ax.errorbar(
        xx[lo], 2*ee[lo], 1*ee[lo],
        uplims=lo[lo],
        label='uplims',
        capsize = 2.0, 
        capthick=2.0,
        elinewidth=2.0,
        marker='none', 
        linestyle='none',
        zorder=3,
        color = color,
        **kwargs
    )
    return Pdata, Puplims

def sed_obs ( ll, ff, ee, lo, redshift = None, frame = 'rest', ax = None, ax_kwargs = {}, return_uplim_label=False ) :

    ax = sed_layout(redshift, frame, ax, **ax_kwargs)

    Pdata, Puplims = _errorbar_with_uplims( ll, ff, ee, lo, ax=ax )
    
    return Pdata

def sed_components ( ll, components,
                     redshift = None, frame = 'rest', ax = None, ax_kwargs = {} ) :
    ax = sed_layout(redshift, frame, ax, **ax_kwargs)
    Pcomp = []
    for k,L in components.items() :
        Pcomp.append(ax.plot(ll, L, '-', lw = 1.5, label=k)[0])
    return Pcomp

def sed_flux ( ll, flux,
               color = 'black', label = 'total', 
               redshift = None, frame = 'rest', ax = None, ax_kwargs = {} ) :
    ax = sed_layout(redshift, frame, ax, **ax_kwargs)
    return ax.plot( ll, flux, ls='-', color=color, lw=1.5, label=label)[0]

def sed_1sigma2sigma ( ll, flux, err, center = False, color = 'gray', center_label = 'mean', 
                       redshift = None, frame = 'rest', ax = None, ax_kwargs = {} ) :
    ax = sed_layout( redshift, frame, ax, **ax_kwargs )
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

def sed_flux_res ( res,
                   model = None,
                   observation = None,
                   plot_observation = False,
                   plot_components = False, 
                   plot_contours = True,
                   frame = 'both',
                   show_legend = True,
                   legend_kwargs = {},
                   ax = None, 
                   ax_kwargs = {}) :
    """ Plots the formatted SED flux from a galapy.sampling.Results instance.
    Optional specs that can be activated or de-activated are the observed points,
    1- and 2-sigma contours around the mean of the sampled SEDs and the different
    physical components contributing to the total best-fit SED.
    
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
        to plot the wavelength-axis (default is ``'both'``)
    show_legend : bool
        (Optional)
    legend_kwargs : dict
        (Optional)
    ax : matplotlib.axes.Axes
        (Optional) an instance of matplotlib axes
    ax_kwargs : dict
        (Optional) keyword arguments to pass to the function ``sed_layout``
    
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

    # Set the wavelength (lambda) axis:
    # If 'redshift' is a free parameter the local variable
    # will be set to the best-fit value, otherwise the
    # spectroscopic redshift will be extracted from the model.
    redshift = None
    lx = None
    obs_scale = frame in set(['obs', 'both'])
    if 'galaxy.redshift' in res.get_sampling_params().par_free :
        redshift = res.get_bestfit( 'params' )['redshift']
        lx = model.wl()
        if obs_scale :
            lx *= ( 1 + redshift )
    else :
        redshift = model.redshift
        lx = model.wl( obs = obs_scale )
        
    # Set the layout
    ax = sed_layout( redshift, frame, ax = ax, **ax_kwargs )
    
    ##################################################################
    # Here the plots start
    
    legend_primary   = []
    legend_secondary = []
    if show_legend :
        legend_kw = dict( _specs_default['sed_legend_kw'] )
        legend_kw.update( legend_kwargs )
        
    
    # Plot the observational dataset
    if plot_observation :
        ll, ff, ee, lo = ( observation.pms.lpiv,
                           observation.fluxes,
                           observation.errors,
                           observation.uplims )
        legend_primary.append(
            sed_obs( ll, ff, ee, lo,
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
        legend_secondary += sed_components( lx, components,
                                            redshift = redshift,
                                            frame = frame,
                                            ax=ax )
        
    if flux_bestf is None :
        flux_bestf = res.get_bestfit( 'SED' )
    legend_primary.append( sed_flux( lx, flux_bestf, label = 'best-fit',
                                     redshift = redshift,
                                     frame = frame,
                                     ax=ax ) )
        
    if plot_contours :
        flux_mean = res.get_mean( 'SED' )
        flux_err  = res.get_std( 'SED' )
        legend_primary += sed_1sigma2sigma( lx, flux_mean, flux_err,
                                            redshift = redshift,
                                            frame = frame,
                                            ax=ax )
    
    ##################################################################
    # Finalise and return
    
    # LEGEND:
    if show_legend :
        lpri = ax.legend(handles = legend_primary, **legend_kw['l1'] )
        ax.add_artist(lpri)
        if len(legend_secondary)>0 :
            lsec = ax.legend(handles = legend_secondary, **legend_kw['l2'] )
            ax.add_artist(lsec)
        
    return(ax)


def sed_residuals_res ( res,
                        frame = 'both', 
                        plot_contours = False, 
                        plot_chi2 = True,
                        ax = None,
                        text_kwargs = {},
                        ax_kwargs = {}) :
    """
    Plots the formatted residuals with respect to the best-fit model from a sampling run.
    
    Parameters
    ----------
    res : Results instance
        A ``Results`` instance from a sampling run.
    frame : string
        One among {``'rest'``, ``'obs'``, ``'both'``}, choose the frame, observed- or rest-frame,
        to plot the wavelength-axis (default is ``'both'``)
    plot_contours : bool
        Whether to plot the 1- and 2-sigma contours around the mean SED (dafault is ``False``).
    plot_chi2 : bool
        Whether to plot a text box with the best-fit value of the reduced 
        chi2 (dafault is ``True``).
    ax : matplotlib.axes.Axes
        (Optional) an instance of matplotlib axes
    text_kwargs : dict
        (Optional) keyword arguments to pass to the matplotlib.offsetbox.AnchoredText class
        (regulates shape and text-formatting of the text-box for the eventual chi2 plot)
    ax_kwargs : dict
        (Optional) keyword arguments to pass to the function ``sed_layout``
    
    Returns
    -------
    : matplotlib.axes.Axes
        An instance of ``matplotlib.axes.Axes`` with the plotted fluxes
    """
    
    # Get sampled model
    pbest = res.get_bestfit('params')
    pgxy = res.get_model()
    pgxy.set_parameters(**pbest)
    redshift = pgxy.redshift

    # Get reference observation
    obs = res.get_observation()
    ll, ff, ee, lo = ( obs.pms.lpiv,
                       obs.fluxes,
                       obs.errors,
                       obs.uplims )
    lsed = pgxy.wl(obs=(frame != 'rest'))

    # compute residuals
    chi = res.get_residuals()
    
    # Set axes
    ax_kw = dict(_specs_default['residuals_ax_kw'])
    ax_kw.update(ax_kwargs)
    ax = sed_layout(redshift, frame=frame, ax = ax, **ax_kw )
    
    #
    if plot_contours :
        sedstd = res.get_std('SED')
        wnstd0 = sedstd > 0.0
        _ = sed_1sigma2sigma( lsed[wnstd0], 
                              (res.get_mean('SED')-res.get_bestfit('SED'))[wnstd0]/sedstd[wnstd0], 
                              numpy.ones_like(lsed[wnstd0]),
                              center=True, redshift=redshift, frame=frame, ax = ax)
        
    _ = sed_flux( lsed, numpy.zeros_like(lsed), redshift=redshift, frame=frame, ax = ax )
    Pdata, *_ = _errorbar_with_uplims(
        ll[~lo], chi[~lo],
        numpy.zeros_like(ee[~lo]),
        numpy.zeros_like(lo[~lo]),
        ax = ax,
    )
    obs_color = Pdata.lines[0].get_color()
    _ = _errorbar_with_uplims(
        ll[lo],
        numpy.zeros_like(chi[lo]),
        numpy.zeros_like(ee[lo]),
        numpy.zeros_like(lo[lo]),
        ax = ax,
        color = obs_color,
        marker = 'x'
    )
    
    # Plot eventual text-box
    if plot_chi2 :
        from matplotlib.offsetbox import AnchoredText
        text_kw = dict(_specs_default['textbox_kw'])
        text_kw.update(text_kwargs)
        
        textstr = f'$\\chi^2_\\mathrm{{red}} = {numpy.sum(chi**2)/(len(ll)-res.ndim):.2f}$'
        text = AnchoredText(textstr, **text_kw)
        _ = ax.add_artist(text)
        
    return ax

######################################################################################
# Functions for plotting the photometric system

def photometric_system ( obj, colors = None, ax = None, ax_kwargs = {} ) :
    """ Generates a visualisation of the photometric system bandpass transmissions.
    
    Parameters
    ----------
    obj : object
        An instance of ``galapy.PhotometricSystem.PMS`` or any other object with a 
        ``pms`` attribute which is an instance of ``galapy.PhotometricSystem.PMS``
        (e.g. ``galapy.Galaxy.GXY`` or ``galapy.Observation.OBS``)
    colors : iterable
        (Optional) a list or iterable of valid colors. One color per bandpass is necessary
    ax : matplotlib.axes.Axes
        (Optional) an instance of matplotlib axes
    ax_kwargs : dictionary
        (Optional) the eventual keyword arguments to pass to the ``set`` function
        of ``ax``.

    Returns
    -------
    : matplotlib.axes.Axes
    an instance of matplotlib axes
    """
    from galapy.PhotometricSystem import PMS
    
    # Check input object (it should be or contain a PMS instance)
    if not isinstance(obj, PMS) :
        if hasattr(obj, 'pms') and isinstance( obj.pms, PMS ) :
            pms = obj.pms
        else :
            raise AttributeError( 'Input object ``obj`` must be an instance of '
                                  'either ``galapy.PhotometricSystem.PMS``, '
                                  'or any other object with a ``pms`` attribute '
                                  'which an instance of ``galapy.PhotometricSystem.PMS`` '
                                  '(e.g. ``galapy.Galaxy.GXY`` or ``galapy.Observation.OBS``)' )
    else :
        pms = obj
    
    # generate or check color iterable
    if colors is None :
        colors = plt.cm.plasma(numpy.linspace(0.,1.,len(pms.keys)))
    elif len(colors) != len(pms) :
        raise AttributeError( 'Attribute ``color`` must be a list or iterable of valid colors '
                              'with the same length of the photometric system' )
            
    # If not passed, attach to the latest axis or generate new one
    if ax is None :
        ax = plt.gca()
        
    # Set axes properties
    default_ax_kw = dict( _specs_default['pms_ax_kw'] )
    default_ax_kw.update(ax_kwargs)
    _ = ax.set(**default_ax_kw)
    xax2 = ax.secondary_xaxis('top', functions=(lambda x : x, lambda x : x))
    _ = xax2.set( xscale = default_ax_kw['xscale'], xlabel = 'wavelength $[\\AA]$' )

    # Set background color
    ax.set_facecolor('white')
    
    # Plot
    xticks = []
    xticklabels = []
    for k, c in zip( pms.keys, colors ) :
        v = pms.bpt[k]
        xticks += [v.get_lpiv()]
        xticklabels += [shorten_string(k)]
        ax.plot( v.get_xaxis(), numpy.asarray( v.get_xaxis() )*numpy.asarray( v.get_yaxis() ), color = c)
        ax.axvline( v.get_lpiv(), color = c, ls = '--', lw = 2.)
    _ = ax.set_xticks(xticks)
    _ = ax.set_xticklabels(xticklabels, rotation=90, fontsize = 11)

    return ax

######################################################################################
# Functions for plotting posteriors

def corner_res ( res, handler = None, which_params = None, getdist_settings = None,
                 param_limits = 'auto', plot_titles = True, mark = 'bestfit',
                 titles_kw = {}, triangle_kw = {}, marker_kw = {} ) :
    """ Plot the triangle plot with the 2D posteriors of a sampling run 
    and 1D marginals on the diagonal.
    
    Parameters
    ----------
    res :  Results instance
        A ``Results`` instance from a sampling run.
    handler : ModelParameters instance
        (Optional) The ``ModelParameters`` corresponding to the given sampling run
    which_params : str or sequence of str
        Either a single string or a sequence of strings. Name of the parameters to show on the 
        triangle plot. Also accepts wildcards (e.g. ``which_params = 'sfh*'`` will show all the
        parameters that contain the sub-string ``sfh``).
    getdist_settings : dict
        a dictionary of custom analysis settings to pass to the ``getdist.MCSamples`` 
        settings argument 
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
          sampling run and the values are sequences of length = 2 with the limits
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
    # necessary to fix the hard-coded matplotlib.use('Agg') forced by getdist
    plt.switch_backend(_mpl_backend)

    ############################################################################
    # Set default kwarg-dictionaries
    
    default_getdist_settings = dict( _specs_default['getdist_settings'] )
    default_titles_kw        = dict( _specs_default['titles_kw'] )
    default_marker_kw        = dict( _specs_default['marker_kw'] )
    default_triangle_kw      = dict( _specs_default['triangle_kw'] )
    
    ############################################################################
    # Check inputs
    
    if not isinstance( res, Results ) :
        raise AttributeError( 
            'Attribute "res" should be an instance of type ``Results``'
        )
    
    if handler is None :
        handler = res.get_sampling_params()
    else :
        if not isinstance( handler, gp.Handlers.ModelParameters ) :
            raise AttributeError(
                'Attribute "handler" should be an instance of type ``ModelParameters``'
            )
    
    if which_params is None :
        which_params = handler.par_free
    else :
        try :
            which_params = filter_strings( handler.par_free, which_params )
        except ValueError :
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


def corner_derived ( res, which_keys = None, log_scale = None,
                     getdist_settings = None,
                     param_limits = 'auto', plot_titles = True, mark = 'bestfit',
                     titles_kw = {}, triangle_kw = {}, marker_kw = {} ) :
    """ Triangle plot of derived physical quantities from a sampling run.

    Mirrors ``corner_res`` but operates on the derived quantities stored in the
    ``Results`` instance (``Mstar``, ``Mdust``, ``Mgas``, ``Zstar``, ``Zgas``,
    ``SFR``, ``TMC``, ``TDD``) rather than the free sampling parameters.
    Samples for which any selected quantity is non-finite are silently excluded.

    Parameters
    ----------
    res : Results instance
        A ``Results`` instance from a sampling run.
    which_keys : sequence of str, optional
        Names of the derived quantities to include.  Defaults to all quantities
        in ``_derived_quantity_meta`` that are present in ``res``.
    log_scale : sequence of str, optional
        Keys to plot on a log10 axis.  Defaults to the per-quantity setting in
        ``_derived_quantity_meta`` (masses and SFR are log by default).
    getdist_settings : dict, optional
    param_limits : str, sequence or dict, optional
        'auto' (default), a list of (lo, hi) pairs, or a dict keyed by quantity name.
    plot_titles : bool, optional
        Show median and 68%% credible interval above each diagonal panel (default True).
    mark : str, optional
        'bestfit' (default), 'mean', or 'median'.
    titles_kw, triangle_kw, marker_kw : dict, optional

    Returns
    -------
    fig : matplotlib.pyplot.Figure
    axes : numpy.ndarray of matplotlib.axes.Axes
    """
    from getdist import plots, MCSamples
    plt.switch_backend( _mpl_backend )

    default_getdist_settings = dict( _specs_default['getdist_settings'] )
    default_titles_kw        = dict( _specs_default['titles_kw'] )
    default_marker_kw        = dict( _specs_default['marker_kw'] )
    default_triangle_kw      = dict( _specs_default['triangle_kw'] )

    if not isinstance( res, Results ) :
        raise AttributeError(
            'Attribute "res" should be an instance of type ``Results``'
        )

    # select keys present in the Results instance
    available = [ k for k in _derived_quantity_meta if hasattr( res, k ) ]
    if which_keys is None :
        which_keys = available
    else :
        which_keys = [ k for k in which_keys if k in available ]
        if len( which_keys ) == 0 :
            raise ValueError(
                'None of the requested keys are available in this Results instance.'
            )

    # log-scale flags
    if log_scale is None :
        use_log = { k : _derived_quantity_meta[k]['log'] for k in which_keys }
    else :
        use_log = { k : ( k in log_scale ) for k in which_keys }

    # mask out non-finite sentinel samples and zeros for log-scaled quantities
    good = numpy.ones( res.size, dtype = bool )
    for k in which_keys :
        vals = getattr( res, k )
        good &= numpy.isfinite( vals )
        if use_log[k] :
            good &= ( vals > 0 )
    if not good.any() :
        raise RuntimeError(
            'No finite samples found for the selected derived quantities.'
        )

    weights = res.weights[good]
    logl    = res.logl[good]

    # samples matrix — log10-transform where requested
    mat = numpy.column_stack( [
        numpy.log10( getattr( res, k )[good] ) if use_log[k]
        else getattr( res, k )[good]
        for k in which_keys
    ] )

    # axis labels
    labels = []
    for k in which_keys :
        lab = _derived_quantity_meta[k]['label']
        if use_log[k] :
            lab = r'\log_{10}\!\left(' + lab + r'\right)'
        labels.append( lab )

    # marker position
    markers = {}
    if mark is not None :
        if mark == 'bestfit' :
            idx = logl.argmax()
            markers = { k : float( mat[idx, i] ) for i, k in enumerate( which_keys ) }
        elif mark == 'mean' :
            markers = {
                k : float( numpy.average( mat[:, i], weights = weights ) )
                for i, k in enumerate( which_keys )
            }
        elif mark == 'median' :
            markers = {
                k : float( quantile_weighted( mat[:, i], 0.5, weights = weights ) )
                for i, k in enumerate( which_keys )
            }
        else :
            warnings.warn( "Argument ``mark`` should be one among 'bestfit', 'median', 'mean'." )

    # update kw dicts
    if getdist_settings is not None :
        default_getdist_settings.update( getdist_settings )
    default_triangle_kw.update( triangle_kw )
    default_titles_kw.update( titles_kw )
    default_marker_kw.update( marker_kw )

    # display limits
    if isinstance( param_limits, (list, tuple) ) :
        display_limits = dict( zip( which_keys, param_limits ) )
    elif isinstance( param_limits, MM ) :
        display_limits = dict( param_limits )
    else :
        display_limits = {}   # 'auto' or unrecognised → let getdist decide

    # sampler type hint for getdist
    sampler = None
    if res.sampler == 'dynesty' :
        sampler = 'nested'
    elif res.sampler == 'emcee' :
        sampler = 'mcmc'

    mcsamples = MCSamples(
        samples  = mat,
        loglikes = -logl,
        weights  = weights,
        names    = which_keys,
        labels   = labels,
        sampler  = sampler,
        settings = default_getdist_settings,
    )

    g = plots.get_subplot_plotter()
    g.settings.axis_tick_max_labels  = 3
    g.settings.num_plot_contours     = len( default_getdist_settings['contours'] )
    g.triangle_plot(
        [mcsamples], params = which_keys,
        filled       = True, title_limit = 0,
        param_limits = display_limits,
        markers      = markers,
        marker_args  = default_marker_kw,
        **default_triangle_kw
    )

    fig, axes = g.fig, g.subplots

    # replicate marker styling on diagonal panels
    if len( markers ) > 0 :
        for i in range( len( axes ) ) :
            axes[i, i].get_lines()[1].set( **default_marker_kw )

    if plot_titles :
        digits = default_titles_kw.get( 'digits', 2 )
        for i, ( k, lab ) in enumerate( zip( which_keys, labels ) ) :
            q16, q50, q84 = quantile_weighted(
                mat[:, i], (0.16, 0.5, 0.84), weights = weights
            )
            lo = q50 - q16
            hi = q84 - q50
            title = ( f'$ {lab} = '
                      f'{q50:.{digits}f}^{{+{hi:.{digits}f}}}_{{-{lo:.{digits}f}}} $' )
            axes[i, i].set_title( title, fontsize = g.settings.axes_fontsize )

    format_axes_ticks( fig, labelsize = g.settings.axes_fontsize )

    return fig, axes

######################################################################################
