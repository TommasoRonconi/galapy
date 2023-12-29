.. _param_file:

Build and modify the parameter file
===================================

The parameter file can be generated in whatever position in the filesystem by calling from a terminal

.. code-block:: console

    $ galapy-genparams [--name/-n NAME|--SFH_model/-sfh SFH_MODEL]

Where

* the argument ``NAME`` can include any path in the filesystem (default is ``galapy_hyper_parameters``).
  Since the generated parameter file has python syntax, we automatically append a ``.py`` extension to the custom file name.
  If no argument is passed the parameter file will be generated in the current directory.

* the argument ``SFH_MODEL`` allows the user to choose a SFH model, if no choice is made, the value should be set by modifying the content of the generated parameter file.

* the ``galapy-genparams`` command can be called with the argument ``--help/-h`` to show an help message.
  
The generated file should be self-explanatory and has to be modified according to the fit the user has to perform.
In what follows we provide a thorough explanation of all the entries that are present in the generated file.

.. _import_obs_data:

Import the observational data
..............................

This first block of parameters defines the photometric observation that we are going to fit the model to.

First, the observational dataset has to be defined, this is done by setting the following python iterables of values:
   
.. code-block:: python
   
   bands = None
   fluxes = None
   errors = None
   uplims = None

The ``None`` values in ``bands``, ``fluxes`` and ``errors`` are placeholders that should be modified by the user.
In particular, the sequences of values can be a ``numpy array``, a python ``list`` or a ``tuple``, it is nevertheless necessary that they all have all the same dimensions.

1. ``bands``: a sequence of strings listing the bandpass transmission names associated to the observation, e.g.

   .. code-block:: python

      bands = ( 'GOODS.i', 'GOODS.v', 'GOODS.z', 'GOODS.b' )

2. ``fluxes``: the observed flux within the bandpass filters listed in the hyperparameter ``bands`` (floating point values)
3. ``errors``: errors associated to the transmissions listed in the hyperparameter ``fluxes`` (floating point values)
4. ``uplims``: a list of booleans (``True`` or ``False``) identifying which of the fluxes listed in hyperparameter ``fluxes`` are non-detections.
   (note that, in case all the measurements are detections, this value can be left to its default value, i.e. ``None``)

Then, we define the photometric system through the two hyperparameters

.. code-block:: python
   
   filters = list(bands)
   filters_custom = None

The default value of the ``filters`` assumes that all the transmissions listed in hyperparameter ``bands`` are present in the database.
If this is true it can be left as is.

The second parameter is optional and can be left to ``None`` if all the filters are present in the database.
In case this is not true there are 2 possibilities: either none of the filters is present in the database, or just a subset is present.
In the first case set ``filters = ()`` and follow the instructions in point 2 of the following enumerated list.
In the second case you will have to set both hyperparameters ``filters`` and ``filters_custom`` manually:

1. List the transmission filters already available in the database.
   In this case the ``filters`` hyperparameter should be set to a list or a tuple listing the sub-set of filters already present:

   .. code-block:: python

	filters = ( 'GOODS.i', 'GOODS.v', 'GOODS.z', 'GOODS.b' )
  
2. List the custom set of filters through the ``filters_custom`` hyperparameter.
   It has to receive a nested dictionary properly formatted. This means each element in the dictionary should have a
   key that defines the chosen name for the filter and an associated value that should be
   itself a dictionary with two keys: ``'wavelengths'`` and ``'photons'``;
   the first key should provide an array-like listing a set of wavelengths while the
   second the transmissions in units of photons corresponding to each of the wavelengths
   in the first array:
   
   .. code-block:: python

	filters = {
	    'custom1' : {
	        'wavelengths' : [ 999., 1000., 1001., 1499., 1500., 1501. ],
		'photons' : [ 0., 0., 1., 1., 0., 0. ]
	    },
	    'custom2' : {
	        'wavelengths' : [ 1999., 2000., 2001., 2499., 2500., 2501. ],
		'photons' : [ 0., 0., 1., 1., 0., 0. ]
	    }	    
	}

   in the example above the photometric system will contain two filters, ``'custom1'`` and
   ``'custom2'``, the first will be a top-hat function in the interval
   :math:`\lambda \in (1000, 1500)\,\mathring{A}`, while the second a top-hat function in the
   interval :math:`\lambda \in (2000, 2500)\,\mathring{A}`
	
These parameters will be used to build an object of type ``galapy.PhotometricSystem.PMS``.

.. tip::
   
   To check what are the filters available in the galapy database:

   .. code-block:: python

	from galapy.PhotometricSystem import print_filters
	print_filters()

   Note that the function also accepts arguments, for filtering by experiment
   (e.g. ``print_filters('Herschel')`` will only print on screen filters used
   in the Herschel experiment)

The last hyperparameter to set in this block defines how to treat eventual upper-limits present in the dataset

.. code-block:: python
   
   method_uplims = 'chi2'

If all the measurements can be treated as detections, this parameter will be ignored.
Otherwise the user can choose between 3 possible behaviours, which will translate in an additional term to the gaussian likelihood used to sample the parameter space.
The three possible methods are:

1. ``'chi2'`` (default): non-detections are treated exactly as detections with a large error;

2. ``'simple'``: a simple step-wise function setting the log-likelihood to :math:`-\infty` (i.e., zero probability) when the model predicts a flux larger than observed and
   to zero (i.e., probability equal to one) when the predicted flux is lower than the limit:

   .. math::

      f\left[\overline{S}_j, \overline{S}_j(\theta),\sigma_j\right] = \left\{
      \begin{aligned}
      &\ \text{-}\infty & \overline{S}_j(\theta) > \overline{S}_j\\
      &\\
      &\ 0& \text{otherwise}
      \end{aligned}
      \right.

3. ``'S12'``: Sawicki (2012) proposes a modification of the :math:`\chi^2` that consists of the integral of the probability of some observation up to the given proposed model.
   If the errors on data are Gaussian, this integral provides the following analytical expression for the corresponding log-likelihood:

   .. math::

      f\left[\overline{S}_j, \overline{S}_j(\theta),\sigma_j\right] =
      -2 \ln \left\{\sqrt{\frac{\pi}{2}} \sigma_j \left[1 + \text{erf}\left(\frac{\overline{S}_j - \overline{S}_j(\theta)}{\sqrt{2}\sigma_j}\right)\right]\right\}~.

   Even though it can be argued that using the expression above is the most formally correct way of accounting for upper limits when errors are Gaussian,
   the combination of logarithm and error function is particularly risky in computational terms.
   Specifically, it tends to hit the numerical limit of floating point numbers representation accuracy really fast, leading to undefined behaviour.

.. _define_model:

Define the physics of the galaxy model
......................................

In this block of hyperparameters the user can define the modelling strategies to adopt, both in terms of galaxy model and noise model.
All the parameters in this section can be mantained to their default value, this will not affect the possibility to start sampling.
Nonetheless, the user should set these values according to the needs of the examined source, as the scientific result of the sampling depends on this set of parameters.

1. ``cosmo`` (default ``= 'Planck18'``): this parameter defines what cosmological model to adopt in computing distances and ages.
   It is used to convert the emitted energy into flux and to check the age of the source against the age of the Universe at the given redshift.
   There are several pre-computed models, including ``'Planck15'``, ``'Planck18'``, ``'WMAP7'``, ``'WMAP9'``.
   Nonetheless it is possible to use user-defined models by passing a dictionary of pre-computed values of luminosity distance and age as a function of redshift for the chosen cosmology:

   .. code-block:: python
      
      cosmo = {
          'redshift' = [ ... ], # the redshift grid of chosen size N
	  'luminosity_distance' = [ ... ], # luminosity distance values defined on the N-size redshift grid
	  'age' = [ ... ] # age of the Universe values defined on the N-size redshift grid 
      }

   The values provided will be linearly interpolated within the provided redshift grid and linearly extrapolated from the extremes outside of the redshift grid.

2. ``sfh_model`` (default ``= 'insitu'``) the star-formation history model of choice. The available parameterised models are listed below:

   * ``'constant'``:

     .. math::

	\psi(t) = \psi_0
	
     where :math:`\psi_0` is a constant floating-point value expressed in units of :math:`M_\odot/\text{yr}` (parameter ``sfh.psi`` in the :ref:`tunable_params` table).

   * ``'delayedexp'``, a generalised version of the delayed exponential SFR:

     .. math::

	\psi(t)\propto \tau^{\kappa}\, \exp{(-\tau/\tau_\star)}
	
     where :math:`\tau_\star` is the characteristic star-formation timescale (parameter ``sfh.tau_star``) and :math:`\kappa` (parameter ``sfh.k_shape``) is a shape parameter for the early evolution;
     :math:`\kappa=0` corresponds to a pure exponential, while :math:`\kappa=1` to the standard delayed exponential,
     the model above is normalized by a factor defined in the free-parameter ``sfh.psi_norm``.

   * ``'lognormal'``:

     .. math::

	\psi(t)\propto \dfrac{1}{\tau}\, \dfrac{1}{\sqrt{2\pi\sigma_\star^2}}\, \exp\left[-\dfrac{\ln^2(\tau/\tau_\star)}{2\,\sigma_\star^2}\right]

     where :math:`\tau_\star` (parameter ``sfh.tau_star``) and :math:`\sigma_\star` (parameter ``sfh.sigma_star``) control the peak age and width.
     Also in this case, the above model is multiplied by a free-parameter ``sfh.psi_norm``.

   * ``'insitu'``, is our default, physically motivated model. It provides a SFR with shape:

     .. math::

	\psi(t)\propto e^{-x}-e^{-s\gamma\, x}
	
    where :math:`x\equiv\tau/s\,\tau_\star` with :math:`s \approx 3` a parameter related to gas condensation,
    while :math:`\gamma` is a parameter including gas dilution, recycling and the strength of stellar feedback.
    Its value is given by :math:`\gamma \equiv 1 + \mathcal{R} + \epsilon_\text{out}`, where :math:`\mathcal{R}`
    is the recycled gas fraction and :math:`\epsilon_\text{out}\approx 3[\psi_\text{max}/M_\odot\text{yr}^{-1}]^{-0.3}`
    is the mass loading factor of the outflows from stellar feedback.
    Therefore, the parameter :math:`\gamma` is completely determined in terms of the free parameter :math:`\psi_\text{max}` and,
    eventually, by the age of the galaxy :math:`\tau`.
    The free parameters in this model are ``sfh.psi_max``, entering the relation for :math:`\epsilon_\text{out}` and normalising the
    model, and ``sfh.tau_star``.

   * ``'interpolated'``: the user can provide a grid of values for the age and SFR at each time, the SFH will then be interpolated over said grid (not for sampling, see the API documentation).

3. ``ssp_lib`` (default ``= 'parsec22.NT'``) defines which Simple Stellar Population table to use. There are several available a complete list can be printed on screen
   by calling

   .. code-block:: python

      from galapy.CompositeStellarPopulation import print_ssp_libs
      print_ssp_libs()

   A thorough description of the differences between the several SSP libraries can be found in the presentation paper or in :ref:`/notebooks/choose_ssp_lib.ipynb`.

4. ``do_Xray`` (default ``= False``) by setting this boolean to ``True`` the SED will be extended down to wavelength :math:`\lambda = 1\,\mathring{A}`.
   The choice affects both the stellar continuum (by adding the X-Ray emission due to stars in binary systems) and the eventual AGN emission (if ``do_AGN = True``).

5. ``do_Radio`` (default ``= False``) this parameter controls the eventual addition of Radio emission due to the stellar continuum. In particular
   
   * Super-Nova Synchrotron
   * Nebular Free-Free emission
     
   The system authomatically accounts for the SSP library chosen, as libraries of the ``parsec22.NT`` family already include synchrotron and libraries of the ``parsec22.NTL``
   family include both synchrotron and nebular emission (free-free and lines, see :ref:`/notebooks/choose_ssp_lib.ipynb` for further details).

   .. note::

      If a library from the ``parsec22.NTL`` family is chosen, setting ``do_Radio = True`` does not have any effect as all the sources of stellar Radio continuum are already
      included on top of the SSPs. This also means that the radio contribute will be accounted for in the energy balance computation (see presentation paper for details).

6. ``do_AGN`` (default ``= False``), if this parameter is set to ``True``, templated emission from the AGN will be added to the final SED model.
   The templates we use are those produced by `Fritz et al., 2006 <https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2966.2006.09866.x>`_ .
   Further description on the templated models and of the possible combinations of parameters can be found in the documentation page of the :py:mod:`galapy.ActiveGalacticNucleus` module.
   Note that, while the choice of template is defined on a discrete distribution of parameters, its contribute to the overall emission is given in terms of some reference emission,
   through the free parameter ``fAGN``.
   In particular, by default we use the total InfraRed luminosity due to the non-AGN component of the emission as a reference.

7. ``lstep`` (default ``= None``) this parameter regulates the eventual sub-sampling to be applied to the wavelength grid, it accepts

   * an integer number: consider a wavelength grid point every ``lstep`` points;
   * a sequence of integers, listing the indices of the elements of the wavelength grid to take into account;
   * a boolean mask with the same size of the wavelength grid, where ``False`` values are masked elements.

   .. warning::

      Since the time required for computing models is tightly related to the resolution of the wavelength grid, this hyper-parameter allows (by sub-sampling the default grid) to
      speed up the computation. The wavelength grid is defined by the SSP library chosen, sub-sampling it means to have a lower resolution on all the quantities derived from the fit.
      In particular, it will directly result in an under-estimation of the energy balanced between attenuated and re-radiated flux. The safest choice is to avoid any sub-sampling.

8. ``noise_model`` (default ``= None``), choose a model for treatment of noise. Currently available models are:

   * ``'calibration_error'``: adds a systematic unknown source of error to the measurement depending on a single additional **nuissance** parameter ``f_sys``, modifying the measured uncertainties as
     :math:`\tilde{\sigma}_i^2(\theta, f_\text{sys}) = \sigma_i^2 + f_\text{sys}^2 \bar{S}_i^2(\theta)`, consistently, the log-likelihood takes on the form

     .. math::

	\ln\mathcal{L}(\bar{S}| \theta, f_\text{sys} ) \equiv - \dfrac{1}{2}\sum_i\biggl\{\dfrac{[\bar{S}_i - \bar{S}_i(\theta)]^2}{\tilde{\sigma}_i^2(\theta, f_\text{sys})} + \ln\bigl[2\pi \tilde{\sigma}_i^2(\theta, f_\text{sys})\bigr]\biggr\}

     where :math:`\bar{S}_i` is the i-th measured flux and :math:`\bar{S}_i(\theta)` the model prediction for the i-th measured flux
     (see presentation paper for further details, some information on this can also be found in the
     `emcee online documentation <https://emcee.readthedocs.io/en/stable/tutorials/line/#maximum-likelihood-estimation>`_).

   * ``None``: no-particular action is performed to model eventual systematics in the error measurements,
     it is assumed that all the sources of uncertainty have been already accounted for when estimating errors on the observed fluxes.
     (note that this might not be the safest choice, especially when etherogeneous datasets are being considered).

9. ``noise_kwargs`` (default ``= {}``) eventual additional arguments to be passed to the noise constructor class (see documentation of the :py:mod:`galapy.Noise` module.


.. _fixed_and_free_parameters:

Choose the fixed and free parameters
....................................

Define the chosen parameterisation of the galaxy and noise models in this block.
The user has to define the content of two python dictionaries: ``galaxy_parameters`` and ``noise_parameters``.
A complete list of the available tunable parameters is provided in :ref:`tunable_params`.
The logic used is the following:

* if a parameter does not appear in the dictionary as a ``key : value`` pair, it will be set to its default value (see last column of the table in :ref:`tunable_params`).
* if a parameter appears in the dictionary there are two possibilities:
  
  - fix the parameter to a value ``float_value`` by passing a ``key : float_value`` pair, where the ``float_value`` is some floating point number.
    This choice will not add a dimension to the sampled parameter space but will change the default value of the target parameter.
  - Set a free parameter to vary within some uninformative uniform prior whose limits are user defined.
    This is achieve by passing a tuple as value of the ``key : value`` pair. In particular the tuple will be something like ``( [lower_limit, upper_limit], boolean )``,
    where ``lower_limit`` and ``upper_limit`` are the lower and upper limits of the uniform uninformative prior, respectively, while the ``boolean`` states whether the
    given free-parameter has to be considered logarithmic or linear.

    + a **logarithmic** free-parameter, e.g. the age of the galaxy, will look something like: ``'age' : ( [6., 11.], True )``,
      meaning that samples will be drawn from the interval :math:`10^6\ \text{years} < \text{age} < 10^{11}\ \text{years}`.
    + a **linear** free-parameter, e.g. the redshift of the galaxy, will look something like: ``'redshift' : ( [0., 10.], False )``,
      meaning that samples will be drawn from the interval :math:`0 < z < 10`.

    Each free-parameter adds a dimension to the parameter space.

.. note::

   The default parameter file is built with all the available tunable parameters already inserted in the dictionary, all of them set as free parameters with a large prior.
   We also made a priori some of these parameters logarithmic, as sampling in the logarithmic space for some of these of parameters is most of the times (but not always!)
   recommended, when the prior volume span different orders of magnitude. Nevertheless we strongly recommend the user to tune the size of the prior and the number of free
   and fixed parameters to taylor the specific use-case.
   Furthermore, when reading the parameter file, the system will authomatically ignore the free parameters that have been included in the two dictionaries but are not used
   in the chose model (i.e. if ``do_AGN = False`` all the AGN-parameters will be ignored, even though included in the ``galaxy_parameters`` dictionary).
   Nevertheless, leaving them in the dictionary will trigger a warning that will be print on screen when the fitting command is called (i.e. ``galapy-fit parameter_file.py``).
   In order not to pollute the STDERR, we recommend to erase from the dictionary all the parameters that are not used in the chosen galaxy model.
   The reason why we have chosen this strategy over a silently ignoring the redundant parameters, is to raise attention in the user in eventual human-driven-mistakes in the
   choice of parameters.

   Another relevant case in which more than the necessary parameters are included by default when the parameter file is generated concerns the SFH model.
   Different SFH models provide different parameterisations and when the command ``galapy-genparams`` is called **without** the ``--SFH_model/-sfh SFH_MODEL`` option, all
   the possible parameterisations from all the available models are included in the ``galaxy_parameters`` dictionary. If the unnecessary parameters are not removed from the dictionary,
   all the parameters that are not used in the SFH model chosen will raise a Warning message.
   By calling, e.g.,

   .. code-block:: bash

	$ galapy-genparams -sfh delayedexp

   the parameter file will be generated with ``sfh_model = 'delayedexp'`` (see point 2 of :ref:`define_model`) and only the parameters relevant to the Delayed Exponential SFH model
   will be present in the ``galaxy_parameters`` dictionary.
    
.. _sampling_and_output:

Sampling and output format choices
..................................

This last section allows to choose a sampler, regulates its behaviour and the choices for the output format.
The parameters to set are the following:

1. ``sampler`` (default ``= 'emcee'``) choose between

   * `Emcee sampler <https://emcee.readthedocs.io>`_ provides an implementation of the Markov-Chain Monte Carlo (MCMC) technique.
     Specifically, it implements an ensemble sampler with affine invariance `Goodman & Weare, 2010 <https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract>`_ that,
     by instantiating many test particles (*walkers*) in the parameter space, builds first order Markov sequences of proposals that are tested against the likelihood.
     The dynamics of this system of particles is regulated by the requirement that, at each new step, a better estimate of the parameters is drawn.

   * `Dynesty sampler <https://dynesty.readthedocs.io>`_ implements Dynamic Nested Sampling `(Higson et al., 2019) <https://ui.adsabs.harvard.edu/abs/2019S%26C....29..891H/abstract>`_,
     a generalised version of nested sampling `(Skilling, 2004) <https://ui.adsabs.harvard.edu/abs/2004AIPC..735..395S/abstract>`_ where the number of test particles
     (here *live points*) is dynamically increased in regions of the posterior where a higher accuracy is required.
     The parameter space is modelled as a nested set of iso-likelihood regions that are sampled until the overall evidence reaches a stopping criterion set by the user.
     In our default hyper-parameters set-up we provide an 80%/20% posterior/evidence split and we model the posterior space with multiple ellipsoids,
     as we expect to have multiple peaks and correlations when sampling high dimensional parameter spaces.
     We use the default stopping function

     .. math::

	\mathcal{S}(f_p, s_p, s_{\mathcal{Z}}, n) \equiv f_p \times \frac{\mathcal{S}_p(n)}{s_p} + (1 - f_p) \times \frac{\mathcal{S}_\mathcal{Z}(n)}{s_{\mathcal{Z}}} < 1
	
     where :math:`f_p` is the fractional importance we place on posterior estimation (20%, as mentioned above), :math:`\mathcal{S}_p` is the posterior stopping function,
     :math:`\mathcal{S}_\mathcal{Z}` is the evidence stopping function, :math:`s_p` is the posterior "error threshold", :math:`s_\mathcal{Z}` is the evidence error threshold,
     and :math:`n` is the total number of Monte Carlo realisations, used to generate the posterior/evidence stopping values.

   When sampling high-dimensional large volumes the degeneracy between parameters can easily generate a complex posterior topology,
   such as multiple peaks on some parameters or non-linear correlations.
   Our suggestion for an optimal usage is to rely on dynamic nested sampling in this case.
   As an empirical rule of thumb, we can recommend to rely on nested sampling when the number of free parameters is larger than :math:`~5` and when it is not necessary
   to include extremely complex priors (as this, even though feasible, is not trivial in nested sampling).

   On the other hand, MCMC provides a more straightforward interface to the inclusion of sophisticated priors and proves to be efficient and to possibly converge faster when
   working on smaller and well-behaved volumes, i.e. when multiple peaks and complex correlations among parameters are not to be expected.

2. ``nwalkers`` (default ``= 64``) and ``nsamples`` (default ``= 4096``). These two hyperparameters are used if and only if the emcee sampler is chosen (see previous point).
   They regulate the number of walkers that will be used and the total number of samples that will be extracted per each walker.

   .. note::

      The total number of samples extracted with the default values will be ``nwalkers`` :math:`\times` ``nsamples`` :math:`= 2^{18}`, corresponding to an expected total execution
      time of about :math:`\approx 15\div20` minutes. Even though, in MCMC sampling with *emcee* using our default hyperparameters, the convergence is not guaranteed, this number
      of samples should be enough for most of well-behaved posteriors.
      More sophisticated convergence criteria to stop the sampler can be included by modifying the ``sampler_kw`` and ``sampling_kw`` hyperparameters.

   .. tip::

      The number of walkers (``nwalkers``) should **always** be :math:`\ge 4\times N_\text{dim}`, where :math:`N_\text{dim}` is the number of free-parameters. 

3. ``sampler_kw`` (default ``= {}``) additional keyword arguments to be passed to the constructor of the chosen sampler.

   * for **emcee**: see documentation of the `Ensamble Sampler <https://emcee.readthedocs.io/en/stable/user/sampler/>`_.
     The default values we have set internally (see :py:mod:`galapy.sampling.Sampler`) are

     .. code-block:: python

	emcee_default_sampler_kw = { 'moves' : None, 'args' : None, 'kwargs' : None }

   * for **dynesty**: see documentation of the `Dynamic Nested Sampler <https://dynesty.readthedocs.io/en/stable/api.html>`_.
     The default values we have set internally (see :py:mod:`galapy.sampling.Sampler`) are

     .. code-block:: python

	dynesty_default_sampler_kw = { 'bound' : 'multi', 'sample' : 'rwalk',
		                       'update_interval' : 0.6, 'walks' : 25,
                                       'bootstrap' : 0 }

4. ``sampling_kw`` (default ``= {}``) additional keyword arguments to be passed to the function running the sampling with the chosen sampler.

   * for **emcee**: see documentation of the ``run_mcmc`` method of the `Ensamble Sampler <https://emcee.readthedocs.io/en/stable/user/sampler/>`_.
     The default values we have set internally (see :py:mod:`galapy.sampling.Sampler`) are

     .. code-block:: python

	emcee_default_sampling_kw = { 'log_prob0' : None, 'rstate0' : None, 'blobs0' : None,
                                      'tune' : False, 'skip_initial_state_check' : False,
                                      'thin_by' : 1, 'thin' : None, 'store' : True,
                                      'progress' : True, 'progress_kwargs' : None }
	
   * for **dynesty**: see documentation of the ``run_nested`` method of the `Dynamic Nested Sampler <https://dynesty.readthedocs.io/en/stable/api.html>`_.
     The default values we have set internally (see :py:mod:`galapy.sampling.Sampler`) are

     .. code-block:: python

	dynesty_default_sampling_kw = { 'nlive_init' : 1024, 'maxiter_init' : None,
                                        'maxcall_init' : None, 'dlogz_init' : 0.02,
                                        'nlive_batch' : 1024, 'maxiter_batch' : None,
                                        'maxcall_batch' : None,
                                        'maxiter' : sys.maxsize, 'maxcall' : sys.maxsize,
                                        'maxbatch' : sys.maxsize,
                                        'n_effective' : None, 'print_progress' : True,
                                        'stop_kwargs' : { 'target_n_effective' : 10000 } }

5. ``output_directory`` (default ``= ''``) choose the position in the filesystem where to store the results. The default empty string (``''``) means "store in current location",
   results will be saved in the directory from which the ``galapy-fit`` command has been called.

6. ``run_id`` (default ``= ''``) choose a common name to associate to all the output files. The empty string will trigger the usage of a *date-time* string with format
   ``'YYYYMMDDhhmm'`` where ``YYYY`` = four digits for the year (we are scheptical that our tool will still be used after year 9999 A.D.), ``MM`` two digits for the month, ``DD`` two
   digits for the day, ``hh`` two digits for the hour (in 24h-format), ``mm`` two digits for the minutes.

7. ``store_method`` (default ``= 'hdf5'``) two possible output formats are currently available in GalaPy:

   * ``'pickle'`` the standard Python serialisation protocol.  Results objects (see documentation of the :py:mod:`galapy.sampling.Results` module) are computed at the end of a sampling run
     then serialised and stored in non-volatile memory. The typical size of the output file can reach up to :math:`\sim 1` GB.
   * ``'hdf5'`` the Hierarchical Data Format (see, e.g., `Folk et al., 2010 <https://dl.acm.org/doi/10.1145/1966895.1966900>`_), a widespread method for storing heterogeneous data.

   .. tip::

      For almost all use-cases, the HDF5 format is a better choice than pickle as it is safer to distribute and backward/forward compatibility is guaranteed.
      Pickle should be used only for internal usage.

8. ``store_lightweight`` (default ``= False``) boolean available only if the HDF5 output format is chosen.
   
   * ``True``: store only samples coordinates, likelihood values and weights along with minimal additional information to re-build the models used in the sampling (typical size :math:`\sim 10` MB);
   * ``False``: along with the information available also stored with this option set as ``True`` option,
     all the additional derived quantities computed when building the ``Results`` object are stored
     (typical size up to :math:`\sim 1` GB, see documentation of the  :py:mod:`galapy.sampling.Results` module).

   .. note::

      The amount of information stored by choosing ``store_method = 'pickle'`` is equivalent to the combination ``store_method = 'hdf5'`` and ``store_lightweight = False``)

9. ``pickle_raw`` (default ``= False``) whether to pickle the sampler raw results, no analysis on the outputs is done, these data are not sufficient to reproduce the models used.
   (might be useful for analyzing the run statistics).

10. ``pickle_sampler`` (default ``= False``) whether to pickle the sampler at the end-of-run state. It is necessary to set this hyperparameter to ``True`` if the user wants to restart the run.

