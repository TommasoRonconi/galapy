.. _param_file:

The parameter file
==================

The parameter file can be generated in whatever position in the filesystem by calling from a terminal

.. code-block:: console

    $ galapy-genparams [--name/-n NAME]

Where the argument ``NAME`` can include any path in the filesystem.
Since the generated parameter file has python syntax, we automatically append a ``.py`` extension to the custom file name.
If no argument is passed the parameter file will be generated in the current directory:

.. code-block:: console

    $ ls
    galapy_hyper_parameters.py

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
In particular, the sequences of values can be a ``numpy array``, a python ``list``, a ``set`` or a ``tuple``, it is nevertheless necessary that they all have all the same dimensions.

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
	
These parameter will be used to build an object of type ``galapy.PhotometricSystem.PMS``.

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

   A thorough description 

4. ``do_Xray`` (default ``= False``)

5. ``do_Radio`` (default ``= False``)

6. ``do_AGN`` (default ``= False``)

7. ``lstep`` (default ``= None``)

8. ``noise_model`` (default ``= None``)

9. ``noise_kwargs`` (default ``= {}``)


.. _fixed_and_free_parameters:

Choose the fixed and free parameters
....................................

.. _sampling_and_output:

Sampling and output format choices
..................................
