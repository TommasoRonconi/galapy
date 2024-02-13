GalaPy: Spectral modelling tool for Galaxies in Python
======================================================

.. |test-badge| image:: https://github.com/TommasoRonconi/galapy/actions/workflows/tests.yml/badge.svg
   :alt: Testing Status

.. |build-badge| image:: https://github.com/TommasoRonconi/galapy/actions/workflows/build-wheels.yml/badge.svg
   :alt: Building

.. |docs-badge| image:: https://readthedocs.org/projects/galapy/badge/?version=latest
   :target: https://galapy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

|test-badge| |build-badge| |docs-badge|
	 
GalaPy is an open source, extensible API for modelling and fitting the Spectral Energy Distribution (SED) of galaxies from the X-ray to the radio band,
as well as the evolution of their components and dust attenuation/reradiation.
It provides functions, classes and terminal commands for Bayesian inference of galaxy properties from panchromatic photometric data,
as well as for the analysis and simulation of galaxies spectral properties.

.. image:: https://raw.githubusercontent.com/TommasoRonconi/galapy_database/main/images/GalaPy_Example.png
   :width: 100%
   :alt: If the image is not directly shown in the text, it can be found in the subdirectory `logo/GalaPy_Example.png`

GalaPy provides an easy-to-use Python user interface while the number-crunching is done with compiled, high-performance, object-oriented C++.

The library is currently optimized for fitting photometric datasets, nevertheless, its simulation capabilities are way more flexible than this,
as it allows for the extraction of many other physical properties of galaxies, such as attenuation curves, matter content evolution histories (divided by component),
metallicity evolution, star formation histories and more.

Galapy enables instantiating multi-component parameterized galaxy objects with a high level of customization.
It produces SEDs in a matter of milliseconds on a single core with a minimal memory consumption.
It has been developed with the aim of providing a fast SED simulator and, thus, to aid research in Galaxy Formation and Evolution,
both from the perspective of observational Astrophysics and Cosmology (thanks to its Bayesian statistical framework) as well as from the perspective of
theoretical and computational researchers interested in a modern modelling tool.

+-----------------------+-------------------------------------------+
| **Free software**     | GPLv3 license                             |
+-----------------------+-------------------------------------------+
| **GitHub repository** | https://github.com/TommasoRonconi/galapy  |
+-----------------------+-------------------------------------------+
| **Python versions**   | >=3.7                                     |
+-----------------------+-------------------------------------------+
| **Dependencies**      | ``setuptools``, ``numpy``, ``scipy``,     |
|                       | ``emcee``, ``dynesty``, ``matplotlib``,   |
|                       | ``getdist``, ``requests``                 |
+-----------------------+-------------------------------------------+

TL;DR
-----

Galaxies are extremely complex astrophysical objects resulting from the interaction of baryonic matter which has collapsed within a Dark Matter halo.
Their formation and evolution strongly depend on the interplay of several factors, including their matter reservoir and accretion history,
the environment they reside and the interactions with their neighbouring objects and, ultimately,
the large scale structure of the Universe and the physics regulating it on cosmological scales.
By studying the properties of individual galaxies, such as their luminosity, chemical composition, and star formation rate,
we can learn about how galaxies form and evolve over time as well as the cosmological conditions that lead to their assembly.

The Spectral Energy Distribution (SED) of a galaxy describes the distribution of its light across different wavelengths, from gamma rays to radio waves,
literally shedding light over the baryonic components and processes that contribute to the overall emission.
Modelling this emission is one of the primary tools of extra-galactic astronomy to constrain models of galaxy formation and evolution,
which are an essential part of our understanding of the Universe as a whole.

Install
.......

The preferred method to install the package is through :code:`pip` as it will install the most recent stable release:
  
.. code-block:: console
     
   $ pip install galapy-fit

for further details, please refer to the `installation guide`_.

Fitting through terminal commands
.................................

Sampling the parameter space can be done from the command line in a terminal.
The steps required for running the sampling are just two:
  
1. first we will have to generate a parameter file, this can be done by running
   the utility command

   .. code-block:: console

      $ galapy-genparams [--name/-n NAME | --SFH_model/-sfh MODEL_NAME ]

   The generated file should be self-explanatory and has to be
   modified according to the fit the user has to perform.
   A detailed guide to the generation and modification of the parameter file
   can be found in `param_file`_.
  
2. Once the parameter file has been generated and properly modified, we can run

   .. code-block:: console

      $ galapy-fit parameter_file.py [--serial/-s | --multiprocessing/-mp NCPU]

   which will run the sampling and authomatically store the results, as specified
   by the user in the parameter file.
   NOTE THAT the two optional arguments regulate whether to run the sampling
   serially or using shared-memory parallelism.
   The default behaviour is to run parallely on all the available CPUs.
   More details are provided in `photometric_fit`_.

.. note::
   GalaPy, in some of its components (e.g. SSP tables, PAH template), makes use of pre-computed functions that are
   available in the official database (`galapy_database`_). When one of the files in the database is accessed for the
   first time it will authomatically be downloaded into the user's filesystem
   (in the default location :code:`$HOME/.galapy/galapy_database`).
   This will of course require an internet connection and can partially slow down the computations.
   We therefore suggest, prior to first run, to download all the database by running

   .. code-block:: console

	$ galapy-download-database
   
   
Quick API hands-on
..................

The GalaPy API allows to directly access methods and classes modelling the different components
that contribute to the overall emission of a galaxy.
By the interplay of these components the final Spectral Energy Distribution (SED) emerges and
travels towards the observer.

In order to control the aforementioned interplay of components the module ``galapy.Galaxy`` implements classes of
type ``GXY``, from which the intrinsic luminosity and the flux at given distance can be retrieved.
An object of type ``GXY`` is built as follows

.. code-block:: python

   import galapy as gp
   gxy = gp.Galaxy.GXY( age = 1.e+9, redshift = 1.0 )

We have built a galaxy :math:`1 \text{Gyr}` old at redshift :math:`z = 1`.
We can always change the parameters of the galaxy we have built by calling the method
   
.. code-block:: python

   gxy.set_parameters( age = 5.e+9 )

For a complete list of the tunable parameters check the relative documentation page: `Free parameters`_.
To get the intrinsic emission from the galaxy and its flux as arriving at the observer we can call the
following two functions
   
.. code-block:: python

   # Intrinsic luminosity:
   L = gxy.get_emission()

   # Flux:
   F = gxy.SED()

Note that the function :code:`gxy.wl( obs = True/False )` returns the wavelength grid in the
observer's frame (:code:`obs = True`) and at rest frame (:code:`obs = False`). 

If, instead of the full spectrum, we want just the flux integrated within some transmission
bands, we will build a photometric galaxy object, and obtain the photo-SED
  
.. code-block:: python

   pgxy = gp.Galaxy.PhotoGXY( age = 5.e+9, redshift = 1.0 )
   pgxy.build_photometric_system( 'filter1', 'filter2', 'filter3', ... )
   pF = pgxy.photoSED()

Further details on the usage of functions and classes of the API are provided in the `tutorials`_
and in the `API documentation`_. 

.. _installation guide: https://galapy.readthedocs.io/en/latest/general/install_guide.html
.. _param_file: https://galapy.readthedocs.io/en/latest/guides/parameter_file.html
.. _photometric_fit: https://galapy.readthedocs.io/en/latest/guides/photometric_fit.html
.. _tutorials: https://galapy.readthedocs.io/en/latest/tutorials/physics_modules.html
.. _API documentation: https://galapy.readthedocs.io/en/latest/python_doc/api_toctree.html
.. _galapy_database: https://github.com/TommasoRonconi/galapy_database/
.. _Free parameters: https://galapy.readthedocs.io/en/latest/general/free_parameters.html

