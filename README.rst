Welcome to GalaPy!
==================

GalaPy is an extensible API for modelling Galactic emission.
It provides an easy-to-use Python user interface while the number-crunching is done with compiled, high-performance, object-oriented C++.

It is intended for researchers in the fields of Galaxy Formation and Evolution, Observational Astrophysics and Cosmology but also for theoretical and computational researchers interested in a lightning fast Galactic spectra simulator.
Galapy enables instantiating multi-component parameterized galaxy objects with a high level of customization.
It produces Spectral Energy Distributions (SED) in a matter of milliseconds on a single core with a minimal memory consumption. 

Even though the library is optimized for fitting photometric datasets, its simulation capabilities are way more flexible than this, as it allows for the extraction of many other physical properties of galaxies, such as attenuation curves, matter content evolution histories (divided by component), metallicity evolution, star formation histories and more.

* **Free software**: GPLv3 license;
* **GitHub repository**: https://github.com/TommasoRonconi/galapy
* **Python versions**: >=3.6

The only external dependency required is NumPy which is summoned both from the Python interface and from the C-Wrapping layer.

The build system tool used is Setuptools, which is also not part of the Python standard contrary to its ancestor, Distutils, but still ...

Lift-off TL;DR
--------------

* Info about SED, GFE [TBA]
* The preferred method to install the package is through :code:`pip` as it will install the most recent stable release:
  
  .. code-block:: console
     
     $ pip install galapy

  for further details, please refer to the `installation guide`_.

* Build a spectroscopic galaxy object, obtain SED

  .. code-block:: python

     import galapy as gp
     gxy = gp.Galaxy.GXY( age, redshift )
     sed = gxy.SED()

  Build a photometric galaxy object, obtain photo-SED
  
  .. code-block:: python

     import galapy as gp
     pgxy = gp.Galaxy.PhotoGXY( age, redshift )
     pgxy.build_photometric_system( 'filter1', 'filter2', 'filter3', ... )
     psed = pgxy.photoSED()

  Link to `tutorials`_ here

* Sampling the parameter space can be done from the command line in a terminal.

  **[Provided that the database has been correctly installed]**

  The steps required for running the sampling are just two:
  - first we will have to generate a parameter file, this can be done by running
    the utility command

    .. code-block:: bash

       $ galapy-genparams [--name/-n NAME]

    The generated file contains should be self-explanatory, this has to be
    modified according to the fitting the user has to perform.
  - Once the parameter file has been generated and properly modified, we can run

    .. code-block:: bash

       $ galapy-fit parameter_file.py [--serial/-s | --multiprocessing/-mp NCPU]

    which will run the sampling and authomatically store the results, as specified
    by the user in the parameter file.
    NOTE THAT the two optional arguments regulate whether to run the sampling
    serially or using shared-memory parallelism.
    The default behaviour is to run parallely on all the available CPUs.
     
* Link to `API docs`_ here 

.. _installation guide: ...
.. _tutorials: ...
.. _API docs: ...

