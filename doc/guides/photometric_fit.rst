.. _photometric_fit:

Fit a model to photometric data
===============================

In broad band observations of astrophysical sources, the total flux is given by the convolution of the radiation coming from the source
with the response of the instrument.
GalaPy objects of type ``galapy.Galaxy.PhotoGXY`` can model this kind of measurements, when the transmission model of the different
observational bands is provided to the class.

Fitting through the command line
................................

Sampling the parameter space can be done from the command line in a terminal.
This will require to compile a templated parameter file that can be generated by calling

.. code-block:: console

    $ galapy-genparams [--name/-n NAME]

Where the argument ``NAME`` can include any path in the filesystem.
Since the generated parameter file has python syntax, we automatically append a ``.py`` extension to the custom file name.
If no argument is passed the parameter file will be generated in the current directory:

.. code-block:: console

    $ ls
    galapy_hyper_parameters.py

The generated file should be self-explanatory and has to be modified according to the fit the user has to perform.

.. Tip::
   The parameter file, for obvious reasons, needs at least to import the observational photometric dataset for running
   a parameter space sampling.
   Without this fundamental step, the program will crash as it does not know against what to compare the generated models.

   For testing the correct installation of the library, follow these 3 steps:

   * Generate a parameter file with a SFH model of choice:

     .. code-block:: console

	$ galapy-genparams -sfh insitu

   * open the generated parameter file with a text editor (e.g. :code:`$ emacs galapy_hyper_parameters.py`) and modify the
     first 4 parameters with the following minimal set-up:

     .. code-block:: python

	bands  = ['GOODS.b']
	fluxes = [1.0]
	errors = [1.0]
	uplims = [0]

     save and close the file.

   * By running the following command, the parameter-space sampling procedure will start using the Emcee sampler:

     .. code-block:: console

	$ galapy-fit galapy_hyper_parameters.py

There are 4 main blocks in the parameter file:

1. :ref:`import_obs_data`:
   where the user should load its observational data as well as define the photometric system;
2. :ref:`define_model`:
   define what SFH model and what SSP library to use, turn on/off components of the spectrum
   (i.e. radio support, x-ray support, AGN templates);
3. :ref:`fixed_and_free_parameters` (as well as their priors):
   each parameter of the model that will not be listed here will be set to its default value,
   otherwise, the user can here modify default values and choose what parameters to keep free;
4. :ref:`sampling_and_output`:
   decide which sampler to use and define the sampling hyperparameters (the default parameters
   should work most of the times), decide what and where to output.

.. tip::

   Keep in mind that the parameter file can be considered as a python module and, therefore,
   all operations that can be done in python are available within this file.
   It is, e.g., possible to import python packages that could be useful to manipulate the
   datasets or load them from system.
   It is also possible to run functions and define classes or whatever else the user might find
   useful for setting up properly their run.
   Some packages we always find useful to import at the beginning of the parameter file are

   .. code-block:: python

	import numpy as np
	import os, sys
   
Once the parameter file has been generated and properly modified, we can run

.. code-block:: console

   $ galapy-fit parameter_file.py [--serial/-s | --multiprocessing/-mp NCPU]

which will run the sampling and authomatically store the results, as specified
by the user in the parameter file.
NOTE THAT the two optional arguments regulate whether to run the sampling
serially or using shared-memory parallelism.
The default behaviour is to run parallely on all the available CPUs.

A thorough description of all the hyper-parameters that can be addressed when preparing the parameter file is provided in :ref:`param_file` page.

Fitting through the Python API
..............................

Coming soon.
