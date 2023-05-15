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

The first hyperparameter to set defines the photometric system, in the parameter file generated it will appear as

.. code-block:: python
   
   filters = None

The ``None`` value is just a place-holder and should be modified.
There are two possible choices:

1. build the photometric system from transmission filters already available in the database.
   In this case the hyperparameter should be set to a list or a tuple listing all the filters
   chosen:

   .. code-block:: python

	filters = ( 'GOODS.i', 'GOODS.v', 'GOODS.z', 'GOODS.b' )
  
2. use a custom set of filters. In this case the hyperparameter has to receive a nested
   dictionary properly formatted. This means each element in the dictionary should have a
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

Once the photometric system has been set, the observational dataset has to be defined, this is done by setting the following python iterables of values:
   
.. code-block:: python
   
   bands = None
   fluxes = None
   errors = None
   uplims = None

Also in this case, the ``None`` values are placeholders that should be modified by the user.
In particular, the sequences of values can be a ``numpy array``, a python ``list``, a ``set`` or a ``tuple``, it is nevertheless necessary that they all have all the same dimensions.

1. ``bands``: a sequence of 

.. _define_model:

Define the physics of the galaxy model
......................................

.. _fixed_and_free_parameters:

Choose the fixed and free parameters
....................................

.. _sampling_and_output:

Sampling and output format choices
..................................
