Quickstart
----------

.. graphviz::
   :align: center
   :caption: Flux diagram summarising the main steps required to run GalaPy.

    digraph G {
          size = "7.5,7.5";
          compound=true;
	  nodesep=.05;
	  fontname="Helvetica";
	  newrank=true;
	  
          a [shape=box, color="#61C0E8", style=filled, label="Install the library", href="https://galapy.readthedocs.io/en/latest/general/install_guide.html", target="_blank"];
          b [shape=box, color="#61C0E8", style=filled, label="Download the database", href="https://galapy.readthedocs.io/en/latest/general/install_guide.html#after-install", target="_blank"];
          c [shape=box, color="#B2C600", style=filled, label="Generate the parameter file", href="https://galapy.readthedocs.io/en/latest/guides/parameter_file.html", target="_blank" ];
          d [shape=box, color="#B2C600", style=filled, label="Run the fit", href="https://galapy.readthedocs.io/en/latest/guides/photometric_fit.html", target="_blank"];
          e [shape=box, color="#E78A61", style=filled, label="Results analysis", href="https://galapy.readthedocs.io/en/latest/notebooks/results_analysis.html", target="_blank"];
	  
	  install_command [shape=box, color="#61C0E8", label="pip install galapy", fontname="Verbatim", fontsize=12];
	  download_command [shape=box, color="#61C0E8", label="galapy-download-database", fontname="Verbatim", fontsize=12];
	  genparams_command [shape=box, color="#B2C600", label="galapy-genparams", fontname="Verbatim", fontsize=12];
	  fit_command [shape=box, color="#B2C600", label="galapy-fit", fontname="Verbatim", fontsize=12];

	  subgraph cluster0 {
	     color="black";
	     style=dashed
             a -> b [maxlen=1.];
	     label="GalaPy Library Set-Up";
	  }
	  
          subgraph cluster1 {
             color=black;
	     c -> d -> e;
	     label="GalaPy Workflow";
          }
	  
	  { rank = same; "a"; "install_command"; }
	  { rank = same; "b"; "download_command"; }
	  { rank = same; "c"; "genparams_command"; }
	  { rank = same; "d"; "fit_command"; }


          b -> c [ltail=cluster0,lhead=cluster1, minlen=2., style=dashed];

       }

.. tip::

   Click (or tap) on the filled boxes in the flux diagram above to be re-directed to 
   the relative documentation page with in-detail explanation of the step.

The GalaPy workflow can be summarised in two blocks of steps, the first one (dubbed "GalaPy library set-up") has to
be performed only once:

1. Install the library through `pip`:

   .. code:: console

      $ pip install galapy

2. Download the database through the dedicated terminal command:

   .. code:: console

      $ galapy-download-database

   Note that the official `GalaPy database <https://github.com/TommasoRonconi/galapy_database>`_ occupies approximately
   300 MB and will take tens of seconds up to few minutes to download (depending on internet connection).

This is enough to have a working installation of GalaPy.

In order to fit a photometric dataset the steps to follow are the ones detailed in the "GalaPy" workflow block above:

1. Generate the parameter file and modify it accordingly to the requirements of your dataset.
   As a term of comparison, to generate a parameter file to fit panchromatic data with an In-Situ SFH model,
   call the command

   .. code:: console

      $ galapy-genparams --name insitu_params --SFH_model insitu

   which will generate a file called ``insitu_params.py``.

   By opening this text file change the entries for the following parameters

   .. code:: python

      bands  = None
      fluxes = None
      errors = None

   with the values of your observed source (see :ref:`import_obs_data` for more details).

   .. note::
      Fluxes have to be given in :math:`mJy` (milli-Jansky), see :ref:`physical_units`
      for further details on the physical units assumed in the library.

   .. tip::
      A complete list of the bands available in the database is print on screen by
      calling function :py:func:`galapy.PhotometricSystem.print_filters`.

   Another section that it could be useful to adapt to the user's requirements is the
   :ref:`sampling_and_output` section of the parameter file, where the sampler and run
   specifications are chosen.

2. Run the fitting algorithm by calling the terminal command, using the parameter file
   generated at the previous step:

   .. code:: console

      $ galapy-fit insitu_params.py

   This will start the fitting procedure that will take (depending on the
   statistical properties of the sample, installation and system) from a bunch
   of minutes to several tens of minutes.
   When the program finishes its run, a file (or several files)
   with extension ``.galapy.hdf5`` is generated in the location specified in the
   parameter file (if no choice is made it will be saved in the working directory).

3. Analyse the dataset by following the tutorial in
   `Results analysis <https://galapy.readthedocs.io/en/latest/notebooks/results_analysis.html>`_

For further details on each of these steps, please follow the links in the filled boxes of the
flux diagram at the beginning of this page and in the following text.
