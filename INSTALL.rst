.. _install_guide:

Installation Guide
==================

The preferred method to install is through the Python Package Index.

It is anyways possible for developers to install from source cloning the repository from GitHub or by downloading a specific version tarball.

.. note::

   Windows installation is not supported yet, neither with pip nor from source.
   We accept the contribution of resourceful developers who might want to try to implement this extension.

From the Python Package Index (PyPI)
------------------------------------

GalaPy is available on PyPI for Linux and MacOS.
With an internet connection, install the library by calling

.. code-block:: console

   $ pip install galapy-fit

or

.. code-block:: console

   $ python -m pip install galapy-fit

PyPI will automatically install all the necessary pre-requisites.

After Install
-------------
  
GalaPy, in some of its components (e.g. SSP tables, PAH template), makes use of pre-computed functions that are
available in the official database (`galapy_database`_). When one of the files in the database is accessed for the
first time it will authomatically be downloaded into the user's filesystem
(in the default location :code:`$HOME/.galapy/galapy_database`).
This will of course require an internet connection and can partially slow down the computations.
We therefore **strongly recommend**, prior to first run, to download all the database by running

.. code-block:: console

   $ galapy-download-database

This operation has to be run only once.

From source
-----------

Installation from source can be done either by cloning the GitHub repository of the project or by downloading the
specified release (this allows a higher control on the version to install).

* Clone the repository and enter in the directory

  .. code-block:: console
     
     $ git clone git@github.com:TommasoRonconi/galapy.git
     $ cd galapy

* Download the tarball for version ``vX.Y.Z`` (substitute the ``X.Y.Z`` string with the version chosen), enter in the directory

  .. code-block:: console

     $ curl -L https://github.com/TommasoRonconi/galapy/archive/refs/tags/vX.Y.Z.tar.gz -o vX.Y.Z.tar.gz
     $ tar -xzvf vX.Y.Z.tar.gz
     $ cd galapy-X.X.X

Installation in both cases can be done using ``pip`` locally, from inside the directory containing the source of the project
(i.e. the directory where the ``setup.py`` and ``pyproject.toml`` files are located), call

.. code-block:: console

   $ pip install .

or

.. code-block:: console

   $ python -m pip install .

By using ``pip`` instead of the classic (and shortly deprecated) ``python setup.py install`` the correct installation
of the library prerequisites is also guaranteed.

.. note::

   Even when installed from source, GalaPy still requires downloading the database (see previous section).

.. tip::

   To test the installation, run

   .. code::

      $ pytest tests/*.py

   and check that all the tests have passed.

Pre-requisites
--------------

A table with the dependencies of the library follows. We also provide a short description of the package and the relevant link for download. 

+------------------------------------------------------+------------------------------------------------------------------------------------+
| `NumPy <https://pypi.org/project/numpy/>`_           | Fundamental package for array computing in Python                                  |
+------------------------------------------------------+------------------------------------------------------------------------------------+
| `SciPy <https://pypi.org/project/scipy/>`_           | Fundamental algorithms for scientific computing in Python                          |
+------------------------------------------------------+------------------------------------------------------------------------------------+
| `emcee <https://pypi.org/project/emcee/>`_           | The Python ensemble sampling toolkit for MCMC                                      |
+------------------------------------------------------+------------------------------------------------------------------------------------+
| `dynesty <https://pypi.org/project/dynesty/>`_       | A dynamic nested sampling package for computing Bayesian posteriors and evidences. |
+------------------------------------------------------+------------------------------------------------------------------------------------+
| `setuptools <https://pypi.org/project/setuptools/>`_ | Easily download, build, install, upgrade, and uninstall Python packages            |
+------------------------------------------------------+------------------------------------------------------------------------------------+
| `matplotlib <https://pypi.org/project/matplotlib/>`_ | Python plotting package                                                            |
+------------------------------------------------------+------------------------------------------------------------------------------------+
| `requests <https://pypi.org/project/requests/>`_     | Python HTTP for Humans.                                                            |
+------------------------------------------------------+------------------------------------------------------------------------------------+
| `pytest <https://pypi.org/project/pytest/>`_         | simple powerful testing with Python                                                |
+------------------------------------------------------+------------------------------------------------------------------------------------+

.. _galapy_database: https://github.com/TommasoRonconi/galapy_database
