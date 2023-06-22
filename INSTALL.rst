Installation Guide
==================


From the Python Package Index (PyPI)
------------------------------------

.. code-block:: console

   $ pip install galapy

PyPI deals with pre-requisites.

From source
-----------

Clone the repository

.. code-block:: console

   $ git clone git@github.com:TommasoRonconi/galapy.git
   $ cd galapy

Download the tarball

.. code-block:: console

   $ wget ...
   $ tar -xzvf ...
   $ cd galapy-X.X.X

Install

.. code-block:: console

   $ python setup.py install

Pre-requisites
^^^^^^^^^^^^^^

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

After Install
-------------
  
.. note::
   GalaPy, in some of its components (e.g. SSP tables, PAH template), makes use of pre-computed functions that are
   available in the official database (`galapy_database`_). When one of the files in the database is accessed for the
   first time it will authomatically be downloaded into the user's filesystem
   (in the default location :code:`$HOME/.galapy/data`).
   This will of course require an internet connection and can partially slow down the computations.
   We therefore suggest, prior to first run, to download all the database by running

.. code-block:: console

   $ galapy-download-database


.. _galapy_database: https://github.com/TommasoRonconi/galapy_database
