GalaPy API Documentation
========================

Here the documentation of modules, sub-packages, classes and functions delivered with the GalaPy API can be found.

.. warning::

   Coverage is not complete yet.

Galaxy Components modules
-------------------------

These modules wrap the different components that build-up a :code:`galapy.Galaxy.GXY` and derived objects.
Each of the galactic components can be used independently and, eventually, extended with further functionalities.

.. autosummary::
   :toctree: _generated

   galapy.Galaxy
   galapy.StarFormationHistory
   galapy.CompositeStellarPopulation
   galapy.InterStellarMedium
   galapy.ActiveGalacticNucleus
   galapy.NebularFreeFree
   galapy.Synchrotron
   galapy.XRayBinaries
   galapy.InterGalacticMedium
   galapy.Cosmology

Photometric System
------------------
   
.. autosummary::
   :toctree: _generated

   galapy.BandpassTransmission
   galapy.PhotometricSystem

Noise modelling
---------------

.. autosummary::
   :toctree: _generated

   galapy.Noise

Sampling
--------

.. autosummary::
   :toctree: _generated

   galapy.Handlers
   galapy.sampling.Observation
   galapy.sampling.Results
   galapy.sampling.Sampler
   galapy.sampling.Statistics

Analysis
--------

.. autosummary::
   :toctree: _generated

   galapy.analysis.plot
   galapy.analysis.funcs
   
Helpers Functions
-----------------

.. autosummary::
   :toctree: _generated

   galapy.internal.utils
   galapy.internal.constants
   
