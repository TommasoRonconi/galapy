Physical Units (used in functions Input/Output)
===============================================

GalaPy, being an API for making astrophysical computations, assumes some of the quantities passed as parameters and returned from functions are in given physical units.

Even though the assumed physical units are also reported throughout the `documentation`_ and the `examples`_, we provide here a table wrapping-up our physical system.

+-----------------+--------------------------------------------+
| **Quantity**    | **Unit**                                   |
+=================+============================================+
| Wavelength      | Angstroms :math:`[\mathring{A}]`           |
+-----------------+--------------------------------------------+
| Frequency       | Hertz :math:`[\text{Hz}]`                  |
+-----------------+--------------------------------------------+
| Time/Age        | Years :math:`[\text{yr}]`                  |
+-----------------+--------------------------------------------+
| Metallicity     | Metals/total mass                          |
+-----------------+--------------------------------------------+
| Temperature     | Kelvin degrees :math:`[K]`                 |
+-----------------+--------------------------------------------+
| Luminosity      | Solar luminosities :math:`[L_\odot]`       |
+-----------------+--------------------------------------------+
| Energy          | :math:`[\text{erg}\cdot s]`                |
+-----------------+--------------------------------------------+
| Energy Flux     | milli-Jansky :math:`[\text{m}Jy]`          |
+-----------------+--------------------------------------------+
| Mass            | Solar masses :math:`[M_\odot]`             |
+-----------------+--------------------------------------------+


Compatibility with other libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NumPy (and SciPy)
-----------------

All objects and functions in GalaPy are compliant with NumPy and SciPy (i.e. input and output values are ND-arrays and/or can be transformed into ND-arrays).

Taking the ``galapy.Galaxy.GXY`` class as an example:

.. code:: ipython

   In [1]: from galapy.Galaxy import GXY
   In [2]: gxy = GXY(age=1.e+10, redshift=0.1)

When we compute the emission from an object of type ``GXY``, the output luminosity is a NumPy ND-array:

.. code:: ipython

   In [3]: L = gxy.get_emission()
   In [4]: type(L)
   Out[4]: numpy.ndarray

This means that, knowing that luminosities are given in Solar masses, :math:`L_\odot`, if we want to convert the array in some other units (e.g. Watts, :math:`\text{W}`),
we only have to multiply :code:`L` by the conversion factor (i.e. :math:`L_\odot = 3.828\times10^{26} \text{Watt}`).
Since :code:`L` is a NumPy array, multiplication by a constant is automatically broadcasted to the whole array:

.. code:: ipython

   In [5]: L * 3.828e+26
   Out[5]: array([2.2075896e+00, 1.4990568e+10, 8.4263495e+11, ..., 3.2180310e+17,
	          8.5173570e+06, 7.3632465e+06], dtype=float32)

Furthermore, most of the functions (see `documentation`_ for in-detail description), accept both scalars and arrays as inputs.
If we want to, e.g., compute the stellar mass at different ages we can either compute it at some given time:

.. code:: ipython

   In [6]: gxy.sfh.Mstar( 1.e+9 )
   Out[6]: 24756528290.883263

Or we can pass an array of ages:

.. code:: ipython

   In [7]: import numpy
   In [8]: gxy.sfh.Mstar( numpy.lospace( 7, 10, 10 )
   Out[8]: array([1.42748587e+07, 6.24427121e+07, 2.65490510e+08, 1.06947300e+09,
                  3.87568658e+09, 1.15193272e+10, 2.47565283e+10, 3.54007409e+10,
                  3.62805578e+10, 3.33419941e+10])

astropy
-------

Some users might be familiar with `astropy`_ and its functionalities to manage physical quantities, units and conversions.
Even though astropy **is not a GalaPy dependency**, GalaPy can nevertheless communicate easily with this external and popular library.
This interoperability proves particularly useful in the context of units conversion, thanks to the large library of units and constants available in astropy and to its intuitive interface.

We can first convert the luminosity array :code:`L` computed above into an astropy object:
		  
.. code:: ipython

   In [9]: import astropy.units as u
   In [9]: L *= u.L_sun
   In [10]: type(L)
   Out[10]: astropy.units.quantity.Quantity

It is now possible to convert it to whatever other luminosity unit. So if we want to convert the value in Watts, we would just do something like this:

.. code:: ipython

   In [11]: L.to(u.Watt)
   Out[11]:

:math:`[2.2075896, 1.4990568\times10^{10}, 8.4263495\times10^{11}, ..., 3.218031\times10^{17}, 8517357, 7363246.5]\ \text{W}`

Consistently, we can also use astropy to easily pass from wavelengths to frequencies, transforming with

.. math::

   \nu = \dfrac{c}{\lambda}

where :math:`c` is the speed of light, :math:`\lambda` is the wavelength and :math:`\nu` is the frequency.

.. code:: ipython

   In[12]: lambda = gxy.wl()
   In[13]: from astropy.constants import c
   In[14]: nu = c.to(u.AA*u.Hz)/lambda

where we have used the speed of light from astropy module :code:`astropy.constants` and converted it in :math:`[\mathring{A}\cdot\text{Hz}]` units.
   
.. _documentation: ...
.. _examples: ...
.. _astropy: https://docs.astropy.org/en/stable/index.html
