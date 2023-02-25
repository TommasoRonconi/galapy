Physical Units
==============

GalaPy, being an API for making astrophysical computations, assumes some of the quantities passed as paraemeters and returned from functions are in given physical units.

Even though the assumed physical units are also reported throughout the `documentation`_ and the `examples`_, we provide here a table wrapping-up our physical system.

+-----------------+--------------------------------------------+
| **Quantity**    | **Unit**                                   |
+=================+============================================+
| Wavelenght      | Angstroms :math:`[\mathring{A}]`           |
+-----------------+--------------------------------------------+
| Time/Age        | Years :math:`[\text{yr}]`                  |
+-----------------+--------------------------------------------+
| Metallicity     | Metals/Total Fraction                      |
+-----------------+--------------------------------------------+
| Temperature     | Kelvin degrees :math:`[K]`                 |
+-----------------+--------------------------------------------+
| Luminosity      | :math:`[\text{erg}\ s\ \mathring{A}^{-1}]` |
+-----------------+--------------------------------------------+
| Energy          | :math:`[\text{erg}]`                       |
+-----------------+--------------------------------------------+

.. note::
   We plan to extend and make our physical unit system more flexible.
   Therefore this page might experience major modifications in future versions of the API.
   We will try our best to make the eventual future transition painless,
   ensuring backward-compatibility whenever possible.

.. _documentation: ...
.. _examples: ...
