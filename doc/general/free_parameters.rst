.. _tunable_params:

Free Parameters of the Galaxy model
===================================

We provide a summary table collecting all the free parameters available for fitting with the tools delivered with GalaPy (both using the command line and the custom fit functionalities). 

All the parameters that can be tuned by sampling the parameter space or adjusted when choosing a particular model over some other are considered *free parameters*.
Whether to fix a parameter to some particular value or to let it vary (and therefore, to consider it one additional dimension of the parameter space to be sampled),
is decided by modifying the ``galaxy_parameters`` dictionary in the parameter-file (see also :ref:`fixed_and_free_parameters`).

The table below collects all of these free-parameters and, along with the keyword used to univocally refer to them, provides the default value assumed internally by GalaPy,
the mathematical symbol and short description of the parameter, consistently with the paper presenting the library.

.. Tip::

   The table has the purpose to guide the user when choosing what parameters to let vary, and within what ranges.
   If no value nor interval of values is provided to GalaPy, the library will fix the value of the parameter to the
   quantity reported in the *Default Value* column.

+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| **Parameter**       | **Default Value**                | **Math Symbol :**                    **Description**                                                                         |
+=====================+==================================+==============================================================================================================================+
| ``age``             | :math:`10^6\ \text{years}`       | :math:`\tau`                       : Age of the galaxy                                                                       |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``redshift``        | :math:`0`                        | :math:`z`                          : Redshift of the galaxy                                                                  |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| **Star Formation History**                                                                                                                                                            |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.model``       | ``insitu``                       | :math:`\psi(\tau)`                 : SFH model: one among the keywords listed in :py:class:`galapy.StarFormationHistory.SFH` |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.tau_quench``  | :math:`2\cdot10^{20}\ \text{yr}` | :math:`\tau_\text{quench}`         : Age of the abrupt quenching                                                             |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| Constant (``sfh.model = 'constant'``)                                                                                                                                                 |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.psi``         | :math:`1\ M_\odot/\text{yr}`     | :math:`\psi_0`                     : Value of the constant SFR                                                               |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.Mdust``       | :math:`10^8\ M_\odot`            | :math:`M_\text{dust}`              : Total dust mass in galaxy at the given age                                              |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.Zgxy``        | :math:`0.01`                     | :math:`Z_\text{gxy}`               : Metallicity of all phases in galaxy at the given age                                    |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| Delayed Exponential (``sfh.model = 'delayedexp'``)                                                                                                                                    |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.psi_norm``    | :math:`1\ M_\odot/\text{yr}`     | :math:`\psi_\text{norm}`           : Normalisation                                                                           |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.k_shape``     | :math:`0.2`                      | :math:`\kappa`                     : Shape parameter of the early evolution                                                  |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.tau_star``    | :math:`10^8\ \text{yr}`          | :math:`\tau_\star`                 : Characteristic timescale                                                                |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.Mdust``       | :math:`10^8\ M_\odot`            | :math:`M_\text{dust}`              : [same as for constant SFH]                                                              |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.Zgxy``        | :math:`0.01`                     | :math:`Z_\text{gxy}`               : [same as for constant SFH]                                                              |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| Log-Normal (``sfh.model = 'lognormal'``)                                                                                                                                              |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.psi_norm``    | :math:`100\ M_\odot/\text{yr}`   | :math:`\psi_\text{norm}`           : Normalisation                                                                           |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.sigma_star``  | :math:`2`                        | :math:`\sigma_\star`               : Characteristic width                                                                    |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.tau_star``    | :math:`3\cdot10^8\ \text{yr}`    | :math:`\tau_\star`                 : Peak age                                                                                |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.Mdust``       | :math:`10^8\ M_\odot`            | :math:`M_\text{dust}`              : [same as for constant SFH]                                                              |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.Zgxy``        | :math:`0.01`                     | :math:`Z_\text{gxy}`               : [same as for constant SFH]                                                              |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| Interpolated (``sfh.model = 'interpolated'``)                                                                                                                                         |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.Mdust``       | :math:`10^8\ M_\odot`            | :math:`M_\text{dust}`              : [same as for constant SFH]                                                              |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.Zgxy``        | :math:`0.01`                     | :math:`Z_\text{gxy}`               : [same as for constant SFH]                                                              |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| In-Situ (``sfh.model = 'insitu'``)                                                                                                                                                    |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.psi_max``     | :math:`100\ M_\odot/\text{yr}`   | :math:`\psi_\text{max}`            : Normalisation                                                                           |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``sfh.tau_star``    | :math:`3\cdot10^8\ \text{yr}`    | :math:`\tau_\star`                 : Characteristic timescale                                                                |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| **Inter Stellar Medium**                                                                                                                                                              |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.f_MC``        | :math:`0.5`                      | :math:`f_\text{MC}`                : Fraction of dust in the MC phase                                                        |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.norm_MC``     | :math:`10^2`                     | :math:`\mathcal{C}_V^\text{MC}`    : Normalisation of the MC attenuation in the visible band                                 |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.N_MC``        | :math:`10^3`                     | :math:`N_\text{MC}`                : Number of MCs in the galaxy                                                             |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.R_MC``        | :math:`10\ \text{pc}`            | :math:`R_\text{MC}`                : Average radius of a MC                                                                  |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.tau_esc``     | :math:`10^7\ \text{yr}`          | :math:`\tau_\text{esc}`            : Time required by stars to start escaping their MC                                       |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.dMClow``      | :math:`1.3`                      | :math:`\delta_\text{MC}^\text{l}`  : Extinction power-law index at wavelength :math:`\lesssim100 \mu m~(10^6 \mathring{A})`  |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.dMCupp``      | :math:`1.6`                      | :math:`\delta_\text{MC}^\text{u}`  : Extinction power-law index at wavelength :math:`\gtrsim100 \mu m~(10^6 \mathring{A})`   |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.norm_DD``     | :math:`1`                        | :math:`\mathcal{C}_V^\text{DD}`    : Normalisation of the DD extinction in the visible band                                  |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.Rdust``       | :math:`10^3\ \text{pc}`          | :math:`R_\text{DD}`                : Radius of the diffuse dust region embedding stars and MCs                               |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.f_PAH``       | :math:`0.2`                      | :math:`f_\text{PAH}`               : Fraction of the total DD luminosity radiated by PAH                                     |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.dDDlow``      | :math:`0.7`                      | :math:`\delta_\text{DD}^\text{l}`  : Extinction power-law index at wavelength :math:`\lesssim100 \mu m~(10^6 \mathring{A})`  |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``ism.dDDupp``      | :math:`2`                        | :math:`\delta_\text{DD}^\text{u}`  : Extinction power-law index at wavelength :math:`\gtrsim100 \mu m~(10^6 \mathring{A})`   |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| **Nebular Free-Free**                                                                                                                                                                 |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``nff.Zi``          | :math:`1`                        | :math:`Z_i`                        : Average atomic number of ions                                                           |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| **Synchrotron**                                                                                                                                                                       |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``syn.alpha_syn``   | :math:`0.75`                     | :math:`\alpha_\text{syn}`          : Spectral index                                                                          |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``syn.nu_self_syn`` | :math:`0.2\ G\text{Hz}`          | :math:`\nu_\text{self}`            : Self-absorption frequency                                                               |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| **Active Galactic Nucleus**                                                                                                                                                           |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``agn.fAGN``        | :math:`10^{-3}`                  | :math:`f_\text{AGN}`               : AGN fraction                                                                            |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``agn.ct``          | :math:`40\ \text{deg}`           | :math:`\Theta`                     : Torus half-aperture angle                                                               |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``agn.al``          | :math:`0`                        | :math:`\alpha`                     : Density parameter (exponential part)                                                    |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``agn.be``          | :math:`-0.5`                     | :math:`\beta`                      : Density parameter (power-law part)                                                      |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``agn.ta``          | :math:`6`                        | :math:`\tau_{9.7}^\text{AGN}`      : Optical depth at :math:`9.7 \mu m`                                                      |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``agn.rm``          | :math:`60`                       | :math:`R_\text{torus}^\text{AGN}`  : Radial ratio of the torus                                                               |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| ``agn.ia``          | :math:`0.001\ \text{deg}`        | :math:`\Psi_\text{los}^\text{AGN}` : Inclination angle                                                                       |
+---------------------+----------------------------------+------------------------------------------------------------------------------------------------------------------------------+

.. note::

   All the parameters that do not regulate directly the physics of the galaxy model or of the noise model and, in general, all the parameters that do not appear in the
   table above, are instead considered **hyper-parameters**.
