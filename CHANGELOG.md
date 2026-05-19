# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

| Release workflow: |
| --- |
|  **1.** Fill in the [Unreleased] section below. |
|  **2.** Run: bumpver update --patch   (or --minor / --major). This commits the version bump in `__init__.py` and creates a local tag. |
|  **3.** Rename [Unreleased] → [X.Y.Z] - YYYY-MM-DD and add a fresh [Unreleased]. |
|  **4.** Commit the changelog: git commit -m "update CHANGELOG for vX.Y.Z" |
|  **5.** Push: git push origin main && git push origin --tags |

> **Note:** pre-release tags are marked as vX.Y.Z-lw for "light-weight" on GitHub;
> the `-lw` tag suffix was used for pre-stable releases and dropped from v0.4.0 onwards.


## [Unreleased]

### Fixed
- Renamed AGN Fritz+2006 template parameters from `agn.ct`, `agn.al`, `agn.be`,
  `agn.ta`, `agn.rm`, `agn.ia` to `agn.template.ct` etc. — the old keys were
  silently ignored during fitting due to a clash in the nested-dict parsing of
  `ModelParameters`.

### Added
- `tests/Test_Run.py`: end-to-end test suite for the fitting pipeline covering
  `Run.initialize`, `Run.loglikelihood`, `Run.logprob`, all four SFH model
  variants, AGN/Radio/X-ray components, noise model, upper limits, custom
  filters, and the prior-transform utility.

### Changed
- Default sampler changed from `emcee` to `dynesty`.

---

## [0.5.3] - 2024-04-20

### Fixed
- Bug in `get_parameters_summary_statistics` when called with
  `stat_type='bestfit_and_interval'`.

### Added
- Installation instructions for Windows users.
- Extended installation guide in documentation.
- arXiv badge in README.
- Docstrings for the bandpass-transmission pybind11 interface.

### Changed
- Minimum required version of matplotlib raised to 3.6.

---

## [0.5.2] - 2024-02-19

### Fixed
- Missing import in sampling module.

### Added
- Tutorial: derivation of the attenuation curve in the results-analysis notebook.
- Docstrings for bandpass transmission (pybind11 interface).

### Changed
- Default dynesty sampler keywords and sampling keywords revised.
- `PMS` class inputs generalised to accept any iterable (not only numpy arrays); added documentation.

---

## [0.5.1] - 2024-02-15

### Fixed
- Bug in `Results` loading (`method='hdf5'`, `lightweight=False`) for runs that
  did not include a noise model.

### Added
- Tutorial for fitting with an interpolated SFH.
- Docstrings for `Cosmology`, `InterGalacticMedium`, `XRayBinaries`, and
  `Synchrotron` modules.

### Changed
- Photometric fluxes are now computed in the observer's frame.
- CI: GitHub Actions runners updated from v3 to v4.

---

## [0.5.0] - 2024-02-13

### Added
- Complete quick-start documentation and extended API reference.
- macOS classifier in package metadata.
- numpy declared as an explicit build dependency.

### Changed
- Project renamed to `galapy-fit` on PyPI.
- `bumpver` integrated for automated version management.

---

## [0.4.0] - 2024-02-09

### Added
- `LICENSE` file (GPLv3).
- Documentation pages for `CSP`, `GXY`, `ISM`, `NFF`, and `SFH` modules.
- Catalogue-loading tutorial.
- `import galapy` now imports the full public API.

### Changed
- New integration scheme for stellar mass in the CSP module: integration variable
  changed from ψ(t) dt to dM(t), implemented in C++. This changes the internal
  approximation for passive galaxies.

### Fixed
- Bug in interpolated SFH (now compliant with `GXY` interface).
- Rare crash (`'Warning in subtract'`) caused by negative times in the time
  tuple at very early ages.
- Multiprocessing on macOS Apple Silicon (Darwin).
- Bug in `_recursive_list_ssp_libs` introduced by the new SSP-library
  auto-detection.

---

## [0.3.0] - 2023-12-07

> **Note:** released as `v0.3.0-lw`;

### Added
- New treatment for upper limits with credible-interval support.
- Lightweight/heavyweight HDF5 storage modes, integrated in the entry points.
- New module `galapy.io` (depends on `h5py`).
- Abstract base class `GalaPyObject` in `galapy.internal.abc`.
- `ModelParameters.return_nested()` can now be called with no arguments to
  return all internally stored parameters (free + fixed).

### Fixed
- PAH template extrapolation: added `erf` damping to suppress discontinuities.

---

## [0.2.0] - 2023-10-12

> **Note:** released as `v0.2.0-lw`.

### Added
- Interpolated (empirical) SFH support.
- `galapy-genparams` entry point: `--SFH_model` flag to select the SFH model
  when generating the parameter file.
- Functions for reading TOPCAT-format and CSV catalogue files (`cat_to_dict`).
- New helper utilities (`shorten_string`, filter-name listing).

### Fixed
- Bug in stellar-mass integration (linear → logarithmic binning; different
  approximation for passive galaxies).
- Bug in interpolated SFH age integration.
- Long-standing spelling errors propagated across the package.
- Various minor plotting and array-casting fixes.

---

## [0.1.0] - 2023-05-15

> **Note:** released as `v0.1.0-lw`.
> Marks the transition to **pybind11** C++ bindings (replacing the previous
> CPython extension strategy) and the migration to the new internal database
> format. Users upgrading from `v0.0.3-lw` must re-download the database.

### Added
- Complete renewal of `galapy.sampling.Run`: now supports parallel runs via
  `multiprocessing` and both `emcee` and `dynesty` samplers.
- Upper-limit support in SED-plotting functions.
- AGN safety check: `fAGN >= 1.0` is now rejected at parameter-setting time.
- PAH suppression for wavelengths beyond `lambda_Wien`.

### Fixed
- Bug in `set_parameters`: now raises an exception when the requested age
  exceeds the age of the universe at the given redshift.
- Bug in SED plotting when redshift is a free parameter.
- Age initialisation conflicting with the selected redshift in the run script.

---

## [0.0.3] - 2023-02-25

> **Note:** released as `v0.0.3-lw`.
> Last release using the CPython C-extension strategy; `v0.1.0-lw` supersedes
> it with pybind11 bindings.

### Added
- Inter-Galactic Medium (IGM) module implementing the Inoue et al. (2014) model.
- Noise modelling (`galapy.Noise`): `CalibrationError` class.
- Abstract interfaces in `galapy.internal.abc`; `Model` base class required for
  all physical models.
- `galapy.Handlers` module (renamed from `GalaxyParameters`); new unified
  `ModelParameters` handler for galaxy + noise.
- Docstrings for the `Galaxy` module; physical-units documentation page.
- `pytest` added to requirements.

### Fixed
- Bug in `Handlers.GXYParameters`.
- Bug in `ibstree_interface` integrate function.
- Likelihood function: NaN from the noise term now returns `-inf` instead of
  propagating.
- Bug in parameter-setting when changing redshift with IGM disabled.

---

## [0.0.2] - 2023-01-18

> **Note:** released as `v0.0.2-lw`.
> Near-complete implementation checkpoint; the noise modelling term was still
> pending and was completed in `v0.0.3-lw`.

### Added
- `Test_GxyPar` test module.
- `cat_to_dict` function for reading TOPCAT-like and CSV files.

### Fixed
- Bug in `galapy-fit` entry point: SFH model was not correctly selected for
  models other than the default.
- Bug in `get_avgAtt()`.
- Various minor array-casting and import fixes.

---

## [0.0.1] - 2022-08-29

> **Note:** released as `v0.0.1-lw` — initial pre-release.

- First public release of the GalaPy spectral modelling library.
