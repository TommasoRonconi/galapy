.. _custom_likelihood:

Using a custom likelihood
=========================

By default ``galapy-fit`` scores every proposed parameter set with the built-in
Gaussian log-likelihood, :func:`galapy.sampling.Run.loglikelihood`. For most
photometric fits this is the right choice: it assumes the flux measurements are
independent and normally distributed around the model, with the quoted
1-:math:`\sigma` errors, and it treats non-detections according to the
``method_uplims`` hyperparameter.

There are cases in which it is not what you want, for instance

* you have an external, independent estimate of a physical property of the
  source (say, the total stellar mass from spectroscopy or dynamics) and
  want the photometric fit to be consistent with it;
* your photometry has outliers and you would rather use a heavier-tailed
  distribution than a Gaussian, so that a single discrepant band does not
  dominate the fit;
* the errors are correlated across bands and you have a covariance matrix.

For these situations the likelihood can be replaced wholesale from the
parameter file, through the ``loglikelihood`` hyperparameter.

.. warning::

   Replacing the likelihood changes the statistical meaning of the whole run:
   the posterior, the errors and the evidence all depend on it. This is an
   advanced feature and it is entirely your responsibility to make sure the
   function you provide is a valid sampling distribution for your data.

Properties of the custom loglikelihood 
--------------------------------------

A log-likelihood usable by galapy is any callable with the signature

.. code-block:: python

   def my_loglikelihood ( par, state, **kwargs ) :
       ...
       return llike

where

``par``
   1-D :class:`numpy.ndarray` with the current values of the free parameters,
   ordered as ``state.handler.par_free``.

``state``
   The :class:`galapy.sampling.Run.PipelineState` of the run. It gives access
   to everything the likelihood may need:

   * ``state.handler`` — the :class:`galapy.Handlers.ModelParameters`
     instance; ``state.handler.return_nested(par)`` converts the flat vector
     into the nested ``{'galaxy': {...}, 'noise': {...}}`` dictionaries
     expected by the models;
   * ``state.model`` — the :class:`galapy.Galaxy.PhotoGXY` galaxy model;
   * ``state.noise`` — the noise model, or ``None`` if the run has none;
   * ``state.data`` — the :class:`galapy.sampling.Observation.Observation`
     being fitted, with the ``fluxes``, ``errors`` and ``uplims`` arrays and
     the photometric system ``pms``.

``**kwargs``
   Extra keyword arguments forwarded by the sampler through ``logl_kw``.
   Currently this is only ``method_uplims``, but **always accept**
   ``**kwargs``: a function without it will raise a :class:`TypeError` as soon
   as sampling starts, and galapy warns about this at load time.

The function **must**

1. return a **scalar** — not a length-1 array, not a ``nan``;
2. return ``-numpy.inf`` whenever the model rejects the parameters, i.e. when
   ``state.model.set_parameters`` raises :class:`RuntimeError` (this happens
   routinely, e.g. when the proposed age exceeds the age of the Universe at
   the proposed redshift, so it is not an error condition — it is how the
   sampler learns the boundary of the allowed region);
3. return ``-numpy.inf`` whenever the computed value is not finite;
4. evaluate the model inside a ``numpy.errstate(all='ignore')`` context —
   overflows and divisions by zero are expected in the far tails of the prior
   and must not be allowed to raise or to spam the output.

Points 2–4 are not optional bookkeeping: a likelihood that raises, or that
returns ``nan``, will either abort a run hours into it or silently corrupt the
posterior.

A complete example: adding a stellar-mass constraint
----------------------------------------------------

Take the first case of the list above: an external estimate of the total
stellar mass, :math:`\log_{10} M_\star = 10.65 \pm 0.15`, that the photometric
fit should account for. The stellar mass is already computed by the galapy
model at every proposed sample (it is one of the default derived quantities of
a :class:`~galapy.sampling.Results.Results` object), so the constraint is just
one extra Gaussian term — in :math:`\log_{10}`, since stellar-mass estimates
carry log-normal uncertainties — on top of the photometric likelihood.

There is no need to re-implement the photometric part: the built-in
:func:`~galapy.sampling.Run.loglikelihood` is an importable function honouring
the same contract, so the custom likelihood can simply call it and add the new
term. This automatically inherits all of the mandatory handling above, plus
the treatment of upper limits and of the noise model. The measured values are
*arguments* of the function, so that the same module can serve any number of
objects — each parameter file binds its own numbers, as shown in the next
section:

.. code-block:: python

   # my_likelihoods.py
   import numpy
   from galapy.sampling.Run import loglikelihood as photometric_loglikelihood

   def mstar_loglikelihood ( par, state, logmstar_obs, logmstar_err, **kwargs ) :
       """Built-in photometric likelihood + external stellar-mass constraint,
       Gaussian in log10( Mstar ).

       logmstar_obs, logmstar_err : the external estimate, log10( Mstar/Msun )
       and its 1-sigma uncertainty in dex; bound per object in the parameter
       file with functools.partial.
       """

       # 1. photometric term: the built-in Gaussian likelihood. It sets the
       #    model parameters, applies the noise model, treats the upper
       #    limits according to method_uplims, and returns -inf for
       #    parameter sets the model rejects.
       llike = photometric_loglikelihood( par, state, **kwargs )
       if not numpy.isfinite( llike ) :
           return -numpy.inf

       # 2. stellar-mass term: the call above has already set the proposed
       #    parameters on the model, so it can be queried directly.
       with numpy.errstate( all = 'ignore' ) :
           logM = numpy.log10( state.model.sfh.Mstar( state.model.age ) )
           chi  = ( logM - logmstar_obs ) / logmstar_err
           llike += -0.5 * chi * chi

       # 3. never return a non-finite value other than -inf
       return llike if numpy.isfinite( llike ) else -numpy.inf

The same composition pattern works for any external constraint that can be
written as a function of the model state — any quantity reachable from
``state.model`` can enter the extra term.

To instead change the *metric* itself (a Student-t distribution for outlier
resistance, a full covariance matrix, ...), the photometric part has to be
re-implemented: start from the source of
:func:`galapy.sampling.Run.loglikelihood`, which is deliberately the reference
implementation of the contract, and see
:func:`galapy.sampling.Statistics.gaussian_loglikelihood` for how the three
``method_uplims`` strategies handle non-detections.

Where to put the function
-------------------------

.. important::

   The likelihood **must live in a module that Python can import by name.**

When the sampling is parallelised, the likelihood is sent to the worker
processes by pickling it, and pickling a function stores a *reference* to it —
its module and qualified name — rather than its code. The worker then imports
that path to rebuild the function. Consequently the following **do not work**
in parallel runs:

* lambdas;
* functions defined inside another function (closures);
* functions defined directly in the parameter file — the parameter file is
  loaded by path as a throw-away module named ``hyper_parameters``, which a
  worker process cannot import.

galapy detects all three cases and emits a warning when the parameter file is
read, but it cannot detect every possibility, so the rule to follow is simply:
put the likelihood in a normal ``.py`` module that is either installed or
reachable through ``PYTHONPATH``, and import it in the parameter file.

A :func:`functools.partial` is the exception that makes per-object data easy:
it is pickled by value, so a partial built *anywhere* — including the
parameter file — works, as long as the function it wraps is importable (the
checks above are applied to the wrapped function).

Since ``PYTHONPATH`` is inherited by spawned workers, the simplest working
setup is

.. code-block:: console

   $ export PYTHONPATH=/path/to/my/module:$PYTHONPATH
   $ galapy-fit my_parameter_file.py

Wiring it in the parameter file
-------------------------------

In the sampler section of the parameter file generated by ``galapy-genparams``
you will find

.. code-block:: python

   loglikelihood = None

``None`` selects the built-in Gaussian likelihood. To use your own, import it,
bind the per-object data, and assign the result:

.. code-block:: python

   from functools import partial
   from my_likelihoods import mstar_loglikelihood

   loglikelihood = partial( mstar_loglikelihood,
                            logmstar_obs = 10.65, logmstar_err = 0.15 )

Assign a callable — do **not** call it. Binding the data with
:func:`functools.partial` *in the parameter file* is safe, in parallel runs
too: unlike a function, a partial is pickled **by value** — it stores a
reference to the wrapped function plus the bound arguments — so the only thing
the workers need to import is ``my_likelihoods`` itself. This is what makes
one likelihood module reusable: to fit four objects with known stellar masses,
write four parameter files binding four different measurements to the same
function.

Forgetting to bind the data is caught early: the bare ``mstar_loglikelihood``
cannot be called with ``(par, state)`` alone, so galapy raises a
:class:`TypeError` when the parameter file is read, before any sampling
starts.

For a likelihood with no per-object data, assign the function itself; fixed
quantities can equally well live as module-level constants or default argument
values in your own module.

What galapy checks for you
--------------------------

When the parameter file is read, the value of ``loglikelihood`` is validated
once in the main process, so that mistakes surface immediately instead of
thousands of likelihood evaluations later:

* a :class:`TypeError` is raised if the value is neither ``None`` nor a
  callable, or if it cannot be called as ``loglikelihood(par, state)``;
* a warning is emitted if it does not accept ``**kwargs``;
* a warning is emitted if it looks unpicklable (lambda, closure, or defined in
  the parameter file itself).

Everything else — returning a scalar, returning ``-inf`` on rejection,
suppressing floating-point errors — is up to you.

Custom likelihoods and model comparison
---------------------------------------

``loglikelihood`` is a **run-wide** setting. It cannot be overridden per entry
of the ``models`` list of a catalogue parameter file, and attempting to do so
raises a :class:`ValueError`.

The reason is that the point of running several model variants on the same
source is to compare their evidences,

.. math::

   Z_k = \int \mathcal{L}_k(\theta)\, \pi_k(\theta)\, \mathrm{d}\theta ,

through the Bayes factor :math:`B_{12} = Z_1 / Z_2`. That ratio says something
about the *models* only if :math:`\mathcal{L}` is the same function of the data
in both runs. The likelihood is the term that carries the data: if one variant
used a Gaussian and another a Student-t, the ratio would largely reflect which
of the two distributions assigns more probability mass to the dataset — an
effect that has nothing to do with the galaxy physics being compared, and one
that is typically much larger than the difference between the models. Since the
result would look perfectly reasonable while being meaningless, galapy refuses
the configuration outright rather than warning about it.

The legitimate uses of a custom likelihood — an external constraint like the
stellar-mass term above (which is one extra *datum*), a heavier-tailed noise
distribution, a different treatment of upper limits — are properties of the
dataset and of your noise assumptions, not of the SED model, and are therefore
meant to apply uniformly to all variants. If you genuinely
want to compare two likelihoods on fixed physics, run them from two separate
parameter files, where nothing invites you to read the ratio as a Bayes factor
about the models.

As a safety net, the results file records which likelihood produced it, as the
``module.qualified_name`` string of the callable used (a
:func:`functools.partial` wrapper is transparent here: the marker identifies
the wrapped function, while the bound values are not recorded):

.. code-block:: python

   >>> res.loglikelihood_name
   'my_likelihoods.mstar_loglikelihood'

and :func:`galapy.analysis.model_comparison.bayes_factor` warns when it is
handed two :class:`~galapy.sampling.Results.Results` objects whose markers
disagree:

.. code-block:: python

   from galapy.analysis.model_comparison import bayes_factor, jeffreys_scale

   log_bf = bayes_factor( res1, res2 )   # warns if the likelihoods differ
   print( jeffreys_scale( log_bf ) )

This catches the cross-parameter-file version of the mistake, which the
restriction on ``models`` cannot.

Performance
-----------

The likelihood is called in the innermost loop of the sampler, often millions
of times per run, and in a typical fit it is dominated by the SED computation
rather than by the metric. Still, avoid doing anything per call that could be
done once: no imports inside the function, no re-reading files, no allocating
large temporaries. Keep the arithmetic vectorised with :mod:`numpy`, and do not
draw random numbers: a stochastic likelihood breaks the assumptions of every
sampler galapy supports.
