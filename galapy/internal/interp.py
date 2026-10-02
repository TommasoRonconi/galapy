import numpy as np

try:
    _trapz = np.trapezoid   # numpy >= 2.0
except AttributeError:
    _trapz = np.trapz

_FLOOR = 1.e-300   # minimum positive value used before taking logs


class lin_interp:
    """Piecewise-linear interpolator with linear extrapolation beyond the grid.

    Pure-Python replacement for the former ``galapy.internal.interp``
    compiled extension.  Uses ``numpy.interp`` for in-range evaluation and
    propagates the edge slopes for out-of-range points.  The public API is
    identical to the old C++ class so all call sites are unaffected.

    Parameters
    ----------
    xv : array-like
        Strictly increasing x-coordinates of the grid.
    fv : array-like
        Function values at each grid point (same length as *xv*).
    """

    def __init__(self, xv, fv):
        self._xv = np.array(xv, dtype=float)
        self._fv = np.array(fv, dtype=float)
        self._slope_left  = ((self._fv[1]  - self._fv[0]) /
                             (self._xv[1]  - self._xv[0]))
        self._slope_right = ((self._fv[-1] - self._fv[-2]) /
                             (self._xv[-1] - self._xv[-2]))

    # ------------------------------------------------------------------
    # Evaluation

    def __call__(self, x):
        """Evaluate the interpolant at *x* (scalar or array)."""
        scalar = np.ndim(x) == 0
        x = np.atleast_1d(np.asarray(x, dtype=float))
        result = np.interp(x, self._xv, self._fv)
        left  = x < self._xv[0]
        right = x > self._xv[-1]
        result[left]  = (self._fv[0] +
                         self._slope_left  * (x[left]  - self._xv[0]))
        result[right] = (self._fv[-1] +
                         self._slope_right * (x[right] - self._xv[-1]))
        return result[0] if scalar else result

    # ------------------------------------------------------------------
    # Grid accessors

    def get_x(self):
        """Return a copy of the x-grid."""
        return self._xv.copy()

    def get_y(self):
        """Return a copy of the y-grid."""
        return self._fv.copy()

    # ------------------------------------------------------------------
    # Integration

    def integrate(self, aa, bb):
        """Trapezoidal integral from *aa* to *bb*.

        Interior grid points within (aa, bb) are included; the boundary
        values are obtained by evaluating the interpolant (which extrapolates
        linearly when *aa* or *bb* lie outside the grid).
        """
        mask = (self._xv > aa) & (self._xv < bb)
        xs = np.concatenate([[aa], self._xv[mask], [bb]])
        fs = np.concatenate([[self(aa)], self._fv[mask], [self(bb)]])
        return _trapz(fs, xs)


class log_interp:
    """Piecewise power-law interpolator with power-law extrapolation.

    Performs interpolation in log-log space (linear interpolation of
    ``log f`` vs ``log x``), which is exact for power-law functions and
    significantly more accurate than linear interpolation for SEDs that
    span many orders of magnitude.

    Zero or negative values in *fv* are replaced by ``_FLOOR`` (1e-300)
    before taking logarithms; ``get_y()`` always returns the original
    un-floored values.

    Parameters
    ----------
    xv : array-like
        Strictly positive, strictly increasing x-coordinates of the grid.
    fv : array-like
        Function values at each grid point (same length as *xv*).
        May contain zeros; values ≤ 0 are floored for log arithmetic.
    """

    def __init__(self, xv, fv):
        self._xv = np.array(xv, dtype=float)
        self._fv = np.array(fv, dtype=float)
        self._log_xv = np.log(self._xv)
        fv_safe = np.where(self._fv > 0, self._fv, _FLOOR)
        self._log_fv = np.log(fv_safe)
        self._slope_left  = ((self._log_fv[1]  - self._log_fv[0]) /
                             (self._log_xv[1]  - self._log_xv[0]))
        self._slope_right = ((self._log_fv[-1] - self._log_fv[-2]) /
                             (self._log_xv[-1] - self._log_xv[-2]))

    # ------------------------------------------------------------------
    # Evaluation

    def __call__(self, x):
        """Evaluate the interpolant at *x* (scalar or array)."""
        scalar = np.ndim(x) == 0
        x = np.atleast_1d(np.asarray(x, dtype=float))
        log_x = np.log(x)
        result = np.exp(np.interp(log_x, self._log_xv, self._log_fv))
        left  = x < self._xv[0]
        right = x > self._xv[-1]
        result[left]  = np.exp(self._log_fv[0] +
                               self._slope_left  * (log_x[left]  - self._log_xv[0]))
        result[right] = np.exp(self._log_fv[-1] +
                               self._slope_right * (log_x[right] - self._log_xv[-1]))
        return result[0] if scalar else result

    # ------------------------------------------------------------------
    # Grid accessors

    def get_x(self):
        """Return a copy of the x-grid."""
        return self._xv.copy()

    def get_y(self):
        """Return a copy of the original (un-floored) y-grid."""
        return self._fv.copy()

    # ------------------------------------------------------------------
    # Integration

    def integrate(self, aa, bb):
        """Integral from *aa* to *bb* using the log-space trapezoid rule.

        Applies the change of variable u = log(x), giving

            integral ≈ Σ  (f_i x_i + f_{i+1} x_{i+1}) / 2 · log(x_{i+1}/x_i)

        which is the standard trapezoid rule in log(x) space.  This is
        more accurate than the linear-space rule for functions on
        log-spaced grids such as SSP wavelength arrays.

        Boundary values are obtained by evaluating the interpolant
        (power-law extrapolation when *aa* or *bb* lie outside the grid).
        """
        mask = (self._xv > aa) & (self._xv < bb)
        xs = np.concatenate([[aa], self._xv[mask], [bb]])
        fs = np.concatenate([[self(aa)], self._fv[mask], [self(bb)]])
        return 0.5 * np.sum((fs[:-1]*xs[:-1] + fs[1:]*xs[1:]) * np.log(xs[1:]/xs[:-1]))
