import numpy as np

try:
    _trapz = np.trapezoid   # numpy >= 2.0
except AttributeError:
    _trapz = np.trapz


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
