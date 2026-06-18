
#------------------------------------------------------------------------------#
#           Test of galapy.internal.interp (pure-Python lin_interp class).
#------------------------------------------------------------------------------#

import numpy as np
import pytest

from galapy.internal.interp import lin_interp, log_interp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ac(arr):
    """Return a C-contiguous float64 copy."""
    return np.ascontiguousarray(arr, dtype=np.float64)


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

class TestConstruction:

    def test_constructs_without_error(self):
        lin_interp(_ac([0., 1., 2.]), _ac([0., 1., 2.]))

    def test_get_x_roundtrip(self):
        x = _ac([0., 1., 2., 3.])
        f = lin_interp(x, _ac([0., 1., 4., 9.]))
        np.testing.assert_allclose(f.get_x(), x)

    def test_get_y_roundtrip(self):
        y = _ac([0., 1., 4., 9.])
        f = lin_interp(_ac([0., 1., 2., 3.]), y)
        np.testing.assert_allclose(f.get_y(), y)


# ---------------------------------------------------------------------------
# Interpolation accuracy
# ---------------------------------------------------------------------------

class TestInterpolation:

    def test_exact_on_linear_function(self):
        # y = 2x + 1; interpolation on a linear function must be exact
        xl = _ac([0., 1., 2., 3.])
        yl = _ac([1., 3., 5., 7.])
        f = lin_interp(xl, yl)
        pts = _ac([0.5, 1.5, 2.5])
        np.testing.assert_allclose(f(pts), [2., 4., 6.], rtol=1e-12)

    def test_exact_at_grid_nodes(self):
        x = _ac([0., 1., 2., 3., 4.])
        y = _ac(x**2)
        f = lin_interp(x, y)
        np.testing.assert_allclose(f(x), y, rtol=1e-12)

    def test_quadratic_midpoints_reference(self):
        # y = x^2 on [0,1,2,3,4]; linear interpolation at midpoints
        x = _ac([0., 1., 2., 3., 4.])
        y = _ac(x**2)
        f = lin_interp(x, y)
        pts = _ac([0.5, 1.5, 2.5, 3.5])
        # midpoint of chord (i, i+1) for y=x^2: (i^2+(i+1)^2)/2
        expected = np.array([0.5, 2.5, 6.5, 12.5])
        np.testing.assert_allclose(f(pts), expected, rtol=1e-12)

    def test_accuracy_improves_with_more_points(self):
        # integral of x^2 from 0 to 1 = 1/3; lin_interp uses trapezoid internally
        # finer grid -> smaller error
        x_coarse = _ac(np.linspace(0., 1., 10))
        x_fine   = _ac(np.linspace(0., 1., 1000))
        f_coarse = lin_interp(x_coarse, _ac(x_coarse**2))
        f_fine   = lin_interp(x_fine,   _ac(x_fine**2))
        mid = _ac(np.linspace(0.05, 0.95, 50))
        err_coarse = np.max(np.abs(f_coarse(mid) - mid**2))
        err_fine   = np.max(np.abs(f_fine(mid)   - mid**2))
        assert err_fine < err_coarse


# ---------------------------------------------------------------------------
# Extrapolation
# ---------------------------------------------------------------------------

class TestExtrapolation:

    @pytest.fixture
    def linear_interp(self):
        # y = 2x, defined on [0, 3]
        x = _ac([0., 1., 2., 3.])
        y = _ac([0., 2., 4., 6.])
        return lin_interp(x, y)

    def test_extrapolates_left(self, linear_interp):
        # linear extension: y(-1) = -2
        result = linear_interp(_ac([-1.]))
        np.testing.assert_allclose(result, [-2.], rtol=1e-12)

    def test_extrapolates_right(self, linear_interp):
        # linear extension: y(4) = 8
        result = linear_interp(_ac([4.]))
        np.testing.assert_allclose(result, [8.], rtol=1e-12)


# ---------------------------------------------------------------------------
# integrate method
# ---------------------------------------------------------------------------

class TestIntegrate:

    def test_constant_function(self):
        x = _ac([0., 1., 2., 3., 4.])
        f = lin_interp(x, _ac(np.ones(5)))
        assert f.integrate(0., 4.) == pytest.approx(4.0)

    def test_linear_function(self):
        # integral of x from 0 to 4 = 8
        x = _ac([0., 1., 2., 3., 4.])
        f = lin_interp(x, x)
        assert f.integrate(0., 4.) == pytest.approx(8.0)

    def test_quadratic_trapezoid_approx(self):
        # integral of x^2 from 0 to 4: analytic=64/3≈21.33
        # trapezoid on 5 equidistant points overestimates -> 22.0
        x = _ac([0., 1., 2., 3., 4.])
        f = lin_interp(x, _ac(x**2))
        assert f.integrate(0., 4.) == pytest.approx(22.0, rel=1e-6)

    def test_integrate_converges_to_analytic(self):
        # with many points the trapezoidal error is small
        x = _ac(np.linspace(0., 4., 10000))
        f = lin_interp(x, _ac(x**2))
        assert f.integrate(0., 4.) == pytest.approx(64./3., rel=1e-6)


# ---------------------------------------------------------------------------
# Input variants
# ---------------------------------------------------------------------------

class TestInputVariants:

    def test_non_contiguous_input_accepted(self):
        # lin_interp silently accepts non-contiguous arrays
        x_stride = np.array([0., 0.5, 1., 1.5, 2., 2.5, 3.])[::2]
        y_stride = np.array([0., 0.5, 1., 1.5, 2., 2.5, 3.])[::2]
        assert not x_stride.flags['C_CONTIGUOUS']
        f = lin_interp(x_stride, y_stride)
        np.testing.assert_allclose(f(_ac([1.0])), [1.0], rtol=1e-10)

    def test_two_point_grid(self):
        f = lin_interp(_ac([0., 1.]), _ac([0., 1.]))
        np.testing.assert_allclose(f(_ac([0.5])), [0.5], rtol=1e-12)

    def test_moderate_grid(self):
        # verify correctness on a moderately large grid (10_000 points)
        n = 10_000
        x = _ac(np.linspace(0., 1., n))
        f = lin_interp(x, _ac(np.sin(x)))
        pts = _ac(np.linspace(0.1, 0.9, 100))
        result = f(pts)
        assert result.shape == (100,)
        assert not np.any(np.isnan(result))
        np.testing.assert_allclose(result, np.sin(pts), atol=1e-5)


# ===========================================================================
# log_interp tests
# ===========================================================================

def _lac(arr):
    """Log-spaced C-contiguous float64 array."""
    return np.ascontiguousarray(arr, dtype=np.float64)


class TestLogInterpConstruction:

    def test_constructs_without_error(self):
        log_interp(_lac([1., 2., 4.]), _lac([1., 4., 16.]))

    def test_get_x_roundtrip(self):
        x = _lac([1., 2., 4., 8.])
        f = log_interp(x, _lac(x ** 2))
        np.testing.assert_allclose(f.get_x(), x)

    def test_get_y_roundtrip(self):
        x = _lac([1., 2., 4., 8.])
        y = _lac(x ** 2)
        f = log_interp(x, y)
        np.testing.assert_allclose(f.get_y(), y)

    def test_get_y_preserves_zeros(self):
        x = _lac([1., 2., 4., 8.])
        y = _lac([0., 1., 4., 16.])
        f = log_interp(x, y)
        assert f.get_y()[0] == 0.0


class TestLogInterpInterpolation:

    def test_exact_on_power_law(self):
        # log_interp is exact for any power law f = x^alpha
        x = _lac(np.logspace(0, 3, 10))
        y = _lac(x ** 2.5)
        f = log_interp(x, y)
        pts = _lac(np.logspace(0.1, 2.9, 50))
        np.testing.assert_allclose(f(pts), pts ** 2.5, rtol=1e-12)

    def test_exact_at_grid_nodes(self):
        x = _lac(np.logspace(0, 4, 20))
        y = _lac(x ** 1.7)
        f = log_interp(x, y)
        np.testing.assert_allclose(f(x), y, rtol=1e-12)

    def test_more_accurate_than_lin_interp_for_power_law(self):
        # on a coarse grid, log_interp should beat lin_interp for a power law
        x_coarse = _lac(np.logspace(0, 2, 5))
        y_coarse  = _lac(x_coarse ** 3)
        f_log = log_interp(x_coarse, y_coarse)
        f_lin = lin_interp(x_coarse, y_coarse)
        pts = _lac(np.logspace(0.1, 1.9, 30))
        err_log = np.max(np.abs(f_log(pts) - pts ** 3))
        err_lin = np.max(np.abs(f_lin(pts) - pts ** 3))
        assert err_log < err_lin


class TestLogInterpExtrapolation:

    @pytest.fixture
    def power_law_interp(self):
        # f = x^2 on [1, 8]
        x = _lac([1., 2., 4., 8.])
        return log_interp(x, _lac(x ** 2))

    def test_extrapolates_left_as_power_law(self, power_law_interp):
        # slope in log-log space is 2; f(0.5) = 0.25
        result = power_law_interp(_lac([0.5]))
        np.testing.assert_allclose(result, [0.25], rtol=1e-12)

    def test_extrapolates_right_as_power_law(self, power_law_interp):
        # f(16) = 256
        result = power_law_interp(_lac([16.]))
        np.testing.assert_allclose(result, [256.], rtol=1e-12)

    def test_extrapolation_returns_positive(self, power_law_interp):
        assert power_law_interp(_lac([0.1]))[0] > 0
        assert power_law_interp(_lac([100.]))[0] > 0


class TestLogInterpZeroHandling:

    def test_zero_in_fv_does_not_raise(self):
        x = _lac([1., 2., 4., 8.])
        y = _lac([0., 1., 4., 16.])
        log_interp(x, y)   # must not raise

    def test_interpolation_with_zeros_returns_finite(self):
        x = _lac([1., 2., 4., 8.])
        y = _lac([0., 1., 4., 16.])
        f = log_interp(x, y)
        pts = _lac([1.5, 3., 6.])
        assert np.all(np.isfinite(f(pts)))

    def test_scalar_input(self):
        x = _lac([1., 2., 4., 8.])
        f = log_interp(x, _lac(x ** 2))
        result = f(2.0)
        assert np.ndim(result) == 0
        assert pytest.approx(result, rel=1e-12) == 4.0


class TestLogInterpIntegrate:

    def test_converges_to_analytic_on_fine_grid(self):
        # integral of x^2 from 1 to 4 = 64/3 - 1/3 = 21
        x = _lac(np.logspace(0, np.log10(4), 10_000))
        f = log_interp(x, _lac(x ** 2))
        assert f.integrate(1., 4.) == pytest.approx(21., rel=1e-4)

    def test_finer_grid_more_accurate(self):
        analytic = 21.
        x_c = _lac(np.logspace(0, np.log10(4), 20))
        x_f = _lac(np.logspace(0, np.log10(4), 2000))
        f_c = log_interp(x_c, _lac(x_c ** 2))
        f_f = log_interp(x_f, _lac(x_f ** 2))
        assert abs(f_f.integrate(1., 4.) - analytic) < abs(f_c.integrate(1., 4.) - analytic)

    def test_integrate_over_extrapolated_range_is_finite(self):
        # Verifies integrate works when boundaries lie outside the stored grid.
        # The log-space trapezoid is an approximation, not exact, so only
        # finiteness and positivity are checked here.
        x = _lac([2., 4., 8.])
        f = log_interp(x, _lac(x ** 2))
        result = f.integrate(1., 16.)
        assert np.isfinite(result) and result > 0
