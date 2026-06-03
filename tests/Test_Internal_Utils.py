
#------------------------------------------------------------------------------#
#           Test of galapy.internal.utils module.
#------------------------------------------------------------------------------#

import math
import tempfile
import os

import numpy as np
import pytest

from galapy.internal.utils import (
    trap_int,
    find_nearest,
    powerlaw_exp_cutoff,
    poly_N,
    unwrap_keys,
    unwrap_values,
    unwrap_items,
    set_nested,
    get_nested,
    func_scalar_or_array,
    quantile_weighted,
    get_credible_interval,
    filter_strings,
    shorten_string,
    FlagVal,
    cat_to_dict,
)


# ---------------------------------------------------------------------------
# trap_int
# ---------------------------------------------------------------------------

class TestTrapInt:

    def test_linear_two_points(self):
        assert trap_int([0., 1.], [0., 1.]) == pytest.approx(0.5)

    def test_linear_three_points(self):
        assert trap_int([0., 0.5, 1.], [0., 0.5, 1.]) == pytest.approx(0.5)

    def test_constant(self):
        assert trap_int([0., 1., 2.], [1., 1., 1.]) == pytest.approx(2.0)

    def test_non_uniform_grid(self):
        # y=x, non-uniform: [0,1,3] -> 0.5*1 + 0.5*(1+3)*2 = 4.5
        assert trap_int([0., 1., 3.], [0., 1., 3.]) == pytest.approx(4.5)

    def test_quadratic_accuracy(self):
        # integral of x^2 from 0 to 1 = 1/3; trapezoid overestimates
        x = np.linspace(0., 1., 10000)
        result = trap_int(x, x**2)
        assert result == pytest.approx(1./3., rel=1e-5)

    def test_numpy_array_input(self):
        x = np.array([0., 1., 2.])
        y = np.array([2., 2., 2.])
        assert trap_int(x, y) == pytest.approx(4.0)


# ---------------------------------------------------------------------------
# find_nearest
# ---------------------------------------------------------------------------

class TestFindNearest:

    def test_scalar_exact(self):
        assert find_nearest([0, 1, 2, 3, 4], 2) == 2

    def test_scalar_below_midpoint(self):
        assert find_nearest([0, 1, 2, 3, 4], 2.3) == 2

    def test_scalar_above_midpoint(self):
        assert find_nearest([0, 1, 2, 3, 4], 2.7) == 3

    def test_scalar_below_range(self):
        assert find_nearest([0, 1, 2, 3, 4], -1.) == 0

    def test_scalar_above_range(self):
        assert find_nearest([0, 1, 2, 3, 4], 10.) == 4

    def test_array_input(self):
        result = find_nearest([10, 20, 30], [12, 27])
        np.testing.assert_array_equal(result, [0, 2])

    def test_returns_scalar_for_scalar_input(self):
        idx = find_nearest([1., 2., 3.], 2.1)
        assert np.ndim(idx) == 0


# ---------------------------------------------------------------------------
# powerlaw_exp_cutoff
# ---------------------------------------------------------------------------

class TestPowerlawExpCutoff:

    def test_scalar_known_value(self):
        # E=1, gamma=1.8, Ecut=300: 1^1.2 * exp(-1/300)
        expected = 0.9966722160545233
        assert powerlaw_exp_cutoff(1.0, 1.8, 300.0) == pytest.approx(expected, rel=1e-6)

    def test_scalar_second_value(self):
        expected = 15.329342132927152
        assert powerlaw_exp_cutoff(10.0, 1.8, 300.0) == pytest.approx(expected, rel=1e-6)

    def test_array_input(self):
        El = np.array([1., 10., 100.])
        expected = np.array([0.99667222, 15.32934213, 179.98452768])
        np.testing.assert_allclose(powerlaw_exp_cutoff(El, 1.8, 300.0), expected, rtol=1e-6)

    def test_decreasing_at_high_energy(self):
        # well above Ecut the exponential suppresses output
        low  = powerlaw_exp_cutoff(1.e1,  1.8, 10.)
        high = powerlaw_exp_cutoff(1.e3,  1.8, 10.)
        assert high < low

    def test_positive_output(self):
        El = np.logspace(-1, 3, 50)
        assert np.all(powerlaw_exp_cutoff(El, 1.8, 300.) > 0)


# ---------------------------------------------------------------------------
# poly_N
# ---------------------------------------------------------------------------

class TestPolyN:

    def test_degree_zero(self):
        assert poly_N(3., [5.]) == pytest.approx(5.0)

    def test_degree_one(self):
        # 2*4 + 3 = 11
        assert poly_N(4., [2., 3.]) == pytest.approx(11.0)

    def test_degree_two(self):
        # x^2 - 2x + 1 at x=3: 9 - 6 + 1 = 4
        assert poly_N(3., [1., -2., 1.]) == pytest.approx(4.0)

    def test_array_input(self):
        # x^2 - 1 at [1, 2, 3]: [0, 3, 8]
        result = poly_N(np.array([1., 2., 3.]), [1., 0., -1.])
        np.testing.assert_allclose(result, [0., 3., 8.], atol=1e-12)

    def test_matches_numpy_polyval(self):
        coeffs = [3., -1., 2., 0.5]
        x = np.linspace(-2., 2., 20)
        np.testing.assert_allclose(poly_N(x, coeffs), np.polyval(coeffs, x), rtol=1e-12)


# ---------------------------------------------------------------------------
# unwrap_keys / unwrap_values / unwrap_items
# ---------------------------------------------------------------------------

_NESTED = {'a': 1, 'b': {'c': 2, 'd': 3}}

class TestUnwrap:

    def test_keys_flat_entry(self):
        keys = list(unwrap_keys(_NESTED))
        assert ['a'] in keys

    def test_keys_nested_entries(self):
        keys = list(unwrap_keys(_NESTED))
        assert ['b', 'c'] in keys
        assert ['b', 'd'] in keys

    def test_keys_count(self):
        assert len(list(unwrap_keys(_NESTED))) == 3

    def test_values(self):
        assert sorted(unwrap_values(_NESTED)) == [1, 2, 3]

    def test_items_keys_and_values(self):
        items = list(unwrap_items(_NESTED))
        keys_from_items   = [k for k, _ in items]
        values_from_items = [v for _, v in items]
        assert sorted(values_from_items) == [1, 2, 3]
        assert ['a'] in keys_from_items
        assert ['b', 'c'] in keys_from_items

    def test_depth_one_dict(self):
        d = {'x': 10, 'y': 20}
        assert sorted(unwrap_values(d)) == [10, 20]


# ---------------------------------------------------------------------------
# set_nested / get_nested
# ---------------------------------------------------------------------------

class TestNestedDictOps:

    def test_set_creates_path(self):
        d = {}
        set_nested(d, ['x', 'y', 'z'], 99)
        assert d == {'x': {'y': {'z': 99}}}

    def test_set_extends_existing(self):
        d = {}
        set_nested(d, ['x', 'y', 'z'], 99)
        set_nested(d, ['x', 'y', 'w'], 100)
        assert d['x']['y'] == {'z': 99, 'w': 100}

    def test_set_overwrites(self):
        d = {}
        set_nested(d, ['a'], 1)
        set_nested(d, ['a'], 2)
        assert d['a'] == 2

    def test_get_nested(self):
        d = {'x': {'y': {'z': 42}}}
        assert get_nested(d, ['x', 'y', 'z']) == 42

    def test_get_nested_depth_one(self):
        d = {'key': 7}
        assert get_nested(d, ['key']) == 7

    def test_roundtrip(self):
        d = {}
        set_nested(d, ['p', 'q'], 3.14)
        assert get_nested(d, ['p', 'q']) == pytest.approx(3.14)

    def test_get_missing_raises(self):
        with pytest.raises(KeyError):
            get_nested({'a': 1}, ['b'])


# ---------------------------------------------------------------------------
# func_scalar_or_array
# ---------------------------------------------------------------------------

class TestFuncScalarOrArray:

    def test_scalar_input_returns_scalar(self):
        result = func_scalar_or_array(4.0, math.sqrt)
        assert np.ndim(result) == 0
        assert result == pytest.approx(2.0)

    def test_array_input_returns_array(self):
        result = func_scalar_or_array(np.array([1., 4., 9.]), math.sqrt)
        np.testing.assert_allclose(result, [1., 2., 3.])

    def test_passes_extra_args(self):
        result = func_scalar_or_array(2., pow, 3)
        assert result == pytest.approx(8.0)

    def test_output_shape_matches_input(self):
        x = np.linspace(1., 4., 10)
        result = func_scalar_or_array(x, math.sqrt)
        assert result.shape == x.shape


# ---------------------------------------------------------------------------
# quantile_weighted
# ---------------------------------------------------------------------------

class TestQuantileWeighted:

    def test_uniform_weights_median_matches_numpy(self):
        vals = np.array([1., 2., 3., 4., 5.])
        result = quantile_weighted(vals, [0.5])
        assert result[0] == pytest.approx(np.percentile(vals, 50), rel=1e-6)

    def test_uniform_weights_quartiles_match_numpy(self):
        vals = np.array([1., 2., 3., 4., 5.])
        result = quantile_weighted(vals, [0.25, 0.5, 0.75], old_style=True)
        expected = np.percentile(vals, [25, 50, 75])
        np.testing.assert_allclose(result, expected, rtol=1e-6)

    def test_invalid_quantile_raises(self):
        with pytest.raises(ValueError):
            quantile_weighted([1., 2., 3.], [-0.1])

    def test_invalid_quantile_above_one_raises(self):
        with pytest.raises(ValueError):
            quantile_weighted([1., 2., 3.], [1.1])

    def test_output_shape(self):
        result = quantile_weighted(np.arange(10, dtype=float), [0.1, 0.5, 0.9])
        assert result.shape == (3,)

    def test_heavier_weight_shifts_median(self):
        vals = np.array([1., 2., 3.])
        # heavy weight on first element: median shifts left
        w_left  = np.array([10., 1., 1.])
        w_right = np.array([1., 1., 10.])
        med_left  = quantile_weighted(vals, [0.5], weights=w_left)[0]
        med_right = quantile_weighted(vals, [0.5], weights=w_right)[0]
        assert med_left < med_right


# ---------------------------------------------------------------------------
# get_credible_interval
# ---------------------------------------------------------------------------

class TestGetCredibleInterval:

    def test_gaussian_68_percent(self):
        rng = np.random.default_rng(42)
        samples = rng.normal(0., 1., 10000)
        idmed = np.argmin(np.abs(samples - np.median(samples)))
        lo, hi = get_credible_interval(samples, idmed, percent=0.68)
        assert lo is not None and hi is not None
        assert lo < 0. < hi
        assert lo == pytest.approx(-1.02, abs=0.05)
        assert hi == pytest.approx( 0.99, abs=0.05)

    def test_upper_limit_when_centre_at_minimum(self):
        samples = np.arange(1., 11.)
        lo, hi = get_credible_interval(samples, 0, percent=0.68)
        assert lo is None
        assert hi is not None

    def test_lower_limit_when_centre_at_maximum(self):
        samples = np.arange(1., 11.)
        lo, hi = get_credible_interval(samples, 9, percent=0.68)
        assert lo is not None
        assert hi is None

    def test_interval_contains_centre(self):
        samples = np.linspace(-5., 5., 1000)
        idmed = len(samples) // 2
        lo, hi = get_credible_interval(samples, idmed, percent=0.68)
        assert lo < samples[idmed] < hi

    def test_wider_percent_gives_wider_interval(self):
        samples = np.linspace(-5., 5., 1000)
        idmed = len(samples) // 2
        lo68, hi68 = get_credible_interval(samples, idmed, percent=0.68)
        lo95, hi95 = get_credible_interval(samples, idmed, percent=0.95)
        assert lo95 <= lo68 and hi95 >= hi68


# ---------------------------------------------------------------------------
# filter_strings
# ---------------------------------------------------------------------------

class TestFilterStrings:

    _LIST = ['sfh.psi', 'sfh.tau', 'ism.f_MC', 'ism.norm']

    def test_wildcard_prefix(self):
        assert filter_strings(self._LIST, 'sfh*') == ['sfh.psi', 'sfh.tau']

    def test_exact_match(self):
        assert filter_strings(self._LIST, 'ism.norm') == ['ism.norm']

    def test_no_match(self):
        assert filter_strings(self._LIST, 'agn*') == []

    def test_list_of_patterns(self):
        result = filter_strings(self._LIST, ['sfh*', 'ism.f*'])
        assert sorted(result) == sorted(['sfh.psi', 'sfh.tau', 'ism.f_MC'])

    def test_invalid_fields_type_raises(self):
        with pytest.raises((ValueError, TypeError)):
            filter_strings(self._LIST, 42)


# ---------------------------------------------------------------------------
# shorten_string
# ---------------------------------------------------------------------------

class TestShortenString:

    def test_short_string_unchanged(self):
        assert shorten_string('short', maxlength=10) == 'short'

    def test_long_string_truncated(self):
        result = shorten_string('galaxy.sfh.psi_max', maxlength=10)
        assert len(result) <= 10 + len('.')
        assert result == 'psi_max'

    def test_single_component(self):
        assert shorten_string('a.b.c.d', maxlength=3) == 'd'

    def test_custom_splitchar(self):
        result = shorten_string('a/b/c', maxlength=3, splitchar='/')
        assert result == 'c'


# ---------------------------------------------------------------------------
# FlagVal
# ---------------------------------------------------------------------------

class TestFlagVal:

    def test_value_attribute(self):
        fv = FlagVal(42, 'my_flag')
        assert fv.value == 42

    def test_flag_attribute(self):
        fv = FlagVal(42, 'my_flag')
        assert fv.flag == 'my_flag'

    def test_repr(self):
        fv = FlagVal(42, 'test_flag')
        assert repr(fv) == 'FlagVal(42,"test_flag")'

    def test_str(self):
        fv = FlagVal(42, 'test_flag')
        assert str(fv) == 'test_flag:42'

    def test_non_integer_value(self):
        fv = FlagVal(3.14, 'pi')
        assert fv.value == pytest.approx(3.14)


# ---------------------------------------------------------------------------
# cat_to_dict
# ---------------------------------------------------------------------------

_CAT_CONTENT = """\
id\tband1\tband1_err\tband2\tband2_err\tredshift
src_a\t1.0\t0.1\t2.0\t0.2\t0.5
src_b\t3.0\t0.3\t4.0\t0.4\t1.0
"""

class TestCatToDict:

    @pytest.fixture
    def cat_file(self, tmp_path):
        f = tmp_path / 'catalog.txt'
        f.write_text(_CAT_CONTENT)
        return str(f)

    def test_returns_dict_with_source_keys(self, cat_file):
        result = cat_to_dict(cat_file, meta_fields=['redshift'])
        assert set(result.keys()) == {'src_a', 'src_b'}

    def test_fluxes_parsed(self, cat_file):
        result = cat_to_dict(cat_file, meta_fields=['redshift'])
        np.testing.assert_allclose(result['src_a']['fluxes'], [1.0, 2.0])

    def test_errors_parsed(self, cat_file):
        result = cat_to_dict(cat_file, meta_fields=['redshift'])
        np.testing.assert_allclose(result['src_a']['errors'], [0.1, 0.2])

    def test_bands_parsed(self, cat_file):
        result = cat_to_dict(cat_file, meta_fields=['redshift'])
        np.testing.assert_array_equal(result['src_a']['bands'], ['band1', 'band2'])

    def test_meta_field_parsed(self, cat_file):
        result = cat_to_dict(cat_file, meta_fields=['redshift'])
        assert result['src_a']['redshift'] == pytest.approx(0.5)
        assert result['src_b']['redshift'] == pytest.approx(1.0)

    def test_second_source(self, cat_file):
        result = cat_to_dict(cat_file, meta_fields=['redshift'])
        np.testing.assert_allclose(result['src_b']['fluxes'], [3.0, 4.0])
