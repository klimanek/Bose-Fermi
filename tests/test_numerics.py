"""Numerical cross-checks of the closed forms against direct quadrature.

These tests compare ``bose_fermi_integral`` with ``mpmath.quad`` for both
statistics, which catches sign and coefficient errors in the special
functions (see ROADMAP P0-1 / P0-2 / P1-3).
"""

import mpmath as mp
import pytest
from sympy import Integral, Rational, exp, log, oo, symbols

from bosefermi import bose_fermi_integral

mp.mp.dps = 25

x = symbols("x", positive=True)

BOSE = -1
FERMI = +1


def quadrature(p, a, b, sign):
    """Numerically integrate x^p / (exp(a*x + b) + sign) on [0, inf)."""
    pp, aa, bb = float(p), float(a), float(b)

    def integrand(t):
        return t**pp / (mp.e ** (aa * t + bb) + sign)

    return mp.quad(integrand, [0, 1, 5, 20, 50, mp.inf])


def assert_matches_quadrature(integral, expected):
    result = bose_fermi_integral(integral)
    assert result != integral, f"integral was not recognized: {integral}"
    numeric = complex(result.evalf(20))
    # The integral is real for real parameters; SymPy's analytic continuation
    # of polylog may leave a purely numerical imaginary residue.
    assert abs(numeric.imag) < 1e-15 * max(1.0, abs(numeric.real))
    value = mp.mpf(numeric.real)
    assert mp.almosteq(value, expected, rel_eps=1e-10, abs_eps=1e-12), (
        f"{integral} -> {result} = {value}, quadrature = {expected}"
    )


@pytest.mark.parametrize("sign", [BOSE, FERMI], ids=["bose", "fermi"])
@pytest.mark.parametrize("b", [0, 1, -1])
@pytest.mark.parametrize("a", [1, 2])
@pytest.mark.parametrize("p", [Rational(1, 2), 1, Rational(3, 2), 2])
def test_matches_quadrature(p, a, b, sign):
    if sign == BOSE and b < 0:
        pytest.skip("Bose-Einstein integral diverges for b < 0")
    integral = Integral(x**p / (exp(a * x + b) + sign), (x, 0, oo))
    assert_matches_quadrature(integral, quadrature(p, a, b, sign))


@pytest.mark.parametrize("p", [Rational(3, 2), Rational(5, 2)])
@pytest.mark.parametrize("b", [Rational(1, 2), Rational(3, 2), 2])
def test_matches_quadrature_noninteger_offset(p, b):
    for sign in (BOSE, FERMI):
        integral = Integral(x**p / (exp(x + b) + sign), (x, 0, oo))
        assert_matches_quadrature(integral, quadrature(p, 1, b, sign))


@pytest.mark.parametrize("b", [0, 1, -1, Rational(1, 2)])
def test_fermi_dirac_p_zero(b):
    integral = Integral(1 / (exp(x + b) + 1), (x, 0, oo))
    assert_matches_quadrature(integral, quadrature(0, 1, b, FERMI))


def test_fermi_dirac_p_zero_is_log_two():
    integral = Integral(1 / (exp(x) + 1), (x, 0, oo))
    assert bose_fermi_integral(integral) == log(2)
