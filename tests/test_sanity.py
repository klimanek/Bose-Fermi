"""Sanity tests for the Bose-Fermi integral evaluator.

These are end-to-end checks through ``bose_fermi_integral``; the closed forms
themselves are cross-checked against quadrature in ``test_numerics.py``.
"""

import time

import pytest
from sympy import S, Integral, exp, gamma, log, oo, pi, simplify, symbols, zeta

from bosefermi import bose_fermi_integral

x = symbols("x", positive=True)


@pytest.mark.parametrize(
    "integral, expected",
    [
        (Integral(x**2 / (exp(x) - 1), (x, 0, oo)), 2 * zeta(3)),
        (Integral(x**3 / (exp(x) - 1), (x, 0, oo)), 6 * zeta(4)),
        (Integral(x**2 / (exp(x) + 1), (x, 0, oo)), (1 - S.One / 4) * 2 * zeta(3)),
        (Integral(x**2 / (exp(2 * x) - 1), (x, 0, oo)), gamma(3) * zeta(3) / 8),
    ],
)
def test_basic_integrals(integral, expected):
    result = bose_fermi_integral(integral)
    assert simplify(result - expected) == 0


def test_symbolic_parameters():
    a, p = symbols("a p", positive=True)
    integral = Integral(x**p / (exp(a * x) - 1), (x, 0, oo))
    expected = gamma(p + 1) * zeta(p + 1) / (a ** (p + 1))
    assert simplify(bose_fermi_integral(integral) - expected) == 0


def test_known_special_values():
    # zeta(2) = pi^2/6  ->  Gamma(2) * zeta(2)
    assert bose_fermi_integral(Integral(x / (exp(x) - 1), (x, 0, oo))) == pi**2 / 6
    # zeta(4) = pi^4/90  ->  6 * zeta(4) = pi^4/15
    assert bose_fermi_integral(Integral(x**3 / (exp(x) - 1), (x, 0, oo))) == pi**4 / 15


def test_edge_cases():
    # p = 0 for bosons diverges (zeta(1)), for fermions it is log(2).
    assert bose_fermi_integral(Integral(1 / (exp(x) - 1), (x, 0, oo))) == S.ComplexInfinity
    assert bose_fermi_integral(Integral(1 / (exp(x) + 1), (x, 0, oo))) == log(2)


@pytest.mark.parametrize(
    "integral",
    [
        Integral(x**2, (x, 1, oo)),  # wrong lower limit
        Integral(x**2, (x, 0, 1)),  # finite limits
        Integral(x**2 / (exp(x**2) - 1), (x, 0, oo)),  # non-linear exponent
        Integral(x**-2 / (exp(x) + 1), (x, 0, oo)),  # negative power
    ],
)
def test_unsupported_returns_original(integral):
    assert bose_fermi_integral(integral) == integral


def test_non_integral_input_is_returned_unchanged():
    expression = x**2
    assert bose_fermi_integral(expression) == expression


def test_consistency():
    integral = Integral(x**2 / (exp(x) - 1), (x, 0, oo))
    assert bose_fermi_integral(integral) == bose_fermi_integral(integral)


def test_performance_sanity():
    integral = Integral(x**2 / (exp(x) - 1), (x, 0, oo))
    start = time.perf_counter()
    bose_fermi_integral(integral)
    assert time.perf_counter() - start < 5.0
