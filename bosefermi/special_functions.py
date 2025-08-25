"""
Special function utilities for Bose-Fermi integrals.
"""

from typing import Union

from sympy import Expr, Rational, S, exp, gamma, polylog, zeta

SymPyExpr = Union[Expr, Rational, int, float]


def fermi_dirac_special_function(p: SymPyExpr, b: SymPyExpr) -> SymPyExpr:
    """
    Compute special function for Fermi-Dirac integrals.

    For b=0: returns (1 - 1/2^p) * zeta(p+1) [Dirichlet eta function]
    For b≠0: returns polylog(p+1, -exp(-b))
    """
    if b == 0:
        return (S.One - S.One / (2**p)) * zeta(p + 1)
    else:
        return polylog(p + 1, -exp(-b))


def bose_einstein_special_function(p: SymPyExpr, b: SymPyExpr) -> SymPyExpr:
    """
    Compute special function for Bose-Einstein integrals.

    For b=0: returns zeta(p+1)
    For b≠0: returns polylog(p+1, exp(-b))
    """
    if b == 0:
        return zeta(p + 1)
    else:
        return polylog(p + 1, exp(-b))


def gamma_prefactor(p: SymPyExpr, a: SymPyExpr) -> SymPyExpr:
    """Compute gamma function prefactor: gamma(p+1) / a^(p+1)"""
    return gamma(p + 1) / (a ** (p + 1))
