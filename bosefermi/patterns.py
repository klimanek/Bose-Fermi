"""
Pattern representation for recognized Bose-Fermi integrals.
"""

from sympy import S
from .special_functions import (
    fermi_dirac_special_function,
    bose_einstein_special_function,
    gamma_prefactor,
    SymPyExpr,
)


class IntegralPattern:
    """Represents a recognized integral pattern with its solution."""

    # Class constants for statistics types
    BOSE_EINSTEIN = -1
    FERMI_DIRAC = +1

    def __init__(
        self,
        p: SymPyExpr,
        a: SymPyExpr,
        b: SymPyExpr,
        sign: int,
        const_factor: SymPyExpr = S.One,
    ):
        """
        Initialize integral pattern.

        Args:
            p: Power of x in numerator
            a: Coefficient of x in exponential
            b: Constant term in exponential
            sign: -1 for Bose-Einstein, +1 for Fermi-Dirac
            const_factor: Additional constant multiplier
        """
        self.p = p
        self.a = a
        self.b = b
        self.sign = sign
        self.const_factor = const_factor

        # Validate sign
        if sign not in [self.BOSE_EINSTEIN, self.FERMI_DIRAC]:
            raise ValueError(f"Invalid sign {sign}, must be ±1")

    @property
    def statistics_type(self) -> str:
        """Return human-readable statistics type."""
        return "Bose-Einstein" if self.sign == self.BOSE_EINSTEIN else "Fermi-Dirac"

    def evaluate(self) -> SymPyExpr:
        """Evaluate the integral using the appropriate formula."""
        # Compute gamma prefactor
        base_result = gamma_prefactor(self.p, self.a)

        # Compute special function based on statistics
        if self.sign == self.BOSE_EINSTEIN:
            special_func = bose_einstein_special_function(self.p, self.b)
        else:  # FERMI_DIRAC
            special_func = fermi_dirac_special_function(self.p, self.b)

        return self.const_factor * base_result * special_func

    def __repr__(self) -> str:
        return (
            f"IntegralPattern({self.statistics_type}, "
            f"p={self.p}, a={self.a}, b={self.b})"
        )

    def __str__(self) -> str:
        sign_str = "-" if self.sign == self.BOSE_EINSTEIN else "+"
        return (
            f"∫ x^{self.p} / (exp({self.a}*x + {self.b}) {sign_str} 1) dx "
            f"[{self.statistics_type}]"
        )
