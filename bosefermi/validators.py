"""
Input validation for Bose-Fermi integrals.
"""

from sympy import Integral, oo, sympify
import logging

from ._logging import make_logger


class IntegralValidator:
    """Validates integral expressions for Bose-Fermi evaluation."""

    def __init__(self, debug: bool = False):
        self.debug = debug
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """Setup debug logger."""
        return make_logger("BoseFermi.Validator", self.debug, "VALIDATOR")

    def validate_integral_type(self, expr) -> bool:
        """Check if expression is a SymPy Integral."""
        if not isinstance(expr, Integral):
            self.logger.debug("Input is not a SymPy Integral object")
            return False
        return True

    def validate_single_variable(self, expr: Integral) -> bool:
        """Check if integral has exactly one integration variable."""
        if len(expr.variables) != 1:
            self.logger.debug(f"Expected 1 variable, got {len(expr.variables)}")
            return False
        return True

    def validate_integration_limits(self, expr: Integral) -> bool:
        """Check if integration limits are from 0 to infinity."""
        limits = expr.limits[0]
        if len(limits) != 3:
            self.logger.debug(f"Invalid limits format: {limits}")
            return False

        var, lower, upper = limits
        if lower != 0 or upper != oo:
            self.logger.debug(
                f"Limits must be (x, 0, ∞), got ({var}, {lower}, {upper})"
            )
            return False
        return True

    def validate_power_non_negative(self, power) -> bool:
        """Check that the power is not known to be negative.

        A negative power either diverges at x = 0 or is moved into the
        denominator by ``as_numer_denom`` before this check.  Symbolic powers
        whose sign is unknown are accepted.
        """
        if power is None:
            self.logger.debug("Power could not be determined")
            return False

        power = sympify(power)
        if power.is_negative is True:
            self.logger.debug(f"Negative power {power} not supported")
            return False
        return True

    def validate_coefficient_nonzero(self, coefficient) -> bool:
        """Check if exponential coefficient is non-zero."""
        if coefficient == 0:
            self.logger.debug("Exponential coefficient is zero - division by zero")
            return False
        return True

    def validate_coefficient_positive(self, coefficient) -> bool:
        """Check that the exponential coefficient is not known to be <= 0.

        ``a <= 0`` makes the integral diverge: the integrand does not decay
        (or grows) as x -> infinity.  Unknown signs are accepted.
        """
        coefficient = sympify(coefficient)
        if coefficient.is_positive is False:
            self.logger.debug(f"Coefficient {coefficient} is not positive")
            return False
        return True

    def validate_bose_einstein_offset(self, offset) -> bool:
        """Check that a Bose-Einstein offset is not known to be negative.

        For ``b < 0`` the denominator ``exp(a*x + b) - 1`` has a pole at
        ``x = -b/a > 0`` and the integral diverges.  Fermi-Dirac integrals
        stay finite for any real ``b``.
        """
        offset = sympify(offset)
        if offset.is_negative is True:
            self.logger.debug(f"Negative Bose-Einstein offset {offset}")
            return False
        return True

    def validate_all(self, expr) -> bool:
        """Run all validation checks."""
        checks = [
            self.validate_integral_type(expr),
            self.validate_single_variable(expr)
            if isinstance(expr, Integral)
            else False,
            self.validate_integration_limits(expr)
            if isinstance(expr, Integral)
            else False,
        ]
        return all(checks)
