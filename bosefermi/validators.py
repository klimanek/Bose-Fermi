"""
Input validation for Bose-Fermi integrals.
"""

from sympy import Integral, oo
import logging


class IntegralValidator:
    """Validates integral expressions for Bose-Fermi evaluation."""

    def __init__(self, debug: bool = False):
        self.debug = debug
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """Setup debug logger."""
        logger = logging.getLogger("BoseFermi.Validator")
        if self.debug and not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter("VALIDATOR: %(message)s")
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.DEBUG)
        return logger

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
        """Check if power is non-negative (required for gamma function)."""
        if power is None:
            self.logger.debug("Power could not be determined")
            return False

        if isinstance(power, (int, float)) and power < 0:
            self.logger.debug(f"Negative power {power} not supported")
            return False
        return True

    def validate_coefficient_nonzero(self, coefficient) -> bool:
        """Check if exponential coefficient is non-zero."""
        if coefficient == 0:
            self.logger.debug("Exponential coefficient is zero - division by zero")
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
