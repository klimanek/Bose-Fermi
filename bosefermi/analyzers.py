"""
Pattern analysis and matching for Bose-Fermi integrals.
"""

from sympy import Integral, exp, Symbol, Wild, Pow, Mul, Add, S, logcombine
from typing import Optional, Tuple
import logging

from .patterns import IntegralPattern
from .validators import IntegralValidator
from .special_functions import SymPyExpr


class PatternMatcher:
    """Handles pattern matching for integral expressions."""

    def __init__(self, debug: bool = False):
        self.debug = debug
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """Setup debug logger."""
        logger = logging.getLogger("BoseFermi.PatternMatcher")

        if self.debug and not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter("MATCHER: %(message)s")
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.DEBUG)

        return logger

    def extract_power_and_coefficient(
        self, expr: SymPyExpr, var: "Symbol"
    ) -> Tuple[Optional[SymPyExpr], SymPyExpr]:
        """
        Extract power and coefficient from expressions like c*x**p.

        Returns:
            Tuple of (power, coefficient) or (None, coefficient) if no x**p found
        """
        if expr == var:
            return S.One, S.One
        elif expr.is_number:
            return S.Zero, expr
        elif isinstance(expr, Pow) and expr.base == var:
            return expr.exp, S.One
        elif isinstance(expr, Mul):
            power = None
            coeff_parts = []

            for arg in expr.args:
                if isinstance(arg, Pow) and arg.base == var:
                    power = arg.exp
                elif arg == var:
                    power = S.One
                else:
                    coeff_parts.append(arg)

            coefficient = Mul(*coeff_parts) if coeff_parts else S.One
            return power if power is not None else S.Zero, coefficient
        else:
            # Check if expression contains var at all
            if expr.has(var):
                return None, None  # Complex expression with var that we can't handle
            else:
                return S.Zero, expr  # Constant expression

    def match_exponential_denominator(
        self, denom: SymPyExpr, var: "Symbol"
    ) -> Optional[Tuple[SymPyExpr, SymPyExpr, int]]:
        """
        Match denominator patterns of the form exp(a*x + b) ± 1.

        Returns:
            Tuple of (a, b, sign) or None if no match found
        """
        # Simplify expressions like exp(a)*exp(b) to exp(a+b)
        denom = logcombine(denom, force=True)

        if not isinstance(denom, Add) or len(denom.args) != 2:
            self.logger.debug(f"Denominator is not a sum of two terms: {denom}")
            return None

        exp_term = None
        sign_term = None

        for term in denom.args:
            if term.func == exp:
                exp_term = term
            elif term.is_number and term in [-1, 1]:
                sign_term = term

        if exp_term is None or sign_term is None:
            self.logger.debug(f"Could not find exp(.) ± 1 pattern in {denom}")
            return None

        # Extract and analyze the exponent
        exp_arg = exp_term.args[0]
        self.logger.debug(f"Found exponential with argument: {exp_arg}")

        if not exp_arg.has(var):
            self.logger.debug(
                "Exponential argument does not contain integration variable"
            )
            return None

        # Match exp argument to a*x + b pattern
        a_ = Wild("a", exclude=[var])
        b_ = Wild("b", exclude=[var])

        match = exp_arg.match(a_ * var + b_)
        if not match:
            # Try simpler patterns
            if exp_arg == var:
                a, b = S.One, S.Zero
            elif exp_arg.is_number:
                self.logger.debug("Exponential argument is constant - not suitable")
                return None
            else:
                self.logger.debug(
                    f"Could not match exp argument {exp_arg} to a*x + b pattern"
                )
                return None
        else:
            a = match[a_]
            b = match.get(b_, S.Zero)

        self.logger.debug(
            f"Matched pattern: exp({a}*{var} + {b}) {'+' if sign_term > 0 else '-'} 1"
        )
        return a, b, int(sign_term)


class IntegralAnalyzer:
    """Main analyzer for integral expressions."""

    def __init__(self, debug: bool = False):
        self.debug = debug
        self.validator = IntegralValidator(debug)
        self.matcher = PatternMatcher(debug)
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """Setup debug logger."""
        logger = logging.getLogger("BoseFermi.Analyzer")

        if self.debug and not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter("ANALYZER: %(message)s")
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.DEBUG)

        return logger

    def analyze(self, expr: Integral) -> Integral | IntegralPattern:
        """
        Analyze an integral and return a pattern if recognized.

        Args:
            expr: SymPy Integral to analyze

        Returns:
            IntegralPattern if recognized, None otherwise
        """
        if not self.validator.validate_all(expr):
            return expr

        var = expr.variables[0]
        integrand = expr.function

        self.logger.debug(f"Analyzing integral: {expr}")
        self.logger.debug(f"Integrand: {integrand}")

        # Separate numerator and denominator
        num, denom = integrand.as_numer_denom()

        # Extract power and constant factor from numerator
        power, const_factor = self.matcher.extract_power_and_coefficient(num, var)

        if not self.validator.validate_power_non_negative(power):
            return expr

        self.logger.debug(f"Extracted: power={power}, const_factor={const_factor}")

        # Check if constant factor actually doesn't depend on var
        if const_factor.has(var):
            self.logger.debug("Numerator contains variable in non-power form")
            return expr

        # Match denominator pattern
        pattern_match = self.matcher.match_exponential_denominator(denom, var)
        if pattern_match is None:
            return expr

        a, b, sign = pattern_match

        if not self.validator.validate_coefficient_nonzero(a):
            return expr

        return IntegralPattern(power, a, b, sign, const_factor)
