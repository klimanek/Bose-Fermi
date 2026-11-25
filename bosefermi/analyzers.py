"""
Pattern analysis and matching for Bose-Fermi integrals.
"""

from sympy import Integral, exp, Symbol, Wild, Pow, Mul, Add, S, log, logcombine
from typing import Optional, Tuple
import logging

from ._logging import make_logger
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
        return make_logger("BoseFermi.PatternMatcher", self.debug, "MATCHER")

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
    ) -> Optional[Tuple[SymPyExpr, SymPyExpr, int, SymPyExpr]]:
        """
        Match a denominator of the form ``A*exp(a*x + b) + C``.

        The denominator is normalised to the canonical ``exp(a*x + b) +- 1``
        form by pulling ``C`` in front of the integral:

        * ``A/C > 0``: ``A*exp(a*x+b) + C = C*(exp(a*x + b + log(A/C)) + 1)``
        * ``A/C < 0``: ``A*exp(a*x+b) + C = -C*(exp(a*x + b + log(-A/C)) - 1)``

        This also covers the ``c*exp(a*x) + 1`` shapes SymPy produces by
        auto-expanding exponents with float offsets, e.g. ``exp(2*x + 0.3)``.

        Returns:
            Tuple ``(a, b, sign, scale)`` or ``None`` if no match was found.
            ``sign`` is -1 for bosons and +1 for fermions; ``scale`` is the
            constant ``1/C`` (or ``-1/C``) pulled in front of the integral.
        """
        # Simplify expressions like exp(a)*exp(b) to exp(a+b)
        denom = logcombine(denom, force=True)

        if not isinstance(denom, Add) or len(denom.args) != 2:
            self.logger.debug(f"Denominator is not a sum of two terms: {denom}")
            return None

        exp_term = None
        const_term = None

        for term in denom.args:
            _, rest = term.as_coeff_Mul()
            if rest.func == exp:
                if exp_term is not None:
                    self.logger.debug(f"More than one exponential term in {denom}")
                    return None
                exp_term = term
            elif term.is_number and not term.has(var):
                if const_term is not None:
                    self.logger.debug(f"More than one constant term in {denom}")
                    return None
                const_term = term
            else:
                self.logger.debug(f"Unsupported term {term} in denominator {denom}")
                return None

        if exp_term is None or const_term is None:
            self.logger.debug(f"Could not find A*exp(.) + C pattern in {denom}")
            return None

        coefficient, exp_part = exp_term.as_coeff_Mul()
        if not coefficient.is_number or coefficient.is_zero:
            self.logger.debug(f"Non-numeric exponential coefficient in {exp_term}")
            return None

        # Extract and analyze the exponent
        exp_arg = exp_part.args[0]
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

        # Normalise A*exp(a*x + b) + C to exp(a*x + b + shift) +- 1
        ratio = coefficient / const_term
        if ratio.is_positive:
            sign = 1
            shift = log(ratio)
            scale = S.One / const_term
        elif ratio.is_negative:
            sign = -1
            shift = log(-ratio)
            scale = -S.One / const_term
        else:
            self.logger.debug(f"Cannot normalise denominator {denom}")
            return None

        b = b + shift
        self.logger.debug(
            f"Matched pattern: exp({a}*{var} + {b}) "
            f"{'+' if sign > 0 else '-'} 1, scale={scale}"
        )
        return a, b, sign, scale


class IntegralAnalyzer:
    """Main analyzer for integral expressions."""

    def __init__(self, debug: bool = False):
        self.debug = debug
        self.validator = IntegralValidator(debug)
        self.matcher = PatternMatcher(debug)
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """Setup debug logger."""
        return make_logger("BoseFermi.Analyzer", self.debug, "ANALYZER")

    def analyze(self, expr: Integral) -> Optional[IntegralPattern]:
        """
        Analyze an integral and return a pattern if recognized.

        Args:
            expr: SymPy Integral to analyze

        Returns:
            IntegralPattern if recognized, None otherwise.

        ``None`` is only the internal "not recognised" signal.  The public
        :func:`bosefermi.core.bose_fermi_integral` turns it back into the
        original integral, following SymPy's convention that an unevaluated
        ``Integral`` (or any input) comes back unchanged.
        """
        if not self.validator.validate_all(expr):
            return None

        var = expr.variables[0]
        integrand = expr.function

        self.logger.debug(f"Analyzing integral: {expr}")
        self.logger.debug(f"Integrand: {integrand}")

        # Separate numerator and denominator
        num, denom = integrand.as_numer_denom()

        # Extract power and constant factor from numerator
        power, const_factor = self.matcher.extract_power_and_coefficient(num, var)

        if not self.validator.validate_power_non_negative(power):
            return None

        self.logger.debug(f"Extracted: power={power}, const_factor={const_factor}")

        # Check if constant factor actually doesn't depend on var
        if const_factor.has(var):
            self.logger.debug("Numerator contains variable in non-power form")
            return None

        # Match denominator pattern
        pattern_match = self.matcher.match_exponential_denominator(denom, var)
        if pattern_match is None:
            return None

        a, b, sign, scale = pattern_match

        if not self.validator.validate_coefficient_nonzero(a):
            return None
        if not self.validator.validate_coefficient_positive(a):
            return None
        if (
            sign == IntegralPattern.BOSE_EINSTEIN
            and not self.validator.validate_bose_einstein_offset(b)
        ):
            return None

        return IntegralPattern(power, a, b, sign, scale * const_factor)
