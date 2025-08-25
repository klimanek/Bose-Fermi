"""
Core functionality for Bose-Fermi integral evaluation.
"""

from sympy import Integral, simplify
from .analyzers import IntegralAnalyzer
from .special_functions import SymPyExpr


def bose_fermi_integral(expr: Integral, debug: bool = False) -> SymPyExpr:
    """
    Recognizes and evaluates integrals of the form:

        ∫ x^p / (exp(a*x + b) ± 1) dx from 0 to ∞

    These appear in statistical physics (Bose-Einstein and Fermi-Dirac integrals).
    The result is expressed exactly using Gamma, Zeta, Eta or Polylogarithm functions.

    Parameters:
    -----------
    expr : sympy.Integral
        The definite integral to evaluate.
    debug : bool
        Enable debug output.

    Returns:
    --------
    sympy.Expr or sympy.Integral
        The exact result if recognized, else the original integral.

    Examples:
    ---------
    >>> from sympy import symbols, oo, Integral, exp
    >>> x = symbols('x', real=True, positive=True)
    >>> integral = Integral(x**2 / (exp(x) - 1), (x, 0, oo))
    >>> result = bose_fermi_integral(integral)
    >>> print(result)  # 2*zeta(3)
    """
    analyzer = IntegralAnalyzer(debug=debug)

    try:
        pattern = analyzer.analyze(expr)
        if pattern is None:
            if debug:
                analyzer.logger.debug(
                    "No pattern matched - returning original integral"
                )
            return expr

        if debug:
            analyzer.logger.debug(f"Pattern recognized: {pattern}")

        result = pattern.evaluate()

        if debug:
            analyzer.logger.debug(f"Evaluation successful: {result}")

        return simplify(result)

    except Exception as e:
        if debug:
            analyzer.logger.debug(f"Error during evaluation: {e}")
        return expr
