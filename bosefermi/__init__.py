"""
Bose-Fermi Integral Evaluator for SymPy

A modular system for recognizing and evaluating integrals that appear
in statistical physics, particularly Bose-Einstein and Fermi-Dirac statistics.
"""

from .core import bose_fermi_integral
from .patterns import IntegralPattern
from .analyzers import IntegralAnalyzer
from .validators import IntegralValidator

__version__ = "0.1.0"
__author__ = "Oldřich Klimánek"
__email__ = "oldrich.klimanek@gmail.com"

__all__ = [
    "bose_fermi_integral",
    "IntegralPattern",
    "IntegralAnalyzer",
    "IntegralValidator",
]
