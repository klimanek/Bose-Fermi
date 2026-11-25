"""Logging behaviour (see ROADMAP P1-4).

A ``debug=True`` call used to install a stream handler that stayed attached to
the module logger, so every later ``debug=False`` call kept printing debug
output.
"""

import logging

from sympy import Integral, exp, oo, symbols

from bosefermi import bose_fermi_integral

x = symbols("x", positive=True)
INTEGRAL = Integral(x / (exp(x) - 1), (x, 0, oo))

LOGGERS = (
    "BoseFermi.Analyzer",
    "BoseFermi.PatternMatcher",
    "BoseFermi.Validator",
)


def test_debug_output_is_emitted(capsys):
    bose_fermi_integral(INTEGRAL, debug=True)
    captured = capsys.readouterr()
    assert "ANALYZER:" in captured.err
    assert "MATCHER:" in captured.err


def test_debug_does_not_leak_into_later_calls(capsys):
    bose_fermi_integral(INTEGRAL, debug=True)
    capsys.readouterr()  # drop the debug output of the first call

    result = bose_fermi_integral(INTEGRAL, debug=False)

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""
    assert result != INTEGRAL  # still evaluated correctly
    assert all(
        logging.getLogger(name).level == logging.WARNING for name in LOGGERS
    )
