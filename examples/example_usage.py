"""
Example usage of the modular Bose-Fermi integral evaluator.
"""

from sympy import symbols, oo, Integral, exp, pi
from bosefermi.core import bose_fermi_integral


def main():
    x = symbols("x", real=True, positive=True)

    # Famous integral: ∫ x³/(eˣ-1) dx from 0 to ∞ = π⁴/15
    print("Example: Famous Bose-Einstein integral")
    print("=" * 40)

    integral = Integral(x**3 / (exp(x) - 1), (x, 0, oo))
    result = bose_fermi_integral(integral, debug=True)

    print(f"\nIntegral: {integral}")
    print(f"Result: {result}")
    print(f"Numerical value: {float(result):.6f}")
    print(f"π⁴/15 = {float(pi**4 / 15):.6f}")
    print(f"Match: {abs(float(result) - float(pi**4 / 15)) < 1e-10}")


if __name__ == "__main__":
    main()
