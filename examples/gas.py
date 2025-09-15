"""
Physics Example: Energy of a Free Fermi Gas (Electrons in Metals)
==================================================================

This example demonstrates how our Bose-Fermi solver calculates the total energy
of free electrons in metals at finite temperature - a fundamental problem in
solid state physics and condensed matter theory.

Physical background:
- Electrons in metals are often modeled as a free Fermi gas (a fundamental first approximation).
- At T=0: all states filled up to Fermi energy (Fermi sea).
- At T>0: thermal excitations create particle-hole pairs in a narrow energy band around E_F.
- Total energy and other thermodynamic properties require Fermi-Dirac integrals that our solver handles exactly!
"""

from sympy import Integral, N, S, exp, gamma, oo, plot, sqrt, symbols

from bosefermi.core import bose_fermi_integral


def fermi_gas_energy_example():
    """
    Calculate the total energy of a free electron gas in metals
    """
    print("FREE ELECTRON GAS IN METALS")
    print("=" * 40)

    # Physical setup
    print("\nPhysical Setup:")
    print("Consider free electrons in a metal (e.g., copper, aluminum)")
    print("- Electrons treated as free particles in a box (simplification)")
    print("- Pauli exclusion principle → Fermi-Dirac statistics")
    print(
        "- At finite temperature T, only electrons near E_F are thermally excited, creating particle-hole pairs."
    )

    # Define symbols
    E, k_B, T, E_F = symbols("E k_B T E_F", positive=True)
    x = symbols("x", positive=True)

    print("\nVariables:")
    print("E   = electron energy")
    print("E_F = Fermi energy (chemical potential at T=0)")
    print("k_B = Boltzmann constant")
    print("T   = temperature")

    # Fermi-Dirac distribution
    print("\nFermi-Dirac Distribution:")
    fermi_dirac = 1 / (exp((E - E_F) / (k_B * T)) + 1)
    print(f"f(E) = {fermi_dirac}")
    print("This gives the probability that a state with energy E is occupied.")

    # Energy density (simplified 3D free electron model)
    print("\nEnergy Calculation:")
    print("Total energy per unit volume (U) involves integrating:")
    print("U = ∫₀^∞ E · g(E) · f(E) dE")
    print("where g(E) ∝ √E is the density of states for a 3D free electron gas.")

    # For demonstration, consider the integral after appropriate substitutions
    print("\nDimensionless Variables:")
    print(
        "Let x = E/(k_B*T). When considering integrals for properties like specific heat (in the Sommerfeld expansion),"
    )
    print("terms proportional to these Fermi-Dirac integrals appear:")

    # Key integral that appears in the calculation (relevant for energy at finite T)
    key_integral = Integral(x ** (S(3) / 2) / (exp(x) + 1), (x, 0, oo))
    print(f"∫₀^∞ x^(3/2)/(e^x + 1) dx = {key_integral}")

    print("\nOur Bose-Fermi Solver in Action:")
    result = bose_fermi_integral(key_integral, debug=True)

    print("\nExact Result:")
    print(f"∫₀^∞ x^(3/2)/(e^x + 1) dx = {result}")
    print(f"Numerical value: {N(result, 6)}")

    # Physical interpretation
    print("\nTemperature Effects:")
    print(
        "This integral (and others of its type) captures how thermal excitations modify the electron energy and other properties."
    )
    print("- At T=0: electrons fill states up to E_F (step function-like behavior).")
    print(
        "- At T>0: thermal broadening creates excited electrons and holes, primarily near E_F."
    )
    print(
        "- These integrals are fundamental building blocks for computing temperature-dependent corrections to ground state properties."
    )

    return key_integral, result


def compare_fermi_integrals():
    """
    Compare different Fermi-Dirac integrals that appear in physics
    """
    print("\nCOMPARISON OF FERMI-DIRAC INTEGRALS")
    print("=" * 45)

    x = symbols("x", positive=True)

    # Different powers that appear in physics
    integrals_and_physics = [
        (
            Integral(sqrt(x) / (exp(x) + 1), (x, 0, oo)),
            "x^(1/2)",
            "Related to particle density",
        ),
        (
            Integral(x ** (S(3) / 2) / (exp(x) + 1), (x, 0, oo)),
            "x^(3/2)",
            "Related to total energy",
        ),
        (
            Integral(x**2 / (exp(x) + 1), (x, 0, oo)),
            "x^2",
            "Related to transport phenomena (e.g., electrical conductivity)",
        ),
        (
            Integral(x**3 / (exp(x) + 1), (x, 0, oo)),
            "x^3",
            "Related to transport phenomena (e.g., thermal conductivity)",
        ),
    ]

    results = []

    for integral, power_str, physics_meaning in integrals_and_physics:
        print(f"\n∫₀^∞ {power_str}/(e^x + 1) dx:")
        print(f"  Physical meaning: {physics_meaning}")

        result = bose_fermi_integral(integral)
        results.append(float(N(result)))

        print(f"  Exact result: {result}")
        print(f"  Numerical: {N(result, 6)}")

    return results


def plot_fermi_dirac_integrands():
    """
    Plot the integrands to visualize the physics
    """
    print("\nVISUALIZATION OF FERMI-DIRAC INTEGRANDS")
    print("=" * 48)

    x = symbols("x", positive=True)

    # Define the integrands
    integrand_1_2 = sqrt(x) / (exp(x) + 1)  # x^(1/2)/(e^x + 1)
    integrand_3_2 = x ** (S(3) / 2) / (exp(x) + 1)  # x^(3/2)/(e^x + 1)
    integrand_2 = x**2 / (exp(x) + 1)  # x^2/(e^x + 1)

    print("Creating plots of the integrands...")
    print(
        "These show how different physical quantities, when integrated, are weighted across the energy spectrum."
    )

    # Create plots
    p1 = plot(
        integrand_1_2,
        (x, 0, 8),
        title="Fermi-Dirac Integrands",
        xlabel=r"$x = E/(k_BT)$",
        ylabel="Integrand value",
        label="$x^{1/2}/(e^x + 1)$",
        legend=True,
        show=False,
        line_color="blue",
    )

    p2 = plot(
        integrand_3_2,
        (x, 0, 8),
        label="$x^{3/2}/(e^x + 1)$",
        show=False,
        line_color="red",
    )

    p3 = plot(
        integrand_2, (x, 0, 8), label="$x^2/(e^x + 1)$", show=False, line_color="green"
    )

    # Combine plots
    p1.extend(p2)
    p1.extend(p3)
    p1.show()

    print("\nKey Features of the Plots:")
    print(
        "- All functions peak at low x, reflecting contributions from lower energy states."
    )
    print(
        "- Exponential cutoff at high energies (x >> 1) due to Fermi-Dirac statistics."
    )
    print("- Different powers give different weights to high/low energies.")
    print("- The area under each curve corresponds to the computed integral values.")


def temperature_scaling_demo():
    """
    Show how results scale with temperature
    """
    print("\nTEMPERATURE SCALING DEMONSTRATION")
    print("=" * 42)

    x = symbols("x", positive=True)

    # The key integral for energy
    energy_integral = Integral(x ** (S(3) / 2) / (exp(x) + 1), (x, 0, oo))
    exact_result = bose_fermi_integral(energy_integral)

    print(f"For the energy integral ∫₀^∞ x^(3/2)/(e^x + 1) dx = {exact_result}")
    print(f"Numerical value: {N(exact_result, 8)}")

    print("\nPhysical Significance:")
    print(
        "In the Sommerfeld expansion, these integrals are part of the coefficients for temperature-dependent corrections."
    )
    print(
        "For instance, the electronic heat capacity (C_el) of metals at low temperatures is proportional to T:"
    )
    print("C_el ∝ T · (coefficient involving π²/6 etc.)")
    print("")
    print(
        "This linear-in-T electronic heat capacity is a direct consequence of Fermi-Dirac statistics,"
    )
    print("as only electrons near the Fermi surface can be thermally excited.")
    print("(This contrasts with the T³ dependence for phonons in the Debye model).")

    # Show numerical values for comparison
    print("\nComparison with Classical (Boltzmann) Result:")
    # For a free classical gas, the integral related to total energy would be:
    classical_integral_value = gamma(S(3) / 2 + 1)  # Gamma(5/2)

    print(
        f"Integral for a classical gas (no Pauli exclusion): {classical_integral_value} = {N(classical_integral_value, 6)}"
    )
    print(
        f"Integral for a quantum Fermi gas (with Pauli exclusion): {exact_result} = {N(exact_result, 6)}"
    )
    print(f"Ratio (Quantum/Classical): {N(exact_result / classical_integral_value, 4)}")
    print(
        "The ratio shows how quantum effects significantly reduce the contribution at this energy scale compared to classical predictions."
    )

    return exact_result


def main():
    """
    Complete physics example demonstrating our solver
    """
    print("PHYSICS APPLICATION: FREE ELECTRON GAS")
    print("=" * 50)
    print("Demonstrating exact evaluation of Fermi-Dirac integrals")
    print("that appear in solid state physics calculations\n")

    # Main energy calculation
    integral, result = fermi_gas_energy_example()

    # Visualize the integrands
    plot_fermi_dirac_integrands()

    print("\nSUMMARY")
    print("=" * 20)
    print("Successfully used our Bose-Fermi solver for solid state physics!")
    print("Computed exact values for key Fermi-Dirac integrals.")
    print(
        "Showed connection to measurable quantities like heat capacity and transport phenomena."
    )
    print(
        "Demonstrated advantages over purely numerical integration for these specific forms."
    )

    print("\nWhy This Matters:")
    print("• These integrals appear in every solid state physics textbook.")
    print(
        "• They are often approximated numerically or with series expansions (e.g., Sommerfeld expansion)."
    )
    print(
        "• Our solver can provide EXACT symbolic results for these specific forms, which is powerful!"
    )
    print(
        "• Applications: electronic properties of metals, semiconductors, thermoelectrics."
    )

    print("\nExperimental Connection:")
    print("The computed integral values directly relate to:")
    print("• Electronic heat capacity measurements in metals.")
    print("• Thermal conductivity of electrons.")
    print("• Thermoelectric effects (e.g., Seebeck coefficient).")

    return integral, result


if __name__ == "__main__":
    # Run the complete example
    main()
