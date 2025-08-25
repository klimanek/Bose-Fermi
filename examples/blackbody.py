"""
Physical Example: Planck's Law and Blackbody Radiation
=====================================================

This example demonstrates how our Bose-Fermi solver is used to calculate
the total energy radiated by a blackbody (Stefan-Boltzmann Law).

Physical Background:
- A blackbody emits electromagnetic radiation according to Planck's law
- The total energy is obtained by integrating over all frequencies
- The result is a Bose-Einstein integral (photons are bosons!)
"""

import matplotlib.pyplot as plt
import numpy as np
from sympy import Integral, exp, oo, pi, symbols

from bosefermi import bose_fermi_integral


def planck_blackbody_example():
    """
    Planck's Law for Blackbody Radiation and Stefan-Boltzmann Law
    """
    print("PLANCK'S LAW FOR BLACKBODY RADIATION")
    print("=" * 50)

    # Define symbols
    print("\nPhysical Constants and Variables:")
    h, c, k_B, T = symbols(
        "h c k_B T", positive=True
    )  # Planck constant, speed of light, Boltzmann constant, temperature
    nu = symbols("nu", positive=True)  # frequency

    print("h   = Planck constant")
    print("c   = speed of light")
    print("k_B = Boltzmann constant")
    print("T   = blackbody temperature")
    print("ν   = radiation frequency")

    # Planck's law - spectral energy density
    print("\nPlanck's Law:")
    print("Spectral energy density of a blackbody:")

    planck_formula = 8 * pi * h * nu**3 / c**3 * 1 / (exp(h * nu / (k_B * T)) - 1)
    print(f"u(ν,T) = {planck_formula}")

    # Total energy - integration over all frequencies
    print("\nTotal Energy:")
    print("Integrate over all frequencies from 0 to ∞:")

    total_energy_integral = Integral(planck_formula, (nu, 0, oo))
    print(f"U(T) = ∫₀^∞ u(ν,T) dν = {total_energy_integral}")

    # Substitution to convert to standard form
    print("\nSubstitution:")
    print("Introduce dimensionless variable: x = hν/(k_B·T)")
    print("Then: ν = x·k_B·T/h, dν = k_B·T/h dx")

    x = symbols("x", positive=True)

    # After substitution, we get our standard integral
    print("\nAfter substitution, we obtain:")
    standard_integral = Integral(x**3 / (exp(x) - 1), (x, 0, oo))
    print(f"∫₀^∞ x³/(e^x - 1) dx = {standard_integral}")

    # Our solver in action
    print("\nOur Bose-Fermi Solver:")
    result = bose_fermi_integral(standard_integral)
    print(f"Result: {result}")
    print(f"Numerical: {float(result):.6f}")

    # Stefan-Boltzmann Law
    print("\nStefan-Boltzmann Law:")
    stefan_boltzmann_const = 8 * pi**5 * k_B**4 / (15 * h**3 * c**3)

    print("Total energy: U(T) = σ·T⁴")
    print("where σ = 8π⁵k_B⁴/(15h³c³) = Stefan-Boltzmann constant")
    print(f"Stefan-Boltzmann constant: σ = {stefan_boltzmann_const}")

    return standard_integral, result


def plot_planck_distribution():
    """
    Plots the Planck distribution for different temperatures
    """
    print("\nPLANCK DISTRIBUTION PLOT")
    print("=" * 30)

    # Dimensionless variable x = hν/(k_B*T)
    x = np.linspace(0.1, 10, 1000)

    # Planck function in dimensionless form: x³/(e^x - 1)
    def planck_dimensionless(x_input):
        """
        Planck function in dimensionless form: x³/(e^x - 1)
        Handles the limit at x->0, where the value is 0.
        The function is robust for scalar inputs and NumPy arrays.
        """
        # Convert input to NumPy array to ensure array operations
        is_scalar = np.isscalar(x_input)
        x_vals = np.atleast_1d(x_input)  # Ensures x_vals is at least a 1D array

        # Create result array of the same size as x_vals
        result = np.zeros_like(x_vals, dtype=float)

        # Find indices where x_vals is "very close" to zero
        is_near_zero = np.isclose(x_vals, 0, atol=1e-9)

        # Set result to 0 for those indices (limit value)
        result[is_near_zero] = 0.0

        # Compute the original expression for non-zero indices
        non_zero_indices = ~is_near_zero
        result[non_zero_indices] = x_vals[non_zero_indices] ** 3 / (
            np.exp(x_vals[non_zero_indices]) - 1
        )

        # Return scalar if input was scalar
        if is_scalar:
            return result[0]
        else:
            return result

    # Calculate values
    y = planck_dimensionless(x)

    # Create the plot
    plt.figure(figsize=(12, 8))

    # Main plot (dimensionless)
    plt.subplot(2, 2, 1)
    plt.plot(x, y, "b-", linewidth=2, label="$f(x) = \\frac{x^3}{e^x - 1}$")
    plt.xlabel("$x = h\\nu/(k_B T)$")
    plt.ylabel("Spectral density (dimensionless)")
    plt.title("Planck Distribution (Dimensionless)")
    plt.grid(True, alpha=0.3)
    plt.legend()

    # Show the integral - area under the curve
    plt.fill_between(
        x,
        y,
        alpha=0.3,
        color="lightblue",
        label=f"$\int$ = Area = $\\pi^4/15$ ≈ {float(pi**4 / 15):.3f}",
    )
    plt.legend()

    # Plot for different temperatures in real units
    plt.subplot(2, 2, 2)

    # Frequencies in Hz (for visualization)
    nu_hz = np.linspace(1e13, 1e15, 1000)  # from IR to UV

    # Constants (simplified values for illustration)
    h_val = 6.626e-34  # J·s
    k_B_val = 1.381e-23  # J/K

    temperatures = [3000, 4000, 5000, 6000]  # K (stellar temperatures)
    colors = ["red", "orange", "yellow", "cyan"]

    for T_val, color in zip(temperatures, colors):
        x_real = h_val * nu_hz / (k_B_val * T_val)
        y_real = planck_dimensionless(x_real)
        # Scale for better visualization
        y_real_scaled = y_real * T_val**4 / 6000**4  # Relative scaling
        plt.plot(
            nu_hz / 1e14,
            y_real_scaled,
            color=color,
            linewidth=2,
            label=f"T = {T_val} K",
        )

    plt.xlabel("Frequency ($10^{14}$ Hz)")
    plt.ylabel("Spectral density (scaled)")
    plt.title("Planck Distribution for Different Temperatures")
    plt.grid(True, alpha=0.3)
    plt.legend()

    # Integrand detail - function around the maximum
    plt.subplot(2, 2, 3)
    x_detail = np.linspace(0, 6, 500)
    y_detail = planck_dimensionless(x_detail)

    plt.plot(x_detail, y_detail, "g-", linewidth=3)

    # Mark the maximum (Wien's displacement law)
    x_max = 2.82  # approximate maximum
    y_max = planck_dimensionless(x_max)
    plt.plot(x_max, y_max, "ro", markersize=8, label=f"Maximum at x ≈ {x_max:.2f}")

    plt.xlabel("$x = h\\nu/(k_B T)$")
    plt.ylabel("$x^3/(e^x - 1)$")
    plt.title("Function Detail - Wien's Displacement Law")
    plt.grid(True, alpha=0.3)
    plt.legend()

    # Comparison with approximations
    plt.subplot(2, 2, 4)

    # Rayleigh-Jeans approximation (small x): x³/x = x²
    rayleigh_jeans = x**2

    # Wien approximation (large x): x³·e^(-x)
    wien = x**3 * np.exp(-x)

    plt.plot(x, y, "b-", linewidth=2, label="Planck (exact)")
    plt.plot(
        x[x < 2],
        rayleigh_jeans[x < 2],
        "r--",
        linewidth=2,
        label="Rayleigh-Jeans (low ν)",
    )
    plt.plot(x[x > 4], wien[x > 4], "g--", linewidth=2, label="Wien (high ν)")

    plt.xlabel("$x = h\\nu/(k_B T)$")
    plt.ylabel("Spectral density")
    plt.title("Comparison with Classical Approximations")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.yscale("log")

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.suptitle(
        "Planck's Law - Physical Application of the Bose-Einstein Integral",
        y=0.97,
        fontsize=14,
        fontweight="bold",
        color="dimgray",
    )
    plt.show()


def stefan_boltzmann_temperature_dependence():
    """
    Plot of total energy dependence on temperature (T⁴ law)
    """
    print("\nSTEFAN-BOLTZMANN LAW")
    print("=" * 30)

    temperatures = np.linspace(1000, 8000, 100)  # K

    # Total energy ∝ T⁴ (relative units)
    total_energy = temperatures**4

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Linear plot
    ax1.plot(temperatures, total_energy / 1e16, "r-", linewidth=3)
    ax1.set_xlabel("Temperature T (K)")
    ax1.set_ylabel("Total energy U (× $10^{16}$ a.u.)")
    ax1.set_title("Stefan-Boltzmann Law: U ∝ T⁴")
    ax1.grid(True, alpha=0.3)

    # Mark some stellar temperatures
    stellar_temps = {"Red Dwarf": 3000, "Sun": 5778, "Sirius": 9940}
    for star, temp in stellar_temps.items():
        if 1000 <= temp <= 8000:
            energy = temp**4 / 1e16
            ax1.plot(temp, energy, "o", markersize=8, label=star)
    ax1.legend()

    # Log-log plot to verify T⁴ dependence
    ax2.loglog(temperatures, total_energy, "b-", linewidth=3, label="U(T)")

    # Reference line with slope 4
    ref_line = 1e-10 * temperatures**4
    ax2.loglog(temperatures, ref_line, "k--", alpha=0.7, label="∝ T⁴")

    ax2.set_xlabel("Temperature T (K)")
    ax2.set_ylabel("Total energy U (a.u.)")
    ax2.set_title("Verification of T⁴ Dependence (log-log)")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.suptitle(
        "Temperature Dependence of Blackbody Total Energy",
        y=0.98,
        fontsize=14,
        fontweight="bold",
        color="dimgray",
    )
    plt.show()


def physical_applications():
    """
    Discussion of physical applications and significance
    """
    print("\nPHYSICAL APPLICATIONS")
    print("=" * 30)

    print("\nAstrophysical Applications:")
    print("• Stellar luminosity: L = 4πR²σT⁴")
    print("• Stellar classification: temperature from color")
    print("• Cosmic microwave background: T = 2.725 K")
    print("• Planetary energy balance")

    print("\nTechnological Applications:")
    print("• Thermal imaging and infrared cameras")
    print("• Incandescent light bulb efficiency")
    print("• Solar panel design and efficiency")
    print("• Furnace and high-temperature measurements")

    print("\nQuantum Mechanics Significance:")
    print("• First successful quantum theory (Planck, 1900)")
    print("• Introduced concept of energy quantization")
    print("• Led to photon concept and quantum field theory")
    print("• Foundation for modern solid-state physics")

    print("\nMathematical Beauty:")
    print("• Connects quantum statistics to classical thermodynamics")
    print("• Links microscopic (photons) to macroscopic (temperature)")
    print("• Demonstrates power of statistical mechanics")
    print("• Beautiful appearance of fundamental constants")

    return True


def numerical_verification():
    """
    Numerical verification of the integral result
    """
    print("\nNUMERICAL VERIFICATION")
    print("=" * 30)

    # Known exact result
    exact_result = pi**4 / 15
    print(f"Exact result: π⁴/15 = {float(exact_result):.8f}")

    # Numerical integration using numpy
    from scipy.integrate import quad

    def integrand(x):
        if x < 1e-10:
            return 0  # Handle x->0 limit
        return x**3 / (np.expm1(x))  # expm1(x) = e^x - 1, more accurate for small x

    numerical_result, error = quad(
        integrand, 0, 20, limit=100
    )  # Integrate to 20 (effectively ∞)

    print(f"Numerical integration: {numerical_result:.8f}")
    print(f"Integration error estimate: {error:.2e}")
    print(
        f"Relative difference: {abs(numerical_result - float(exact_result)) / float(exact_result) * 100:.6f}%"
    )

    # Convergence analysis
    print("\nConvergence Analysis:")
    upper_limits = [5, 10, 15, 20, 25]

    for upper_limit in upper_limits:
        result, _ = quad(integrand, 0, upper_limit, limit=50)
        relative_error = abs(result - float(exact_result)) / float(exact_result) * 100
        print(
            f"Upper limit {upper_limit:2d}: {result:.6f} (error: {relative_error:.4f}%)"
        )

    return exact_result, numerical_result


def main():
    """
    Main function - complete physical example
    """
    print("PHYSICAL APPLICATION OF BOSE-FERMI INTEGRALS")
    print("=" * 55)
    print("Calculation of blackbody total energy using our solver")
    print()

    # Theoretical calculation
    integral, result = planck_blackbody_example()

    # Numerical verification
    exact_val, numerical_val = numerical_verification()

    # Physical applications
    physical_applications()

    # Visualization
    plot_planck_distribution()
    stefan_boltzmann_temperature_dependence()

    print("\nSUMMARY:")
    print("=" * 20)
    print("Successfully used our Bose-Fermi solver to derive")
    print("the Stefan-Boltzmann law from Planck's law!")
    print(f"Integral ∫₀^∞ x³/(e^x - 1) dx = π⁴/15 ≈ {float(exact_val):.6f}")
    print("This result is fundamental to understanding blackbody radiation")
    print("Applications: astrophysics, thermodynamics, quantum mechanics")

    print("\nPhysical Significance:")
    print("Our integral directly determines:")
    print("• Total energy radiated by stars")
    print("• Stefan-Boltzmann constant")
    print("• Planck's quantum theory of radiation")
    print("• Foundation of quantum statistical mechanics")

    print("\nHistorical Impact:")
    print("• Planck's 1900 derivation launched quantum mechanics")
    print("• Resolved the 'ultraviolet catastrophe' of classical physics")
    print("• Led to Einstein's photon concept (1905)")
    print("• Foundation for modern understanding of light and matter")

    return {
        "integral": integral,
        "exact_result": exact_val,
        "numerical_result": numerical_val,
        "solver_result": result,
    }


if __name__ == "__main__":
    # Run the complete example
    results = main()

    print("\nFinal Results:")
    for key, value in results.items():
        if hasattr(value, "evalf"):
            print(f"   {key}: {value} = {float(value):.8f}")
        else:
            print(f"   {key}: {float(value):.8f}")
