"""
Astrophysical Example: White Dwarf Structure and Chandrasekhar Limit
===================================================================

This example demonstrates how our Bose-Fermi solver helps calculate white dwarf
structure - the final stages of stellar evolution for stars like our Sun.

Physical Background:
- White dwarfs are supported by degenerate electron gas pressure
- Electrons behave as a Fermi gas (Fermi-Dirac statistics)
- Critical mass limit (Chandrasekhar limit) ≈ 1.4 M☉
- Above this limit, the star collapses into a neutron star or black hole

Key Integrals:
- Degenerate gas pressure: ∫ x^(5/2)/(e^x + 1) dx
- Energy density: ∫ x^(3/2)/(e^x + 1) dx
- Thermal corrections at finite temperature
"""

import matplotlib.pyplot as plt
import numpy as np
from sympy import Integral, N, S, exp, gamma, oo, symbols

from bosefermi.core import bose_fermi_integral


def white_dwarf_structure_example():
    """
    Calculate white dwarf structure using Fermi-Dirac integrals
    """
    print("WHITE DWARF STRUCTURE")
    print("=" * 30)

    print("\nPhysical situation:")
    print("After exhausting nuclear fuel, a Sun-like star:")
    print("1. Burns helium in core → carbon-oxygen white dwarf")
    print("2. Outer layers expelled → planetary nebula")
    print("3. Remains hot, dense core ~ Earth size, Sun mass")
    print("4. Only gravitational support: degenerate electron pressure!")

    # Symbol definitions
    x = symbols("x", positive=True)

    print("\nKey quantities:")
    print("x = (E - μ)/(k_B*T) = dimensionless energy")
    print("μ = electron chemical potential")
    print("T = temperature (often T ~ 10^7 K)")
    print("E_F = Fermi energy ~ μ at T=0")

    print("\nDegenerate electron gas:")
    print("Electrons fill states up to Fermi energy E_F")
    print("Pauli: max 2 electrons per state (spin ↑↓)")
    print("→ Enormous pressure even at T=0!")

    # Key integrals for structure
    print("\nKey integrals:")

    # 1. Particle density (electron number)
    print("\n1. Electron density:")
    density_integral = Integral(x ** (S(1) / 2) / (exp(x) + 1), (x, 0, oo))
    print(f"n_e ∝ ∫₀^∞ x^(1/2)/(e^x + 1) dx = {density_integral}")

    density_result = bose_fermi_integral(density_integral)
    print(f"Result: {density_result}")
    print(f"Numerical: {N(density_result, 6)}")

    # 2. Energy density
    print("\n2. Energy density:")
    energy_integral = Integral(x ** (S(3) / 2) / (exp(x) + 1), (x, 0, oo))
    print(f"ε ∝ ∫₀^∞ x^(3/2)/(e^x + 1) dx = {energy_integral}")

    energy_result = bose_fermi_integral(energy_integral)
    print(f"Result: {energy_result}")
    print(f"Numerical: {N(energy_result, 6)}")

    # 3. Pressure (most important for hydrostatic equilibrium!)
    print("\n3. Degenerate gas pressure:")
    pressure_integral = Integral(x ** (S(5) / 2) / (exp(x) + 1), (x, 0, oo))
    print(f"P ∝ ∫₀^∞ x^(5/2)/(e^x + 1) dx = {pressure_integral}")

    pressure_result = bose_fermi_integral(pressure_integral)
    print(f"Result: {pressure_result}")
    print(f"Numerical: {N(pressure_result, 6)}")

    print("\nHydrostatic equilibrium:")
    print("dP/dr = -G*M(r)*ρ(r)/r²")
    print("Gravity ← vs → Degenerate electron pressure")
    print("Our integral determines P(ρ) - the equation of state!")

    return density_result, energy_result, pressure_result


def chandrasekhar_limit_calculation():
    """
    Calculate the Chandrasekhar limit
    """
    print("\nCHANDRASEKHAR LIMIT")
    print("=" * 25)

    print("\nPhysical idea:")
    print("There exists a critical mass M_Ch ≈ 1.4 M☉:")
    print("• M < M_Ch: stable white dwarf")
    print("• M > M_Ch: gravitational collapse → neutron star/black hole")

    print("\nPressure-density relationship:")
    print("For non-relativistic electrons: P ∝ ρ^(5/3)")
    print("For relativistic electrons: P ∝ ρ^(4/3)")

    # Numerical values for estimate
    print("\nNumerical calculation of Chandrasekhar limit:")

    # Constant in Chandrasekhar formula
    # M_Ch = k * (μ_e)^(-2) where μ_e = mean mass per electron
    # For C-O composition: μ_e ≈ 2 (each nucleus contributes Z/A ≈ 0.5 electrons per nucleon)

    mu_e = 2.0  # Mean mass per electron (C-O mixture)

    # Chandrasekhar constant (accurate value)
    chandrasekhar_constant = 5.83  # M☉, accurate theoretical value
    M_ch = chandrasekhar_constant / (mu_e**2)

    print(f"For μ_e = {mu_e} (C-O white dwarf):")
    print(f"M_Ch = {M_ch:.3f} M☉")
    print("Observed value: ~1.4 M☉")

    # Different compositions
    print("\nComposition dependence:")
    compositions = {
        "He": 2.0,  # Helium nuclei
        "C-O": 2.0,  # Carbon-oxygen
        "O-Ne-Mg": 2.15,  # Heavier elements
        "Fe": 2.15,  # Iron nuclei (theoretical)
    }

    for comp, mu_val in compositions.items():
        M_ch_comp = chandrasekhar_constant / (mu_val**2)
        print(f"{comp:8s}: μ_e = {mu_val:.2f}, M_Ch = {M_ch_comp:.3f} M☉")

    return M_ch


def temperature_effects():
    """
    Effects of finite temperature on structure
    """
    print("\nTEMPERATURE EFFECTS")
    print("=" * 25)

    print("\nAt T = 0 (complete degeneracy):")
    print("All states filled up to E_F")
    print("Fermi-Dirac → step function")

    print("\nAt finite T:")
    print("Thermal excitations create particle-hole pairs")
    print("Fermi-Dirac distribution 'smears' Fermi edge")

    x = symbols("x", positive=True)

    # Compare integrals at different temperatures
    print("\nTemperature corrections to pressure:")

    # Expand for small T/E_F
    print("For T << E_F, we can use Sommerfeld expansion:")
    print("P(T) = P(0) * [1 + (π²/12)*(k_B*T/E_F)² + ...]")

    # Our integral for T=0 (step function limit)
    print("\nComparison:")
    pressure_T0 = Integral(x ** (S(5) / 2) / (exp(x) + 1), (x, 0, oo))
    pressure_result = bose_fermi_integral(pressure_T0)

    # For comparison: classical gas would give different result
    classical_result = gamma(S(7) / 2)  # ∫ x^(5/2) dx

    print(f"Quantum degenerate: {N(pressure_result, 6)}")
    print(f"Classical gas: {N(classical_result, 6)}")
    print(f"Ratio (Quantum/Classical): {N(pressure_result / classical_result, 4)}")

    return pressure_result


def plot_white_dwarf_properties():
    """
    Plots of white dwarf properties
    """
    print("\nGRAPHICAL ANALYSIS")
    print("=" * 25)

    # Create data for plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Fermi-Dirac distribution for different temperatures
    x = np.linspace(0, 8, 1000)

    ax1 = axes[0, 0]
    temperatures = [0.1, 0.5, 1.0, 2.0]  # in units of E_F/k_B
    colors = ["blue", "green", "orange", "red"]

    for T, color in zip(temperatures, colors):
        if T < 0.2:  # T=0 approximation
            fermi = np.where(x < 3, 1, 0)  # Step function
            label = "T ≈ 0 (degenerate)"
        else:
            fermi = 1 / (np.exp(x / T) + 1)
            label = f"T = {T} E_F/k_B"
        ax1.plot(x, fermi, color=color, linewidth=2, label=label)

    ax1.set_xlabel("E/E_F")
    ax1.set_ylabel("f(E) - occupation probability")
    ax1.set_title("Fermi-Dirac Distribution")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 2. Integrands for different quantities
    ax2 = axes[0, 1]

    def safe_fermi_integrand(x_vals, power):
        """Safe calculation of integrand x^p/(e^x + 1)"""
        result = np.zeros_like(x_vals)
        nonzero = x_vals > 1e-10
        result[nonzero] = x_vals[nonzero] ** power / (np.exp(x_vals[nonzero]) + 1)
        return result

    density_integrand = safe_fermi_integrand(x, 0.5)  # n_e
    energy_integrand = safe_fermi_integrand(x, 1.5)  # ε
    pressure_integrand = safe_fermi_integrand(x, 2.5)  # P

    ax2.plot(x, density_integrand, "b-", linewidth=2, label="Density: x^(1/2)/(e^x+1)")
    ax2.plot(x, energy_integrand, "g-", linewidth=2, label="Energy: x^(3/2)/(e^x+1)")
    ax2.plot(
        x, pressure_integrand, "r-", linewidth=2, label="Pressure: x^(5/2)/(e^x+1)"
    )

    ax2.set_xlabel("x = E/(k_B*T)")
    ax2.set_ylabel("Integrand")
    ax2.set_title("Integrands for Different Physical Quantities")
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_yscale("log")

    # 3. Mass-radius relation
    ax3 = axes[1, 0]

    # Theoretical M-R relation for white dwarfs: R ∝ M^(-1/3)
    masses = np.linspace(0.2, 1.4, 100)  # M☉
    radii = 5800 * masses ** (-1 / 3)  # km (empirical fit)

    ax3.plot(masses, radii, "purple", linewidth=3, label="R ∝ M^(-1/3)")

    # Mark some known white dwarfs
    known_wd = {
        "Sirius B": (1.02, 5800),
        "Procyon B": (0.60, 7400),
        "40 Eri B": (0.48, 8100),
    }

    for name, (mass, radius) in known_wd.items():
        ax3.plot(mass, radius, "o", markersize=8, label=name)

    # Chandrasekhar limit
    ax3.axvline(
        x=1.4, color="red", linestyle="--", linewidth=2, label="Chandrasekhar Limit"
    )

    ax3.set_xlabel("Mass (M☉)")
    ax3.set_ylabel("Radius (km)")
    ax3.set_title("Mass-Radius Relation for White Dwarfs")
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # 4. Comparison with other compact objects
    ax4 = axes[1, 1]

    # Different star types (corrected neutron star radii)
    star_types = {
        "Main Sequence": ([0.1, 100], [0.1, 100], "blue"),
        "Red Giants": ([0.5, 10], [10, 1000], "red"),
        "White Dwarfs": ([0.2, 1.4], [3000, 20000], "lightgray"),  # km
        "Neutron Stars": ([1.0, 2.5], [10, 15], "gray"),  # km, corrected
    }

    for star_type, (mass_range, radius_range, color) in star_types.items():
        # Fill region
        masses_fill = [mass_range[0], mass_range[1], mass_range[1], mass_range[0]]
        radii_fill = [
            radius_range[0],
            radius_range[0],
            radius_range[1],
            radius_range[1],
        ]

        if star_type in ["White Dwarfs", "Neutron Stars"]:
            # Convert to solar radii for comparison
            radius_range_solar = [r / 695700 for r in radius_range]  # km to R☉
            radii_fill = [
                radius_range_solar[0],
                radius_range_solar[0],
                radius_range_solar[1],
                radius_range_solar[1],
            ]

        ax4.fill(
            masses_fill,
            radii_fill,
            alpha=0.3,
            color=color,
            label=star_type,
            edgecolor="black" if color == "lightgray" else None,
        )

    ax4.set_xlabel("Mass (M☉)")
    ax4.set_ylabel("Radius (R☉)")
    ax4.set_title("White Dwarfs in Context of Other Stars")
    ax4.set_yscale("log")
    ax4.set_xscale("log")
    ax4.grid(True, alpha=0.3)
    ax4.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.suptitle(
        "White Dwarf Astrophysics: Applications of Fermi-Dirac Integrals",
        y=0.97,
        fontsize=14,
        fontweight="bold",
        color="dimgray",
    )
    plt.show()


def stellar_evolution_context():
    """
    Stellar evolution context
    """
    print("\nSTELLAR EVOLUTION AND WHITE DWARFS")
    print("=" * 40)

    print("\nLife cycle of Sun-type stars:")
    print("1. Gravitational collapse → hydrogen ignition (main sequence)")
    print("2. Core H exhaustion → red giant (He burning)")
    print("3. Instability → thermal pulses, planetary nebula")
    print("4. Remnant white dwarf ~ 0.6 M☉, T ~ 100,000 K")
    print("5. Gradual cooling over billions of years")

    print("\nWhy are white dwarfs hot?")
    print("• Gravitational energy → heat during collapse")
    print("• Initial T ~ 100,000 K (cf. Sun ~ 5,800 K)")
    print("• Cool very slowly: t_cool ~ 10^10 years")

    print("\nWhy are they so dense?")
    print("• Sun's mass in Earth's volume!")
    print("• ρ ~ 10^6 g/cm³ (cf. water: 1 g/cm³)")
    print("• 1 teaspoon ~ 5 tons!")

    print("\nRole of quantum mechanics:")
    print("• Classical physics: collapse to a point")
    print("• Quantum mechanics: Pauli exclusion principle")
    print("• Electrons 'resist' further compression")
    print("• Fermi-Dirac pressure = quantum effect!")

    # Energy comparisons
    print("\nCharacteristic energy comparison:")
    print("Kinetic energy: E_kin ~ (ℏ²/2m)(n)^(2/3)")
    print("Gravitational energy: E_grav ~ GM²/R")
    print("Balance determines final R and stability")

    return True


def main():
    """
    Main astrophysical demonstration
    """
    print("ASTROPHYSICAL APPLICATION: WHITE DWARFS")
    print("=" * 50)
    print("Demonstration of Fermi-Dirac integrals in stellar structure")
    print()

    # Basic structure
    density_res, energy_res, pressure_res = white_dwarf_structure_example()

    # Chandrasekhar limit
    M_ch = chandrasekhar_limit_calculation()

    # Temperature effects
    temperature_effects()

    # Graphical analysis
    plot_white_dwarf_properties()

    # Stellar evolution context
    stellar_evolution_context()

    print("\nSUMMARY OF ASTROPHYSICAL APPLICATIONS")
    print("=" * 45)
    print("Successfully applied our Bose-Fermi solver to:")
    print("   • White dwarf structure")
    print("   • Chandrasekhar limit (~1.46 M☉)")
    print("   • Quantum pressure of degenerate electrons")
    print("   • Temperature corrections to T=0 approximation")

    print("\nKey Fermi-Dirac integrals:")
    print(f"   • Density: ∫ x^(1/2)/(e^x+1) dx = {N(density_res, 4)}")
    print(f"   • Energy: ∫ x^(3/2)/(e^x+1) dx = {N(energy_res, 4)}")
    print(f"   • Pressure: ∫ x^(5/2)/(e^x+1) dx = {N(pressure_res, 4)}")

    print("\nPhysical significance:")
    print("   • These integrals determine the equation of state P(ρ)")
    print("   • Equation of state → hydrostatic equilibrium")
    print("   • Equilibrium → white dwarf radius and stability")
    print("   • Critical mass → stellar fate!")

    print("\nObservational connections:")
    print("   • Sirius B: first discovered white dwarf (1862)")
    print("   • Procyon B, 40 Eridani B: other nearby examples")
    print("   • Gaia: catalog of thousands of white dwarfs")
    print("   • Gravitational waves: LISA will detect white dwarfs")

    print("\nFuture applications:")
    print("   • Cosmology: white dwarfs as 'cosmic clocks'")
    print("   • Dark matter: anomalies in white dwarf motion?")
    print("   • Type Ia supernovae: thermonuclear WD explosions")

    return {
        "density_integral": density_res,
        "energy_integral": energy_res,
        "pressure_integral": pressure_res,
        "chandrasekhar_mass": M_ch,
    }


if __name__ == "__main__":
    # Run complete astrophysical demonstration
    results = main()

    print("\nNumerical results:")
    for key, value in results.items():
        if hasattr(value, "evalf"):
            print(f"   {key}: {value} = {N(value, 6)}")
        else:
            print(f"   {key}: {value}")
