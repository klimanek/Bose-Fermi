# Bose-Fermi Integral Evaluator for SymPy

This module provides functionality to automatically recognize and evaluate a class of definite integrals frequently encountered in statistical physics, especially in the context of Bose-Einstein and Fermi-Dirac statistics. *(Will probably be useful to about two people on the planet.)*

These integrals have the general form:

$$
\int_0^\infty \frac{x^p}{e^{\alpha x + \mu} \pm 1} \, \mathrm{d}x
$$

where:
- $p$ is the power of the integration variable $x$,
- $\alpha > 0$ is a scaling factor (inverse temperature),
- $\mu$ is the chemical potential,
- the $-$ sign corresponds to **bosons** (Bose-Einstein),
- the $+$ sign corresponds to **fermions** (Fermi-Dirac).

The module recognizes such integrals symbolically and returns an exact result using known special functions (`gamma`, `polylog`, etc.), whenever possible.

--
## Example

$$I := \int_0^\infty\frac{x^3}{e^x - 1}\,\mathrm{d}x$$

```python
In [1]: from sympy import symbols, Integral, exp, oo
   ...: from bosefermi import bose_fermi_integral
   ...:
   ...: x = symbols("x")
   ...: I = Integral(x**3 / (exp(x) - 1), (x, 0, oo))
   ...: bose_fermi_integral(I)

Out[1]: pi**4/15

```

Result:

$$I = \frac{\pi^4}{15}$$


## Performance
```python
%%timeit
bose_fermi_integral(expr)

1.66 ms ± 6.1 μs per loop (mean ± std. dev. of 7 runs, 1,000 loops each)
```

---

## Physical Interpretation

This integral corresponds, up to a Gamma prefactor, to the *standard* Bose-Einstein or Fermi-Dirac integral defined in physics:

$$
B_k(\mu) = \frac{1}{\Gamma(k+1)} \int_0^\infty \frac{x^k}{e^{x + \mu} - 1} \, \mathrm{d}x \quad \text{(bosons)}
$$

or

$$
F_k(\mu) = \frac{1}{\Gamma(k+1)} \int_0^\infty \frac{x^k}{e^{x + \mu} + 1} \, \mathrm{d}x \quad \text{(fermions)}
$$


---

## Supported Integrals

The following symbolic forms are supported:

| Integral | Result |
|---------|--------|
| $\int_0^\infty \frac{x^p}{e^x - 1} \, dx$ | $\Gamma(p+1)\cdot\zeta(p+1)$ |
| $\int_0^\infty \frac{x^p}{e^{\alpha x} - 1} \, dx$ | $\Gamma(p+1)\cdot\zeta(p+1)/\alpha^{p+1}$ |
| $\int_0^\infty \frac{x^p}{e^{x + \mu} - 1} \, dx$ | $\Gamma(p+1)\cdot\operatorname{Li}_{p+1}(e^{-\mu})$ |
| $\int_0^\infty \frac{x^p}{e^{\alpha x + \mu} - 1} \, dx$ | $\Gamma(p+1)\cdot\operatorname{Li}_{p+1}(e^{-\mu})/\alpha^{p+1}$ |
| $\int_0^\infty \frac{x^p}{e^x + 1} \, dx$ | $\Gamma(p+1)\cdot\eta(p+1)$ |
| $\int_0^\infty \frac{x^p}{e^{\alpha x} + 1} \, dx$ | $\Gamma(p+1)\cdot\eta(p+1)/\alpha^{p+1}$ |
| ... |

---

## Installation

```bash
# For now, clone and use directly:
git clone [your-repo-url]
cd bose-fermi-integrals
python -m pip install -e .
```


## Contributing & Integration

This module is designed for eventual integration into **SymPy**. Contributions are warmly welcome!

### How to contribute:
- **Bug reports**: Found an integral that should work but doesn't? Please report it!
- **New patterns**: Know of related integrals from physics/mathematics? Let's add them!
- **Performance improvements**: Optimizations in pattern matching are always appreciated
- **SymPy integration**: Help with integrating this into the main SymPy codebase

### Future SymPy integration:
We're actively working on integrating this functionality into SymPy's main integration engine. The goal is to make these evaluations available automatically when you call `.doit()` on supported integrals.

**Contact**: [Your GitHub/email] for collaboration opportunities.
