# Patch pro SymPy: Bose–Fermiho integrály

`bosefermi-sympy.patch` je `git format-patch` proti `sympy/sympy@master`
(stav 2026-11-25, soubor `sympy/integrals/meijerint.py`; stejná struktura jako
v SymPy 1.14). Přidává uzavřené tvary pro

$$
\int_0^\infty \frac{x^p}{e^{ax+b}\mp 1}\,\mathrm{d}x
$$

přímo do `integrate()` / `Integral.doit()`.

Kontext a zdůvodnění (proč to nejde přes tabulku Meijerových G funkcí) je
v [`MEIJERINT_TABLE_EXTENSION.md`](MEIJERINT_TABLE_EXTENSION.md).

---

## Co patch obsahuje

| soubor | změna |
|---|---|
| `sympy/integrals/bosefermi.py` | **nový** modul (~130 řádků): matcher `A*exp(a*x+b)+C` a vyhodnocení přes `gamma`, `zeta`/`dirichlet_eta`, `polylog` |
| `sympy/integrals/meijerint.py` | +1 import a +6 řádků háčku v `_meijerint_definite_2()` před G pipeline |
| `sympy/integrals/tests/test_bosefermi.py` | **nový** testovací soubor (9 testů) |

Diffstat: `3 files changed, 213 insertions(+)`.

---

## Jak aplikovat

```bash
cd /cesta/k/sympy
git am /cesta/k/Bose-Fermi/notes/sympy/bosefermi-sympy.patch
```

Nebo bez vytvoření commitu:

```bash
git apply --stat  /cesta/k/Bose-Fermi/notes/sympy/bosefermi-sympy.patch
git apply --check /cesta/k/Bose-Fermi/notes/sympy/bosefermi-sympy.patch
git apply         /cesta/k/Bose-Fermi/notes/sympy/bosefermi-sympy.patch
```

Kdyby se hunk kvůli posunu řádků nechytil, jsou to jen dvě místa:

1. `sympy/integrals/meijerint.py`, importy — za `from .integrals import Integral`:

   ```python
   from .bosefermi import bose_fermi_definite
   ```

2. `_meijerint_definite_2()`, hned za `if f == 0: return S.Zero, True`:

   ```python
   special = bose_fermi_definite(f, x)
   if special is not None:
       _debug('Recognised a Bose-Fermi integral')
       return special
   ```

Zbytek je nový soubor `bosefermi.py` a test.

---

## Co to umí (ověřené výstupy)

```python
>>> x = symbols('x', positive=True)
>>> integrate(x**3/(exp(x) - 1), (x, 0, oo))
pi**4/15
>>> integrate(1/(exp(x) + 1), (x, 0, oo))
log(2)
>>> integrate(x**2/(exp(x + 1) + 1), (x, 0, oo))
-2*polylog(3, -exp(-1))
>>> integrate(x**2/(exp(2*x + 0.3) + 1), (x, 0, oo)).evalf(12)
0.170957109828
>>> a, p = symbols('a p', positive=True)
>>> integrate(x**p/(exp(a*x) - 1), (x, 0, oo))
a**(-p - 1)*gamma(p + 1)*zeta(p + 1)
```

Podporované tvary:

| tvar | výsledek |
|---|---|
| $\int_0^\infty \frac{x^p}{e^{ax+b}-1}dx$ | $\frac{\Gamma(p+1)}{a^{p+1}}\mathrm{Li}_{p+1}(e^{-b})$ |
| $\int_0^\infty \frac{x^p}{e^{ax+b}+1}dx$ | $-\frac{\Gamma(p+1)}{a^{p+1}}\mathrm{Li}_{p+1}(-e^{-b})$ |
| $b=0$ (bosony / fermiony) | $\Gamma(p+1)\zeta(p+1)$ / $\Gamma(p+1)\eta(p+1)$ |
| $p=0$, $b=0$ (fermiony) | $\ln 2$ |
| $A e^{ax+b}+C$ | normalizace na $C\left(e^{ax+b+\log|A/C|}\pm1\right)$ |

Podmínky: $a>0$; bosony navíc $p>0$ a $b\ge 0$; fermiony $p>-1$ (libovolné
reálné $b$). Když podmínky neplatí (nebo tvar nesedí), háček vrátí `None`
a pokračuje stávající integrační pipeline — `Integral(x**3/(e^x-1))` s
`meijerg=False` tedy zůstane nevyřešený.

---

## Ověření

- `git apply --check` proti `meijerint.py` z master prošel.
- Patch byl dočasně zkopírován do instalovaného SymPy 1.14 v `.venv`
  (záloha + obnovení, kontrolní SHA256 sedí) a spuštěn:

  ```
  uv run python -m pytest sympy/integrals/tests/test_bosefermi.py -q
  9 passed
  ```

- Testy pokrývají: BE/FD s `b=0`, `b≠0`, symbolické `p`/`a`, škálované
  jmenovatele, float offset, odmítnutí divergentních případů,
  `meijerg=False` a `conds='separate'`.

---

## Poznámky

- Patch nepřidává žádnou závislost; vše je uvnitř SymPy.
- Pokud `meijerint.py` mezitím změní okolí `_meijerint_definite_2`, stačí
  háček vložit ručně (viz výše) — je to 6 řádků.
- Výsledky pro fermiony s $b=0$ vycházejí jako `dirichlet_eta(p+1)`; SymPy
  je sám převede na $\zeta$ jen po `expand_func()`/`rewrite(zeta)`.
- Soubor je v tomto repu **untracked** stejně jako zbytek `notes/`.
