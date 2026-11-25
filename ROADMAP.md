# ROADMAP — Bose-Fermi Integral Evaluator

Stav revize: 2026-11-25, větev `main`, HEAD `bf5e9e9` (po sérii oprav dle tohoto dokumentu).
`uv run pytest -q` → 100 passed, 8 skipped (numerické kros-checky proti `mpmath.quad`).

Legenda priorit: **P0** = špatný výsledek, **P1** = testy/ kvalita, **P2** = robustnost, **P3** = balení.

---

## P0 — Chyby s dopadem na správnost výsledku

### P0-1: Chybí znaménko u Fermi–Dirac pro `b ≠ 0`

**Kde:** `bosefermi/special_functions.py:22`

```python
# nyní
return polylog(p + 1, -exp(-b))
# správně
return -polylog(p + 1, -exp(-b))
```

**Proč:** $\int_0^\infty \frac{x^p}{e^{ax+b}+1}dx
= -\frac{\Gamma(p+1)}{a^{p+1}}\mathrm{Li}_{p+1}(-e^{-b})$.
Odvození: $\frac{1}{e^{ax+b}+1}=\sum_{n\ge1}(-1)^{n-1}e^{-n(ax+b)}$, takže součet je $-Li_{p+1}(-e^{-b})$.

**Repro** (mpmath quad vs. sympy):

| integrál | `bose_fermi_integral` | numericky |
|---|---|---|
| `Integral(x**2/(exp(x+1)+1), (x,0,oo))` | −0.7051297586 | +0.7051297586 |
| `Integral(x**3/(exp(3*x-1)+1), (x,0,oo))` | −1.1161851025 | +1.1161851025 |
| `Integral(x**2/(exp(x)+1), (x,0,oo))` (b=0, větev η) | +0.2253856693 | +0.2253856693 ✔ |

Týká se **všech** FD integrálů s chemickým potenciálem (μ ≠ 0). BE větev je v pořádku.
Pozn.: `README.md` ani tabulka v `notes/PROJECT_ANALYSIS.md` variantu FD s `b ≠ 0` neuvádějí — chybí i dokumentace.

- [x] Opravit znaménko
- [x] Doplnit řádek do tabulky podporovaných integrálů v `README.md`
- [x] Regresní test proti `mpmath.quad` (viz P1-3)

### P0-2: FD s `p = 0` vrací `nan` místo `ln 2`

**Kde:** `bosefermi/special_functions.py:20`

```python
return (S.One - S.One / (2**p)) * zeta(p + 1)
```

Pro `p = 0` je to `(1 − 1)·zeta(1)` = `0*zoo` = `nan`, správně má být `log(2)`.
(BE `p = 0` vrací `zoo`, což je pro divergentní integrál obhajitelné.)

- [x] Ošetřit `p == 0` → `log(2)` (nebo `limit`), případně obecně `eta(p+1)` přes `dirichlet_eta`
- [x] Test: `Integral(1/(exp(x)+1), (x,0,oo)) == log(2)`

---

## P1 — Testy a kvalita kódu

### P1-1: `tests/test_sanity.py` nikdy neselže

7 funkcí vrací `bool` místo `assert`; pytest návratovou hodnotu ignoruje
(`PytestReturnNotNoneWarning`). Tím prošel i P0-1.

- [x] Přepsat na `assert` (nebo na `pytest.mark.parametrize`)
- [x] Odstranit `debug=True` z testů, kde jen zanáší výstup

### P1-2: Podmíněné asserty v `tests/test_bosefermi.py`

Místa `if isinstance(result, IntegralPattern): assert ...` (ř. 326, 337, 350, 500)
tiše přeskočí ověření, když match selže.

- [x] Nahradit `assert isinstance(result, IntegralPattern)` + následné asserty
- [x] Přidat test, že `Integral(x**-2/(exp(x)+1))` vrací původní integrál (dnes ano, ale jen
      shodou okolností přes `as_numer_denom` — viz P1-5)

### P1-3: Chybí numerický kros-check

Všechny testy jsou strukturální; žádný neporovnává výsledek s nezávislým výpočtem.

- [x] Parametrizovaný test: symbolický výsledek vs. `mpmath.quad` pro `p ∈ {0, 1, 3/2, 5/2}`,
      `a ∈ {1, 2, 3}`, `b ∈ {−1, 0, 0.3, 1, 2}`, oba `sign` → `tests/test_numerics.py`
- [x] Tolerovat jen relativní odchylku ~1e-10

### P1-4: Debug logger prosakuje mezi voláními

**Kde:** `bosefermi/analyzers.py:21–32`, `145–156`

`_setup_logger` nikdy neodebere handler a `self.logger.debug(...)` se volá bezpodmínečně,
takže po jednom `debug=True` tiskne debug i `debug=False`:

```python
bose_fermi_integral(I1, debug=True)
bose_fermi_integral(I2, debug=False)   # → stále "ANALYZER: ..." / "MATCHER: ..."
```

- [x] Handler na instanci (`logger.addHandler` v `__init__`) + `propagate = False`, nebo
      `logger.disabled`/`NullHandler`, když `debug=False`
- [x] Zvážit jedno sdílené logger nastavení místo tří identických `_setup_logger`

### P1-5: `validate_power_non_negative` je mrtvá kontrola

**Kde:** `bosefermi/validators.py:62`

`isinstance(power, (int, float))` je pro `S(-2)` `False` (SymPy `Integer` není `int`):
`validate_power_non_negative(S(-2))` → `True`, totéž pro symbol s `negative=True`.
Dnes to nekousne jen proto, že `as_numer_denom` přesune zápornou mocninu do jmenovatele.

- [x] `if power.is_negative is True:` / `elif power.is_nonnegative is False:`
- [x] Testy na `S(-2)`, `Rational(-1,2)`, `Symbol('n', negative=True)`

### P1-6: `except Exception` v `core.py:60`

Chyby při evaluaci se tiše spolknou a vrátí se původní integrál.

- [x] Nechat fallback, ale logovat vždy (i bez `debug`), nebo použít `warnings.warn`
- [ ] Zúžit na očekávané výjimky, aby se neskryl bug v pattern matchingu
      — záměrně zůstává široký `except`, ale každá výjimka je nyní vidět jako `RuntimeWarning`

---

## P2 — Robustnost a rozsah

### P2-1: Float v exponentu rozbije matchování

`exp(2*x + 0.3)` se v SymPy automaticky rozvine na `1.3498...*exp(2*x)`,
takže `term.func == exp` v `match_exponential_denominator` selže a integrál zůstane nevyřešen
(u celočíselného `b` to funguje).

- [x] Před matchováním zkusit `nsimplify`/`rationalize` exponentu, nebo rozpoznat `c*exp(a*x)` tvar
      (vyřešeno obecným `A*exp(a*x+b) + C`)
- [x] Test s `b = 0.3`

### P2-2: Nenormalizované tvary jmenovatele

Aktuálně nerozpoznané (vrací původní integrál — bezpečné, ale zbytečné):

| vstup | poznámka |
|---|---|
| `x/(2*(exp(x)−1))` | vytkne se na `2*exp(x) − 2`, `-2 ∉ {±1}` |
| `x/(exp(x)/2 + 1)` | koeficient u `exp` |
| `x/(exp(x+1) − E)` | dva `Mul` členy v `Add` |
| `x**3/(1 − exp(−x))` | ekvivalent `x**3*exp(x)/(exp(x)−1)` |
| `x**3*exp(−x)/(1 − exp(−x))` | numerátor s `exp` |

- [x] Přenést konstantní faktor z celého jmenovatele do `const_factor`
      (obecné `A*exp(a*x+b) + C`, viz P2-1)
- [x] Rozšířit testy o tyto tvary (4 z 5; `x**3/(1−exp(−x))` v `x → ∞` diverguje, zůstává nerozpoznaný)
- [ ] Rewriting numerátoru s `exp(-x)` před analýzou
      (např. `x**3*exp(-x)/(1−exp(−x))` = `x**3/(e^x−1)`) — zbývá

### P2-3: Chybí kontroly konvergence

Dokumentace předpokládá `a > 0`, `b > 0` a `Re p > −1`, ale validátor kontroluje jen `a != 0`
a numericky záporné `p`. Symbolické `a` s neznámým znaménkem projde a vrátí formálně neplatný výsledek.
Negativní `a` je zamítnuto jen mimochodem (přes přepis `exp(−x) → 1/exp(x)`).

- [x] `a.is_positive` — zamítá se jen známé `False` (záporné nebo nulové `a`), neznámé znaménko projde
- [x] `b > 0` pro `polylog` větev: BE s `b < 0` má pól na integrační cestě → zamítnuto;
      FD funguje pro libovolné reálné `b` (analytické pokračování ověřeno proti kvadratuře)
- [x] `Re(p) > −1` pro konvergenci (záporné `p` zamítnuto přes `is_negative`)

---

## P3 — Balení, CI, úklid

- [x] Přesunout `pytest` z `[project].dependencies` do `[dependency-groups].dev`
- [x] `matplotlib`, `numpy`, `scipy` do vlastní skupiny `examples` (+ `[tool.uv] default-groups`),
      runtime zůstává jen `sympy` (knihovna sama potřebuje jen `sympy`)
- [x] `.venv/` bylo prázdné (jen `_virtualenv`, chybí i `sympy`) → vyřešeno: `uv sync` po opravě
      setuptools discovery (viz níže) prostředí znovu vytvořil a nainstaloval 135 balíčků
- [x] Build padal na `Multiple top-level packages discovered in a flat-layout: ['notes', 'bosefermi']`
      → `pyproject.toml` má nově `[tool.setuptools.packages.find] include = ["bosefermi*"]`;
      `uv build` nyní obsahuje jen `bosefermi/`, `notes/` se nebalí
- [x] `bosefermi.egg-info/` smazán z pracovního stromu
- [x] `uv.lock` commitnut (odpovídá novému rozdělení závislostí); `notes/` ponecháno
      autorovi jako untracked; zastaralý `white_dwarf.md` (zachycený stdout staré verze
      `examples/white_dwarf.py`) smazán
- [x] CI: GitHub Actions `.github/workflows/tests.yml` — pytest na 3.10–3.13 + build job
      (osa SymPy 1.12/1.14 zatím není)
- [ ] `notes/PROJECT_ANALYSIS.md` aktualizovat (soubor není commitnutý, necháno na autorovi)

---

## Tipy k architektuře

1. **Jedno místo pro formule.** `special_functions.py` je dobrý nápad, ale rozpadá se do `if b == 0 / else`.
   Zavést jednu funkci `special_function(p, b, sign)` s jasně oddělenými případy a unit testy
   proti známým hodnotám (`Li_s(−1) = −η(s)`, `Li_s(1) = ζ(s)`).
2. **Vracet `IntegralPattern` i s odvozením.** Přidat `pattern.formula()` (repr → lidsky čitelný vzorec)
   usnadní debug i budoucí integraci do SymPy; dnes `__str__` tiskne vstup, ne výsledek.
3. **Nemíchat validaci a extrakci.** `validate_power_non_negative` dostává `None` z matcheru a řeší,
   co matcher neumí; čistší je nechat matcher vracet `None` a validovat až výsledný pattern.
4. **`simplify()` v `core.py`** je nejdražší operace (~1.37 ms z README). Zvážit `cancel`/`together`
   nebo `simplify(..., ratio=...)` a měřit; pro porovnání v testech stačí `simplify(result - expected)`.
5. **`as_numer_denom` nestačí** pro všechny zápisy — normalizace (P2-2) patří do samostatné fáze
   před matchem, ne do matcheru samotného.
6. **Konzistence s poznámkami o SymPy.** `notes/sympy/` plánuje Meijer-G patterny a uvádí jen BE
   větev (`polylog(p+1, exp(-b))`); uzavřený tvar pro FD s `b ≠ 0` tam chybí. Po opravě P0-1/P0-2
   doplnit správné formule i do `SYMPY_INTEGRATION_ANALYSIS.md` a `SYMPY_IMPLEMENTATION_PLAN.md`,
   aby se chyba nepřenesla do integrace.
7. **Symetrie BE/FD.** FD `b=0` je `η(s)`, BE `b=0` je `ζ(s)`; pro `p=0` dej oběma konzistentní
   chování (`ln 2` vs. divergence) a otestuj jako první řádek tabulky — je to nejčastější případ.

---

## Poznámka k ověřování

Prostředí spravuje `uv` (Python 3.13.1, SymPy 1.14; `uv sync` instaluje i skupiny `dev` a `examples`):

```bash
uv sync                                  # runtime + dev + examples
uv run pytest -q                         # 100 passed, 8 skipped
uv run python examples/example_usage.py  # π⁴/15
uv build                                 # sdist + wheel (jen bosefermi/)
```

Numerické kros-checky porovnávají symbolické výsledky s `mpmath.quad` (`tests/test_numerics.py`).

---

## Co ještě zbývá

- Rewriting numerátoru s `exp(-x)` před analýzou (P2-2), např. `x**3*exp(-x)/(1−exp(−x))`.
- Zúžit `except Exception` v `core.py` na očekávané výjimky (P1-6); teď je každá vidět jako `RuntimeWarning`.
- Osa SymPy 1.12/1.14 v CI (matice je zatím jen podle Pythonu).
- Sjednocení formulí do jedné `special_function(p, b, sign)` a `IntegralPattern.formula()` (tipy 1–2).
- Commitnout `ROADMAP.md` a `notes/` (nyní untracked).

> Tento soubor zatím **není commitnutý** (záměrně, na pokyn).
