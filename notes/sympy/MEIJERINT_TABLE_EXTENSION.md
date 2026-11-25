# Rozšíření tabulek v `meijerint.py` o Bose–Fermiho integrál

> **TL;DR** — Do `_create_lookup_table()` (tabulka „funkce → Meijerova G") tenhle integrál
> obecně přidat **nelze**: tabulka umí jen **konečné** kombinace G funkcí, ale Mellinův obraz
> $1/(e^{ax+b}\pm 1)$ je $a^{-s}\Gamma(s)\mathrm{Li}_s(\pm e^{-b})$. Pro $b=0$ jde o
> $\Gamma(s)\zeta(s)$, resp. $\Gamma(s)\eta(s)$; to **není** Meijerova G funkce. Reálné místo
> je háček v `_meijerint_definite_2()` (nebo v `meijerint_definite()`), který vrátí hotový
> uzavřený tvar.
> `meijerint_definite()` už dnes výsledek bez `meijerg` **přijímá a vrací** (větev
> `if _has(res[0], meijerg) … else: return res`), takže se to dá udělat čistě.

Dokument vychází ze SymPy **1.14** (nainstalovaná verze, odkud jsou čísla řádků) a ověřeně stejné
struktury na `master` ([meijerint.py](https://github.com/sympy/sympy/blob/master/sympy/integrals/meijerint.py),
listopad 2026; čísla řádků ±2).

---

## 1. Co je v `meijerint.py` „tabulka"

Soubor nemá datový soubor s tabulkami — tabulka se staví v kódu:

| funkce | řádek | role |
|---|---|---|
| `_create_lookup_table(table)` | ~88 | plní `_lookup_table` dvojicemi `f(x) → G` |
| `add(formula, an, ap, bm, bq, arg, fac, cond, hint)` | ~97 | přidá záznam `f(x) = fac · G(...)` |
| `addi(formula, terms, cond, hint)` | ~102 | totéž pro obecný seznam členů |
| `_mytype(f, x)` | ~301 | klíč tabulky: množina typů funkcí v `f` (např. `(exp,)`) |
| `_rewrite_single(f, x)` | ~1470 | najde v tabulce, jinak zkusí rekurzivní Mellin |
| `_rewrite1` / `_rewrite2` | ~1610 / ~1624 | rozklad na jednu / dvě G funkce |
| `_meijerint_definite_4(f, x)` | ~2010 | Saxenova věta: `_int0oo_1` / `_int0oo` + `hyperexpand` |
| `meijerint_definite(f, x, a, b)` | ~1781 | vstupní bod z `Integral.doit()` |

Záznam v tabulce má tvar:

```python
# _create_lookup_table(): každý záznam = (formula, [(fac, meijerg(...))], cond, hint)
# skutečný záznam z meijerint.py (~ř. 162): e^{-t} = G^{1,0}_{0,1}(t | 0)
add(exp(polar_lift(-1)*t), [], [], [0], [])
```

Tj. **levá strana je funkce `proměnné z` a pravá strana je výhradně `meijerg(...)`** (případně
seznam několika G funkcí). `_rewrite_single` se na tabulku ptá přes `_mytype(f, z)` a zkouší
`f.match(formula, old=True)`; jakmile najde shodu, vrátí `[(C, s, G), …]`.

Pro cestu `∫₀^∞` pak `_meijerint_definite_4` aplikuje Saxenovu větu (`_int0oo_1` pro jednu G,
`_int0oo` pro dvě) a výsledek `hyperexpand`uje. **V samotné G pipeline neexistuje jiná cesta
k výsledku než přes G funkce** — jediná výjimka je háček z §3.1, který hotový výsledek vrátí
ještě před ní.

### Kdy se to vůbec zavolá

`Integral.doit()` v `sympy/integrals/integrals.py` volá `try_meijerg()` (~ř. 568) jen pro
definitivní integrál s alespoň jednou nekonečnou mezí a jen pokud `meijerg is not False`.
Náš integrál se tam dostane, ale `meijerint_definite` dnes vrátí `None`:

```python
>>> from sympy.integrals.meijerint import meijerint_definite
>>> meijerint_definite(x**3/(exp(x) - 1), x, 0, oo) is None
True
>>> integrate(x**3/(exp(x) - 1), (x, 0, oo))
Integral(x**3/(exp(x) - 1), (x, 0, oo))
```

---

## 2. Proč to do tabulky nejde (obecně)

Tabulka vyžaduje, aby bylo možné integrand napsat jako **konečný** součin/součet G funkcí.
Pro náš integrál to neplatí:

1. **Mellinův obraz není G funkce.** Pro $b=0$:
   $\mathcal{M}\{1/(e^x-1)\}(s)=\Gamma(s)\zeta(s)$ a
   $\mathcal{M}\{1/(e^x+1)\}(s)=\Gamma(s)\eta(s)$.
   SymPy to ani nespočítá:

   ```python
   >>> mellin_transform(1/(exp(x) - 1), x, s)
   MellinTransform(1/(exp(x) - 1), x, s)
   >>> inverse_mellin_transform(gamma(s)*zeta(s), s, x, (1, oo),
   ...                          as_meijerg=True, needeval=True)
   IntegralTransformError: Inverse Mellin Transform could not be computed: needeval.
   ```

   `_rewrite_single` končí na `mellin_transform(...)` + `_inverse_mellin_transform(...)`; když
   inverze neumí vrátit G funkci, cesta končí.

2. **Rozvoj je nekonečná řada.** $1/(e^{ax+b}-1)=\sum_{n\ge1} e^{-n(ax+b)}$; `_rewrite_single`
   vrací **konečný** seznam `(C, s, G)`, ne nekonečnou sumu.

3. **Výsledek je mimo třídu G funkcí.** $\Gamma(p+1)\zeta(p+1)$ / $\Gamma(p+1)\mathrm{Li}_{p+1}$
   nejde `hyperexpand`ovat z G funkce; i kdyby se integrand na G přepsal, výsledek by musel projít
   post-processingem (a pro obecné $p$ to možné není).

**Důsledek:** přidat `1/(exp(z)-1)` do `_create_lookup_table` jako další záznam *nejde* — `add()`
očekává `meijerg(...)` a žádné takové parametry neexistují. (Pro speciální $p$, kde $\zeta(2n)$
vyjde jako $\pi^{2n}$, to také neprojde: integrand `x**3/(e^x-1)` sám G funkce není.)

---

## 3. Kde to jde rozšířit

### 3.1 Doporučeno: háček v `_meijerint_definite_2`

`meijerint_definite` pro mez `(0, ∞)` dělá:

```python
# ~ř. 1865
elif (a, b) == (S.Zero, S.Infinity):
    res = _meijerint_definite_2(f, x)
    if res:
        if _has(res[0], meijerg):
            results.append(res)
        else:
            return res        # ← výsledek bez G funkce se vrací rovnou
```

Takže stačí, aby `_meijerint_definite_2` umělo náš tvar vyřešit **dřív**, než začne zkoušet
`_guess_expansion` a G pipeline:

```python
# sympy/integrals/meijerint.py, _meijerint_definite_2 (za x = dummy)
from .bosefermi import bose_fermi_definite

special = bose_fermi_definite(f, x)
if special is not None:
    return special          # (vysledek, podminka)
```

Výhody: minimální zásah, respektuje `meijerg=True/False` i `conds=` z `Integral.doit()`,
nechá `Piecewise` logiku na stávajícím kódu. Nevýhoda: do modulu, který je „o G funkcích",
přidává ne-G větev — upstream to může chtít vidět spíš jako samostatný modul s explicitním
komentářem (viz 3.2).

### 3.2 Háček v `integrals.py` (`try_meijerg`)

Alternativa je zavolat vlastní pravidlo v `Integral.doit()` před `try_meijerg` (~ř. 568–580):

```python
if len(xab) == 3 and (xab[1], xab[2]) == (S.Zero, S.Infinity):
    ret = bose_fermi_definite(function, xab[0])
    if ret is not None:
        function, cond = ret
        ...
```

Je to čistší oddělení, ale musí se ručně ošetřit `conds='piecewise'` a `meijerg=False`,
které jinak řeší `try_meijerg`.

### 3.3 Mellinova transformace v `transforms.py` (samostatně užitečné)

Do tabulky v `sympy/integrals/transforms.py` (`mellin_transform`) lze doplnit

$$
\mathcal{M}\left\{\frac{1}{e^x-1}\right\}(s)=\Gamma(s)\zeta(s),\qquad
\mathcal{M}\left\{\frac{1}{e^x+1}\right\}(s)=\Gamma(s)\eta(s),
$$

s pásmem konvergence $1<\Re s<\infty$ (bosony) resp. $0<\Re s<\infty$ (fermiony). Tím se opraví
`mellin_transform`, ale **meijerint to nepomůže** — inverze $\Gamma(s)\zeta(s)$ na G funkci
nevede (viz §2). Ber to jako doplněk, ne jako řešení.

### 3.4 Co nedělat

- Přidávat do `_create_lookup_table` „falešný" G záznam pro `1/(e^z-1)` — `_int0oo_1`
  a `_check_antecedents_1` by dostaly něco, co není G funkce, a podmínky konvergence by nedávaly smysl.
- Přidávat záznamy pro konkrétní `p` (např. `x**3/(e^x-1)`) — tabulka je o funkcích, ne o
  určitých integrálech, a výsledek $\pi^4/15$ z G funkce nevyjde.

---

## 4. Matematika, kterou má háček vracet

Pro $a>0$; u bosonů navíc $p>0$, u fermionů $p>-1$ (konvergence u $x=0$):

| integrál | výsledek |
|---|---|
| $\int_0^\infty \frac{x^p}{e^{ax+b}-1}dx$ | $\dfrac{\Gamma(p+1)}{a^{p+1}}\,\mathrm{Li}_{p+1}(e^{-b})$, $b\ge 0$ |
| $\int_0^\infty \frac{x^p}{e^{ax+b}+1}dx$ | $-\dfrac{\Gamma(p+1)}{a^{p+1}}\,\mathrm{Li}_{p+1}(-e^{-b})$, libovolné reálné $b$ |
| $b=0$, bosony | $\Gamma(p+1)\zeta(p+1)$ |
| $b=0$, fermiony | $\Gamma(p+1)\eta(p+1)$ (v SymPy `dirichlet_eta(p+1)`) |
| $p=0$, fermiony | $\ln 2$ — `dirichlet_eta(1)` to vrátí sám, `(1-2^{-p})\zeta(p+1)` nikoli |

Poznámky:

- **Škálovaný jmenovatel.** SymPy běžně rozvine `exp(2*x + 0.3)` na `1.3498…*exp(2*x)`, takže
  se hodí matchnout obecné $A e^{ax+b}+C$ s číselnými $A,C$ a převést
  $A e^{ax+b}+C=C\left(e^{ax+b+\log|A/C|}\pm1\right)$; faktor $1/C$ (nebo $-1/C$) jde do předfaktoru.
  Přesně to dělá `PatternMatcher.match_exponential_denominator()` v tomto repu.
- **$b<0$ u fermionů** znamená $\mathrm{Li}_{p+1}$ mimo jednotkový kruh → analytické pokračování.
  SymPy `polylog` to zvládne, ale u numerického `evalf` může zůstat šumová imaginární složka
  (~1e-28); je to jen numerika, ne chyba.
- **$p=0$, $b\ne0$**: `polylog(1, -exp(-b))` se symbolickým `b` nezjednoduší; vyplatí se vrátit
  rovnou `log(1 + exp(-b))`.
- **Symbolické parametry bez předpokladů**: neodmítat, ale vrátit `cond`
  (`And(a > 0, p > -1, …)`), aby si `Integral.doit()` mohl udělat `Piecewise`.

---

## 5. Návrh kódu

Nový modul `sympy/integrals/bosefermi.py` (≈70 řádků), bez závislosti na tomto repu:

```python
from sympy import And
from sympy.core import Add, S, Wild
from sympy.functions import exp, log
from sympy.functions.special.gamma_functions import gamma
from sympy.functions.special.zeta_functions import dirichlet_eta, polylog, zeta


def _match(f, x):
    """f = c*x**p / (A*exp(alpha*x + mu) + C) -> (c, p, alpha, mu, sign, scale) | None."""
    num, den = f.as_numer_denom()
    coeff, num = num.as_coeff_Mul()

    if num == 1:
        p = S.Zero
    elif num == x:
        p = S.One
    elif num.is_Pow and num.base == x:
        p = num.exp
    else:
        return None

    # jmenovatel: právě dva členy, jeden A*exp(...), druhý číselná konstanta C
    terms = Add.make_args(den)
    if len(terms) != 2:
        return None
    exp_term = const = None
    for term in terms:
        if term.as_coeff_Mul()[1].func == exp:
            if exp_term is not None:
                return None
            exp_term = term
        elif term.is_number:
            if const is not None:
                return None
            const = term
        else:
            return None
    if exp_term is None or const is None:
        return None

    A, exp_part = exp_term.as_coeff_Mul()
    if not A.is_number or A.is_zero:
        return None

    # exponent alpha*x + mu
    alpha_, mu_ = Wild('alpha', exclude=[x]), Wild('mu', exclude=[x])
    m = exp_part.args[0].match(alpha_*x + mu_)
    if m:
        alpha, mu = m[alpha_], m.get(mu_, S.Zero)
    elif exp_part.args[0] == x:
        alpha, mu = S.One, S.Zero
    else:
        return None

    # A*exp(...) + C = C*(exp(... + log|A/C|) + sign)
    ratio = A / const
    if ratio.is_positive:
        sign, shift, scale = S.One, log(ratio), S.One / const
    elif ratio.is_negative:
        sign, shift, scale = S.NegativeOne, log(-ratio), -S.One / const
    else:
        return None

    return coeff, p, alpha, mu + shift, sign, scale


def bose_fermi_definite(f, x):
    """Vrátí (vysledek, cond) pro ∫₀^∞ f dx, nebo None."""
    m = _match(f, x)
    if m is None:
        return None
    coeff, p, alpha, mu, sign, scale = m

    cond = And(alpha > 0, p + 1 > 0)
    if sign == S.NegativeOne:                  # bosony: konvergence i u nuly
        cond = And(cond, p > 0, mu >= 0)

    pref = coeff * scale * gamma(p + 1) / alpha**(p + 1)
    if sign == S.One:                          # fermiony
        if mu == 0:
            special = dirichlet_eta(p + 1)     # p = 0 -> log(2) samo
        elif p == 0:
            special = log(1 + exp(-mu))
        else:
            special = -polylog(p + 1, -exp(-mu))
    else:                                      # bosony
        special = zeta(p + 1) if mu == 0 else polylog(p + 1, exp(-mu))

    return pref * special, cond
```

Háček (varianta 3.1):

```python
# sympy/integrals/meijerint.py
def _meijerint_definite_2(f, x):
    ...
    x = dummy
    if f == 0:
        return S.Zero, True

    from .bosefermi import bose_fermi_definite
    special = bose_fermi_definite(f, x)
    if special is not None:
        _debug('Bose-Fermi integral recognised')
        return special

    for g, explanation in _guess_expansion(f, x):
        ...
```

Podmínky se vracejí jako `cond` (ne jako `bool`); `Integral.doit()` je při `conds='piecewise'`
zabalí do `Piecewise`, při `conds='separate'` vrátí jako druhou položku. Pro číselné parametry
se `And(3 > 0, 4 > 0)` samo vyhodnotí na `True`.

---

## 6. Testy a ověření

Do `sympy/integrals/tests/test_meijerint.py` (nebo nového `test_bosefermi.py`):

```python
from sympy import exp, gamma, integrate, log, pi, symbols, zeta
from sympy.abc import x
from sympy.integrals.meijerint import meijerint_definite

a, p = symbols('a p', positive=True)

res, cond = meijerint_definite(x**3/(exp(x) - 1), x, 0, oo)
assert res == pi**4/15
assert cond == True        # číselné parametry -> podmínka se vyhodnotí

assert integrate(x**3/(exp(x) - 1), (x, 0, oo)) == pi**4/15
assert integrate(1/(exp(x) + 1), (x, 0, oo)) == log(2)
assert integrate(x**p/(exp(a*x) - 1), (x, 0, oo)) == gamma(p+1)*zeta(p+1)/a**(p+1)
```

Numericky se to dá ověřit stejně jako v tomto repu (`tests/test_numerics.py`):
symbolický výsledek vs. `mpmath.quad` s tolerancí ~1e-10, pro
`p ∈ {0, 1/2, 1, 3/2, 2, 5/2}`, `a ∈ {1, 2}`, `b ∈ {-1, 0, 1}`, oba znaky.

Spuštění v SymPy: `python -m pytest sympy/integrals/tests/test_meijerint.py -q`
(případně `bin/test integrals`).

---

## 7. Úskalí a doporučení pro upstream

- **Výkon.** `meijerint` je drahý (`_meijerint_definite_2` zkouší několik expanzí). Náš test je
  O(1) a patří **před** `_guess_expansion`, ne za něj.
- **`meijerg=False`.** Pokud se háček dá do `integrals.py`, musí respektovat `meijerg=False`;
  v `_meijerint_definite_2` je to automatické.
- **Podmínky.** Nevracet výsledek pro neplatné parametry (`a ≤ 0`, bosony s `b < 0`, `p ≤ -1`);
  vrátit `None`, ať se pokračuje ostatními metodami. Pro symbolické parametry radši `cond`.
- **Branch cuts.** U fermionů s $b<0$ je potřeba analytické pokračování `polylog`; v testech
  tolerovat malou imaginární složku.
- **Žádná závislost na `bosefermi`.** Do SymPy se portuje ~60–80 řádků; logika matcheru
  (`A*exp(a*x+b)+C`, `a*x+b` přes `Wild`) je v tomto repu odladěná a přenositelná.
- **Kam to patří.** `meijerint.py` je o G funkcích; upstream může preferovat nový modul
  (`sympy/integrals/bosefermi.py`) a explicitní volání z `Integral.doit()` nebo z
  `meijerint_definite` s komentářem „known non-G definite integral". Pro `manualintegrate`
  to není — ten staví primitivní funkce, kdežto tady jde o určitý integrál bez elementární primitivy.
- **Dokumentace.** Doplnit řádek do tabulky v `doc/src/modules/integrals/integrals.rst`
  a do docstringu `integrate()`.

---

## 8. Odkazy

- [sympy/sympy — `sympy/integrals/meijerint.py`](https://github.com/sympy/sympy/blob/master/sympy/integrals/meijerint.py)
- Tabulky v `_create_lookup_table` odkazují na **[P]** = Prudnikov, Brychkov, Marichev,
  *Integrals and Series*, Vol. 3 (sekce 8.4) a **[L]** = Luke, *The Special Functions and
  Their Approximations*.
- Tento repozitář: `bosefermi/analyzers.py` (matcher), `bosefermi/special_functions.py`
  (uzavřené tvary), `tests/test_numerics.py` (kros-checky), `ROADMAP.md` (P0-1/P0-2).
