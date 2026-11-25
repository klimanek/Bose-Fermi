# Analýza Integrace Bose-Fermi do SymPy

## Klíčový Kontext: Vztah k SymPy

### Proč Projekt Vznikl

Bose-Fermi projekt vznikl jako **dočasné řešení** pro specifickou třídu integrálů, které **SymPy aktuálně neumí evaluovat**:

```python
# Příklad integrálu, který SymPy neumí vyřešit:
from sympy import symbols, Integral, exp, oo
x = symbols('x')
integral = Integral(x**3 / (exp(x) - 1), (x, 0, oo))
result = integral.doit()  # → Vrací nevyhodnocený integrál!
```

### Cíl: Integrace do SymPy

**Hlavní cíl projektu je integrace do SymPy**, nikoliv být samostatným balíčkem. Projekt slouží jako:
- **Důkaz konceptu** (proof of concept)
- **Testovací platforma** pro rozpoznávání vzorů
- **Základ pro budoucí implementaci** v SymPy

## Problém s Meijerovou G-funkcí

### Co je Meijerova G-funkce

Meijerova G-funkce je **zobecněná hypergeometrická funkce**, která může reprezentovat širokou třídu speciálních funkcí:


$$G_{p,q}^{m,n}
\left(
  z \, \bigg| \, \begin{matrix} a_1, \dots, a_p \\ b_1, \dots, b_q \end{matrix}
\right)$$


**Důležité**: SymPy používá **stejnou notaci jako Mathematica**, která se liší od některých literárních konvencí. V SymPy je funkce implementována jako `meijerg([[a1, ..., ap], [b1, ..., bq]], [[c1, ..., cm], [d1, ..., dn]], z)` kde parametry odpovídají specifické SymPy/Mathematica konvenci.

### Proč SymPy Vývojáři Preferují Meijerovu G-funkci

1. **Unifikovaný přístup**: Meijerova G-funkce může reprezentovat mnoho speciálních funkcí (zeta, polylog, hypergeometrické funkce) v jedné formě
2. **Konzistence**: SymPy již používá Meijerovu G-funkci pro jiné speciální funkce
3. **Rozšiřitelnost**: Snadnější přidávání nových vzorů a funkcí
4. **Matematická elegance**: Jednotný framework pro speciální funkce

### Aktuální Implementace vs. Meijerova G-funkce

**Aktuální implementace** používá přímé funkce:
- `gamma(p+1)` - Gamma funkce
- `zeta(p+1)` - Riemannova zeta funkce
- `polylog(p+1, z)` - Polylogaritmus
- `eta(p+1)` - Dirichletova eta funkce

**Požadovaná implementace** by měla používat:
- `meijerg([[...], [...]], [[...], [...]], z)` - Meijerova G-funkce

## Technické Detaily

### Jak by Měla Vypadat Integrace

#### 1. Překlad aktuálních funkcí na Meijerovu G-funkci

```python
# Aktuální implementace:
def bose_einstein_special_function(p, b):
    if b == 0:
        return zeta(p + 1)
    else:
        return polylog(p + 1, exp(-b))

# Požadovaná implementace s Meijerovou G-funkcí:
def bose_einstein_special_function(p, b):
    if b == 0:
        # zeta(p+1) = G_{1,2}^{2,1}
        return meijerg([[1], [1]], [[1, p+1], []], 1)
    else:
        # polylog(p+1, z) = G_{1,2}^{2,1}
        return meijerg([[1], [1]], [[1, p+1], []], exp(-b))
```

#### 2. Překlad vzorců

**Aktuální vzorec:**
```
∫₀^∞ x^p / (e^(αx) - 1) dx = Γ(p+1) · ζ(p+1) / α^(p+1)
```

**S Meijerovou G-funkcí:**
```
∫₀^∞ x^p / (e^(αx) - 1) dx = 
  G_{1,2}^{2,1}
  (
    1 \, \bigg| \, \begin{matrix} 1 \\ 1, p+1 \end{matrix}
  ) 
  · 
  G_{0,1}^{1,0}
  (
    α^{p+1} \, \bigg| \, \begin{matrix} - \\ p+1 \end{matrix}
  )
```

### Výhody Použití Meijerovy G-funkce

1. **Unifikace**: Jedna funkce pro všechny speciální případy
2. **SymPy konzistence**: Lepší integrace s existujícími SymPy funkcem
3. **Rozšiřitelnost**: Snadnější přidávání nových vzorů
4. **Matematická přesnost**: Přesná reprezentace všech speciálních funkcí

## Implementační Strategie

### 1. Překlad Speciálních Funkcí

Vytvořit překladovou tabulku:

| Funkce | Meijerova G-funkce |
|--------|-------------------|
| ζ(s) | G_{1,2}^{2,1}(1 | 1; 1, s) |
| Li_s(z) | G_{1,2}^{2,1}(z | 1; 1, s) |
| η(s) | G_{1,2}^{2,1}(1 | 1; 1, s) - 2^{-s} G_{1,2}^{2,1}(1 | 1; 1, s) |
| Γ(s) | G_{0,1}^{1,0}(1 | -; s) |

### 2. Úprava Special_functions.py

```python
from sympy import meijerg

def gamma_prefactor(p, a):
    """Compute gamma function prefactor using Meijer G-function
    
    Note: SymPy uses Mathematica notation for Meijer G-function:
    meijerg([[a1, ..., ap], [b1, ..., bq]], [[c1, ..., cm], [d1, ..., dn]], z)
    """
    return meijerg([[1], []], [[p+1], []], 1) / meijerg([[1], []], [[p+1], []], a**(p+1))

def bose_einstein_special_function(p, b):
    """Compute using Meijer G-function (SymPy/Mathematica notation)"""
    if b == 0:
        # zeta(p+1) = G_{1,2}^{2,1}(1 | 1; 1, p+1) in SymPy notation
        return meijerg([[1], [1]], [[1, p+1], []], 1)
    else:
        # polylog(p+1, z) = G_{1,2}^{2,1}(z | 1; 1, p+1) in SymPy notation
        return meijerg([[1], [1]], [[1, p+1], []], exp(-b))

def fermi_dirac_special_function(p, b):
    """Compute using Meijer G-function (SymPy/Mathematica notation)"""
    if b == 0:
        # eta(p+1) = (1 - 1/2^p) * zeta(p+1)
        return (meijerg([[1], [1]], [[1, p+1], []], 1) - 
                meijerg([[1], [1]], [[1, p+1], []], 1) / (2**p))
    else:
        # polylog(p+1, -z) = G_{1,2}^{2,1}(-z | 1; 1, p+1) in SymPy notation
        return meijerg([[1], [1]], [[1, p+1], []], -exp(-b))
```

### 3. Testování a Validace

Vytvořit testy, které ověřují ekvivalenci s důrazem na SymPy/Mathematica notaci:

```python
def test_meijer_equivalence():
    """Test that Meijer G-function (SymPy notation) gives same results as direct functions"""
    p = S(3)
    
    # Test gamma function using SymPy/Mathematica notation
    gamma_direct = gamma(p+1)
    gamma_meijer = meijerg([[1], []], [[p+1], []], 1)  # G_{0,1}^{1,0}(1 | -; p+1)
    assert simplify(gamma_direct - gamma_meijer) == 0
    
    # Test zeta function using SymPy/Mathematica notation
    zeta_direct = zeta(p+1)
    zeta_meijer = meijerg([[1], [1]], [[1, p+1], []], 1)  # G_{1,2}^{2,1}(1 | 1; 1, p+1)
    assert simplify(zeta_direct - zeta_meijer) == 0
    
    # Test polylog function
    z = S.Half
    polylog_direct = polylog(p+1, z)
    polylog_meijer = meijerg([[1], [1]], [[1, p+1], []], z)  # G_{1,2}^{2,1}(z | 1; 1, p+1)
    assert simplify(polylog_direct - polylog_meijer) == 0
```

## Výzvy a Řešení

### 1. Výkon

**Problém**: Meijerova G-funkce může být pomalejší než přímé funkce

**Řešení**:
- Použít optimalizované implementace
- Cacheování výsledků
- Použít přímé funkce pro jednoduché případy

### 2. Čitelnost

**Problém**: Meijerova G-funkce je méně čitelná než přímé funkce

**Řešení**:
- Vytvořit wrapper funkce pro běžné případy
- Dokumentace s příklady
- Komentáře v kódu

### 3. Kompatibilita

**Problém**: Starší verze SymPy nemusí mít plnou podporu Meijerovy G-funkce

**Řešení**:
- Zajistit kompatibilitu s novějšími verzemi SymPy
- Podmíněné importy
- Fallback na přímé funkce

## Doporučení pro Implementaci

### 1. Postupná Migrace

1. **Fáze 1**: Přidat podporu pro Meijerovu G-funkci vedle stávajících funkcí
2. **Fáze 2**: Vytvořit překladovou vrstvu
3. **Fáze 3**: Postupně nahradit přímé funkce Meijerovou G-funkcí
4. **Fáze 4**: Optimalizace a testování

### 2. Dokumentace

- Vytvořit dokumentaci s příklady použití Meijerovy G-funkce
- Přidat matematické odvození
- Vytvořit tabulky ekvivalence

### 3. Testování

- Vytvořit komplexní testy pro všechny podporované integrály
- Otestovat hraniční případy
- Ověřit numerickou stabilitu

## Závěr

Přechod na Meijerovu G-funkci je **klíčový pro integraci do SymPy**. Tento přechod přinese:

1. **Lepší integrace**: Konzistence s existujícími SymPy funkcem
2. **Rozšiřitelnost**: Snadnější přidávání nových vzorů
3. **Unifikaci**: Jednotný přístup k speciálním funkcím
4. **Budoucnost**: Základ pro další rozvoj

**Doporučení**: Zahájit migraci na Meijerovu G-funkci s důrazem na:
- Postupnou implementaci
- Komplexní testování
- Dokumentaci
- Kompatibilitu s existujícím kódem

Tato migrace je **nezbytná pro úspěšnou integraci do SymPy** a měla by být prioritou pro další vývoj projektu.
