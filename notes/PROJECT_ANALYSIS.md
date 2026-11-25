# Projekt Bose-Fermi: Komplexní Analýza

## Přehled Projektu

Bose-Fermi je specializovaný Python balíček pro symbolickou evaluaci integrálů, které se vyskytují v statistické fyzice, konkrétně v Bose-Einsteinově a Fermi-Diracově statistice. Projekt je navržen pro integraci do SymPy a poskytuje přesné výsledky pro třídu integrálů, které jsou klíčové v kvantové statistické mechanice.

## Základní Informace

- **Název**: Bose-Fermi Integral Evaluator for SymPy
- **Verze**: 0.1.0
- **Autor**: Oldřich Klimánek
- **Licence**: MIT
- **Závislosti**: SymPy, NumPy, Matplotlib, SciPy
- **Python**: 3.10+

## Architektura Projektu

```
bosefermi/
├── __init__.py          # Hlavní rozhraní
├── core.py             # Jádro - hlavní funkce bose_fermi_integral()
├── analyzers.py        # Analýza a rozpoznávání vzorů
├── patterns.py         # Reprezentace integrálních vzorů
├── validators.py       # Validace vstupů
└── special_functions.py # Matematické funkce (gamma, zeta, polylog)

tests/
├── test_bosefermi.py   # Komplexní testy komponent
└── test_sanity.py      # Sanitní testy

examples/
├── example_usage.py    # Základní použití
├── blackbody.py        # Planckův zákon
├── gas.py              # Fermiho plyn
└── white_dwarf.py      # Astrofyzikální aplikace
```

## Klíčové Komponenty

### 1. Jádro (core.py)

Hlavní funkce `bose_fermi_integral(expr, debug=False)`:
- Přijímá SymPy Integral objekt
- Analyzuje a rozpoznává vzory
- Vrací přesný výsledek nebo původní integrál
- Podpora pro debugování

### 2. Analyzátory (analyzers.py)

Třídy pro rozbor integrálů:
- `IntegralAnalyzer`: Hlavní analyzátor
- `PatternMatcher`: Rozpoznávání vzorů
- Podpora pro extrakci mocnin a koeficientů
- Rozpoznávání exponentiálních jmenovatelů

### 3. Vzory (patterns.py)

Reprezentace rozpoznaných integrálů:
- `IntegralPattern`: Třída pro uložení parametrů
- Podpora pro Bose-Einstein (sign = -1)
- Podpora pro Fermi-Dirac (sign = +1)
- Metody pro evaluaci a reprezentaci

### 4. Validátory (validators.py)

Kontrola vstupních podmínek:
- Validace typu integrálu
- Kontrola proměnných a limit
- Ověření mocnin a koeficientů
- Komplexní validace `validate_all()`

### 5. Speciální Funkce (special_functions.py)

Matematické funkce:
- `gamma_prefactor(p, a)`: Gamma funkce
- `bose_einstein_special_function(p, b)`: Zeta/polylog
- `fermi_dirac_special_function(p, b)`: Eta/polylog

## Podporované Integrály

Projekt rozpoznává a evaluuje integrály tvaru:

```
∫₀^∞ x^p / (e^(αx + μ) ± 1) dx
```

Kde:
- `p`: Mocnina integrační proměnné (p ≥ 0)
- `α`: Kladný škálovací faktor (inverzní teplota)
- `μ`: Chemický potenciál
- `-`: Bose-Einstein statistika
- `+`: Fermi-Dirac statistika

### Podporované varianty:

| Integral | Výsledek |
|---------|--------|
| ∫₀^∞ x^p / (e^x - 1) dx | Γ(p+1)·ζ(p+1) |
| ∫₀^∞ x^p / (e^(αx) - 1) dx | Γ(p+1)·ζ(p+1)/α^(p+1) |
| ∫₀^∞ x^p / (e^(x+μ) - 1) dx | Γ(p+1)·Li_{p+1}(e^{-μ}) |
| ∫₀^∞ x^p / (e^(αx+μ) - 1) dx | Γ(p+1)·Li_{p+1}(e^{-μ})/α^(p+1) |
| ∫₀^∞ x^p / (e^x + 1) dx | Γ(p+1)·η(p+1) |
| ∫₀^∞ x^p / (e^(αx) + 1) dx | Γ(p+1)·η(p+1)/α^(p+1) |

## Fyzikální Aplikace

### 1. Planckův Zákon (blackbody.py)
- Výpočet celkové energie vyzařované absolutně černým tělesem
- Stefan-Boltzmannův zákon
- Aplikace v astrofyzice a termodynamice

### 2. Fermiho Plyn (gas.py)
- Energie volných elektronů v kovech
- Teplotní závislost elektronového plynu
- Aplikace v fyzice pevných látek

### 3. Bílé Trpaslíky (white_dwarf.py)
- Struktura a stabilita bílých trpaslíků
- Chandrasekharova mez (~1.4 M☉)
- Degenerovaný elektronový plyn

## Příklady Použití

### Základní příklad:

```python
from sympy import symbols, oo, Integral, exp
from bosefermi import bose_fermi_integral

x = symbols('x', real=True, positive=True)
integral = Integral(x**3 / (exp(x) - 1), (x, 0, oo))
result = bose_fermi_integral(integral)
# Výsledek: π⁴/15
```

### Symbolické parametry:

```python
from sympy import symbols, oo, Integral, exp, gamma, zeta
from bosefermi import bose_fermi_integral

x, a, p = symbols('x a p', positive=True)
integral = Integral(x**p / (exp(a*x) - 1), (x, 0, oo))
result = bose_fermi_integral(integral)
# Výsledek: gamma(p+1)*zeta(p+1)/a^(p+1)
```

## Testování

### Komplexní testy (test_bosefermi.py)
- Testy validátorů
- Testy rozpoznávání vzorů
- Testy analyzátorů
- Testy integrálních vzorů
- Integrační testy

### Sanitní testy (test_sanity.py)
- Základní integrály
- Symbolické parametry
- Hraniční podmínky
- Validace vstupů
- Výkonové testy
- Známé speciální hodnoty
- Konzistence

## Výkon

- Typická doba evaluace: ~1.37 ms
- Rychlé rozpoznávání vzorů
- Efektivní validace
- Optimalizováno pro SymPy integraci

## Budoucí Rozvoj

### Plánované funkce:
- Rozšíření o další vzory integrálů
- Vylepšená podpora pro symbolické parametry
- Integrace do SymPy
- Rozšířená dokumentace
- Vylepšené příklady

### Možnosti rozšíření:
- Podpora pro komplexní integrály
- Rozšíření o další speciální funkce
- Vylepšené rozpoznávání vzorů
- Optimalizace výkonu

## Integrace do SymPy

Projekt je navržen pro integraci do SymPy:
- Kompatibilní rozhraní
- Podpora SymPy objektů
- Symbolické výpočty
- Přesné výsledky

## Závěr

Bose-Fermi je dobře strukturovaný, modulární projekt s jasným rozhraním a podporou pro debugování. Poskytuje přesné výsledky pro důležité třídy integrálů v statistické fyzice a má potenciál pro rozšíření o další vzory a funkce. Projekt je připraven pro integraci do SymPy a může být užitečný pro výzkumníky a studenty v oblasti fyziky a matematiky.
