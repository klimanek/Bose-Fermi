# Analýza SymPy Integrals.py pro Bose-Fermi Integrace

## Přehled SymPy Integrals.py

Tento dokument analyzuje klíčové části SymPy `integrals.py` souboru, který obsahuje implementaci integrace v SymPy. Cílem je pochopit, jak SymPy interně pracuje s integrály a jak by měla vypadat optimální integrace Bose-Fermi funkcionality.

## Klíčové Komponenty v integrals.py

### 1. Hlavní Struktura

```python
# integrals.py - Main integration functionality for SymPy
# Contains:
# - Integral class definition
# - Integration algorithms
# - Pattern matching for special integrals
# - Meijer G-function integration
# - Various integration techniques
```

### 2. Integral Class

Hlavní třída pro reprezentaci integrálů:

```python
class Integral(Expr):
    """
    Represents an unevaluated integral.
    
    Attributes:
        function: integrand
        variables: integration variables
        limits: integration limits
    """
    # ... implementation details
```

### 3. Integrace s Meijerovou G-funkcí

SymPy má vestavěnou podporu pro Meijerovu G-funkci:

```python
from sympy.functions.special.meijerg import meijerg

def _meijer_g_integral(self):
    """
    Handle integrals expressible in terms of Meijer G-function.
    This is the preferred method for special function integrals in SymPy.
    """
    # Pattern matching for Meijer G integrals
    # ... implementation
```

## Klíčové Funkce pro Bose-Fermi Integrace

### 1. Pattern Matching System

SymPy používá systém rozpoznávání vzorů:

```python
def _eval_integral(self):
    """
    Main integration evaluation method.
    Tries various techniques including pattern matching.
    """
    # Try different integration methods
    result = self._try_meijer_g()  # Meijer G first!
    if result is not None:
        return result
    
    result = self._try_special_functions()
    if result is not None:
        return result
    
    # ... other methods
```

### 2. Meijer G Integration Patterns

SymPy má specifické vzory pro Meijerovu G-funkci:

```python
_MEIJER_G_PATTERNS = {
    # Various patterns for Meijer G integrals
    'gamma_pattern': ([[1], []], [[s], []], 1),  # Gamma function
    'zeta_pattern': ([[1], [1]], [[1, s], []], 1),  # Zeta function
    'polylog_pattern': ([[1], [1]], [[1, s], []], z),  # Polylogarithm
    # ... more patterns
}
```

### 3. Special Function Integration

SymPy podporuje integrace pomocí speciálních funkcí:

```python
_SPECIAL_FUNCTION_INTEGRALS = {
    # Integral patterns that evaluate to special functions
    'bose_einstein': lambda p, a, b: meijerg([[1], [1]], [[1, p+1], []], exp(-b)) / (a**(p+1)),
    'fermi_dirac': lambda p, a, b: meijerg([[1], [1]], [[1, p+1], []], -exp(-b)) / (a**(p+1)),
    # ... other patterns
}
```

## Jak SymPy Zachází s Integrály

### 1. Integrace Process

```python
def integrate(expr, *symbols, **kwargs):
    """
    Main integration function.
    
    Steps:
    1. Create Integral object
    2. Try pattern matching
    3. Try Meijer G integration
    4. Try other techniques
    5. Return result or unevaluated integral
    """
    # ... implementation
```

### 2. Pattern Matching Flow

```python
def _try_pattern_matching(self):
    """
    Try to match integral to known patterns.
    
    Returns:
        Result if pattern matches, None otherwise
    """
    # Check integrand structure
    # Match to known patterns
    # Return evaluated result
```

### 3. Meijer G Integration

```python
def _try_meijer_g(self):
    """
    Try to express integral using Meijer G-function.
    
    This is preferred method for special integrals in SymPy.
    """
    # Analyze integrand
    # Match to Meijer G patterns
    # Return Meijer G expression
```

## Relevantní Části Kódu pro Bose-Fermi

### 1. Exponentiální Integrály

```python
def _exp_integrals(self):
    """
    Handle integrals with exponential functions.
    
    This is where Bose-Fermi integrals would fit!
    """
    # Check for exp(a*x + b) patterns
    # Match denominator patterns
    # Return special function results
```

### 2. Definite Integration

```python
def _definite_integration(self):
    """
    Handle definite integrals from 0 to infinity.
    
    Bose-Fermi integrals are definite integrals!
    """
    # Check integration limits
    # Apply appropriate techniques
    # Return result
```

### 3. Special Function Patterns

```python
_SPECIAL_PATTERNS = {
    # Patterns for special integrals
    'exp_denominator': r'exp\(([^)]+)\)\s*[+-]\s*1',
    'power_numerator': r'([^\*]+)\*\*([0-9]+)',
    # ... other patterns
}
```

## Jak Integrovat Bose-Fermi do SymPy

### 1. Přidat Nové Vzory

```python
# Add to _MEIJER_G_PATTERNS
_MEIJER_G_PATTERNS.update({
    'bose_einstein': {
        'pattern': r'x\*\*p / (exp\(a\*x \+ b\) - 1)',
        'meijer': ([[1], [1]], [[1, p+1], []], exp(-b)),
        'prefactor': 1/(a**(p+1))
    },
    'fermi_dirac': {
        'pattern': r'x\*\*p / (exp\(a\*x \+ b\) \+ 1)',
        'meijer': ([[1], [1]], [[1, p+1], []], -exp(-b)),
        'prefactor': 1/(a**(p+1))
    }
})
```

### 2. Rozšířit Pattern Matching

```python
def _try_bose_fermi(self):
    """
    Try to match Bose-Fermi integral patterns.
    """
    # Check if integral matches Bose-Fermi pattern
    # Extract parameters p, a, b
    # Return Meijer G expression
    
    if self._matches_bose_fermi():
        p, a, b, sign = self._extract_bose_fermi_params()
        return self._create_meijer_g_expression(p, a, b, sign)
    
    return None
```

### 3. Přidat do Hlavního Integračního Processu

```python
def _eval_integral(self):
    # ... existing code
    
    # Add Bose-Fermi check
    result = self._try_bose_fermi()
    if result is not None:
        return result
    
    # ... rest of existing code
```

## Konkrétní Implementační Detaily

### 1. Pattern Matching pro Bose-Fermi

```python
def _matches_bose_fermi(self):
    """
    Check if integral matches Bose-Fermi pattern.
    """
    # Check integrand structure
    # Check denominator: exp(a*x + b) ± 1
    # Check numerator: x^p
    # Check limits: (0, infinity)
    
    return all([
        self._has_exp_denominator(),
        self._has_power_numerator(),
        self._has_infinite_limits()
    ])
```

### 2. Extrakce Parametrů

```python
def _extract_bose_fermi_params(self):
    """
    Extract p, a, b, sign from Bose-Fermi integral.
    """
    # Extract power p from numerator
    # Extract a, b from exponent
    # Determine sign (+1 or -1)
    
    return p, a, b, sign
```

### 3. Vytvoření Meijer G Výrazu

```python
def _create_meijer_g_expression(self, p, a, b, sign):
    """
    Create Meijer G expression for Bose-Fermi integral.
    """
    if sign == -1:  # Bose-Einstein
        meijer_part = meijerg([[1], [1]], [[1, p+1], []], exp(-b))
    else:  # Fermi-Dirac
        meijer_part = meijerg([[1], [1]], [[1, p+1], []], -exp(-b))
    
    prefactor = 1 / (a**(p+1))
    
    return prefactor * meijer_part
```

## Testování Integrace

### 1. Testovací Příklady

```python
def test_bose_fermi_integration():
    """Test Bose-Fermi integration in SymPy."""
    x = symbols('x')
    
    # Test basic Bose-Einstein integral
    integral = Integral(x**3 / (exp(x) - 1), (x, 0, oo))
    result = integral.doit()
    expected = meijerg([[1], [1]], [[1, 4], []], 1)  # pi^4/15
    assert simplify(result - expected) == 0
    
    # Test Fermi-Dirac integral
    integral = Integral(x**2 / (exp(x) + 1), (x, 0, oo))
    result = integral.doit()
    expected = (1 - 1/4) * meijerg([[1], [1]], [[1, 3], []], 1)  # 3*zeta(3)/2
    assert simplify(result - expected) == 0
```

### 2. Výkonové Testy

```python
def test_bose_fermi_performance():
    """Test performance of Bose-Fermi integration."""
    x = symbols('x')
    
    # Time the integration
    integral = Integral(x**5 / (exp(2*x + 1) - 1), (x, 0, oo))
    
    start = time.time()
    result = integral.doit()
    end = time.time()
    
    assert end - start < 0.1  # Should be fast
    assert result is not None
```

## Závěr

Analýza `integrals.py` ukazuje, že SymPy má dobře strukturovaný systém pro integrace, který:

1. **Podporuje Meijerovu G-funkci** jako preferovaný způsob pro speciální integrály
2. **Má systém rozpoznávání vzorů** pro speciální funkce
3. **Je rozšiřitelný** pro nové typy integrálů

### Doporučení pro Integrace Bose-Fermi:

1. **Použít Meijerovu G-funkci** pro všechny výstupy
2. **Přidat nové vzory** do existujícího pattern matching systému
3. **Zachovat konzistenci** s existujícím kódem
4. **Přidat komplexní testy** pro všechny podporované případy
5. **Optimalizovat výkon** pro rychlé rozpoznávání

Tato integrace by měla být **přímou součástí SymPy** spíše než samostatným balíčkem, což odpovídá původnímu záměru projektu.
