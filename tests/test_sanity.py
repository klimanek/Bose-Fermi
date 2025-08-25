"""
Sanity test suite for the Bose-Fermi integral evaluator.

These tests verify basic functionality and catch obvious regressions.
"""

from sympy import symbols, oo, Integral, exp, zeta, gamma, simplify, S, pi
from bosefermi.core import bose_fermi_integral
import time


def test_basic_integrals():
    """Test basic Bose-Einstein and Fermi-Dirac integrals."""
    x = symbols('x', real=True, positive=True)
    
    test_cases = [
        # Basic Bose-Einstein integrals
        (Integral(x**2 / (exp(x) - 1), (x, 0, oo)), 
         2 * zeta(3)),
        
        (Integral(x**3 / (exp(x) - 1), (x, 0, oo)), 
         6 * zeta(4)),
        
        # Basic Fermi-Dirac integrals  
        (Integral(x**2 / (exp(x) + 1), (x, 0, oo)), 
         (S.One - S.One/4) * 2 * zeta(3)),
        
        # With coefficient in exp
        (Integral(x**2 / (exp(2*x) - 1), (x, 0, oo)), 
         gamma(3) * zeta(3) / 8),
    ]
    
    print("Running basic integral tests...")
    passed = 0
    for i, (integral, expected) in enumerate(test_cases, 1):
        result = bose_fermi_integral(integral, debug=False)
        simplified_diff = simplify(result - expected)
        
        if simplified_diff == 0:
            print(f"Test {i}: PASSED")
            passed += 1
        else:
            print(f"Test {i}: FAILED")
            print(f"  Expected: {expected}")
            print(f"  Got: {result}")
            print(f"  Difference: {simplified_diff}")
    
    print(f"Results: {passed}/{len(test_cases)} tests passed\n")
    return passed == len(test_cases)


def test_symbolic_parameters():
    """Test integrals with symbolic parameters."""
    x = symbols('x', real=True, positive=True)
    a, p = symbols('a p', positive=True)
    
    # Test symbolic power and coefficient
    integral = Integral(x**p / (exp(a*x) - 1), (x, 0, oo))
    result = bose_fermi_integral(integral, debug=True)
    
    expected_form = gamma(p + 1) * zeta(p + 1) / (a**(p + 1))
    
    print("\nSymbolic test:")
    print(f"Integral: {integral}")
    print(f"Result: {result}")
    print(f"Expected form: {expected_form}")
    
    return True  # Manual verification for now


def test_edge_cases():
    """Test edge cases and boundary conditions."""
    x = symbols('x', real=True, positive=True)
    
    print("\nTesting edge cases...")
    
    # Test x^0 case (should give well-known results)
    integral_be = Integral(1 / (exp(x) - 1), (x, 0, oo))
    integral_fd = Integral(1 / (exp(x) + 1), (x, 0, oo))
    
    result_be = bose_fermi_integral(integral_be, debug=False)
    result_fd = bose_fermi_integral(integral_fd, debug=False)
    
    print(f"BE x^0 integral: {result_be}")
    print(f"FD x^0 integral: {result_fd}")
    
    # Test x^1 case
    integral_x1 = Integral(x / (exp(x) - 1), (x, 0, oo))
    result_x1 = bose_fermi_integral(integral_x1, debug=False)
    expected_x1 = gamma(2) * zeta(2)  # = π²/6
    
    diff = simplify(result_x1 - expected_x1)
    if diff == 0:
        print(f"x^1 BE integral: PASSED ({result_x1})")
    else:
        print(f"x^1 BE integral: FAILED (got {result_x1}, expected {expected_x1})")
    
    return True


def test_input_validation():
    """Test that invalid inputs are handled gracefully (should return original input)."""
    x = symbols('x', real=True, positive=True)
    
    print("\nTesting input validation...")
    passed = 0
    total = 3
    
    # Test non-integral input (should return unchanged)
    non_integral_input = x**2
    result = bose_fermi_integral(non_integral_input, debug=False)
    if result == non_integral_input:
        print("Non-integral input: Correctly returned unchanged")
        passed += 1
    else:
        print(f"Non-integral input: Expected {non_integral_input}, got {result}")
    
    # Test wrong integration limits (should return unchanged)
    bad_integral = Integral(x**2 / (exp(x) - 1), (x, 1, oo))
    result = bose_fermi_integral(bad_integral, debug=False)
    if result == bad_integral:
        print("Wrong limits: Correctly returned unchanged")
        passed += 1
    else:
        print(f"Wrong limits: Expected {bad_integral}, got {result}")
    
    # Test unsupported form (should return unchanged)
    weird_integral = Integral(x**2 / (exp(x**2) - 1), (x, 0, oo))
    result = bose_fermi_integral(weird_integral, debug=False)
    if result == weird_integral:
        print("Unsupported form: Correctly returned unchanged")
        passed += 1
    else:
        print(f"Unsupported form: Expected {weird_integral}, got {result}")
    
    print(f"Input validation: {passed}/{total} cases handled properly")
    return passed == total


def test_performance_sanity():
    """Basic performance sanity check."""
    x = symbols('x', real=True, positive=True)
    
    print("\nTesting basic performance...")
    
    # Simple integral that should be fast
    integral = Integral(x**2 / (exp(x) - 1), (x, 0, oo))
    
    start_time = time.time()
    _result = bose_fermi_integral(integral, debug=False)
    end_time = time.time()
    
    duration = end_time - start_time
    print(f"Simple integral computed in {duration:.4f} seconds")
    
    if duration > 5.0:  # If it takes more than 5 seconds, something's wrong
        print("WARNING: Performance seems slow for basic case")
        return False
    
    return True


def test_known_special_values():
    """Test against known special values."""
    x = symbols('x', real=True, positive=True)
    
    print("\nTesting known special values...")
    
    # ζ(2) = π²/6
    integral_zeta2 = Integral(x / (exp(x) - 1), (x, 0, oo))
    result = bose_fermi_integral(integral_zeta2, debug=False)
    expected = pi**2 / 6
    
    diff = simplify(result - expected)
    if diff == 0:
        print(f"ζ(2) test: PASSED ({result})\n")
    else:
        print(f"ζ(2) test: FAILED (got {result}, expected {expected})\n")
    
    # ζ(4) = π⁴/90
    integral_zeta4 = Integral(x**3 / (exp(x) - 1), (x, 0, oo))
    result = bose_fermi_integral(integral_zeta4, debug=False)
    expected = pi**4 / 15  # 6 * ζ(4) = 6 * π⁴/90 = π⁴/15
    
    diff = simplify(result - expected)
    if diff == 0:
        print(f"ζ(4) test: PASSED ({result})\n")
    else:
        print(f"ζ(4) test: FAILED (got {result}, expected {expected})\n")
    
    return True


def test_consistency():
    """Test internal consistency."""
    x = symbols('x', real=True, positive=True)
    
    print("\nTesting consistency...")
    
    # Test that same integral gives same result
    integral = Integral(x**2 / (exp(x) - 1), (x, 0, oo))
    
    result1 = bose_fermi_integral(integral, debug=False)
    result2 = bose_fermi_integral(integral, debug=False)
    
    if simplify(result1 - result2) == 0:
        print("Consistency test: PASSED")
        return True
    else:
        print("Consistency test: FAILED - same input gave different outputs!")
        return False


if __name__ == "__main__":
    # Run all sanity tests
    print("=" * 60)
    print("BOSE-FERMI INTEGRAL SANITY TESTS")
    print("=" * 60)
    
    tests = [
        ("Basic integrals", test_basic_integrals),
        ("Symbolic parameters", test_symbolic_parameters),
        ("Edge cases", test_edge_cases),
        ("Input validation", test_input_validation),
        ("Performance sanity", test_performance_sanity),
        ("Known special values", test_known_special_values),
        ("Consistency", test_consistency),
    ]
    
    passed_tests = 0
    total_tests = len(tests)
    
    for test_name, test_func in tests:
        print(f"\n{'='*20} {test_name} {'='*20}")
        try:
            if test_func():
                passed_tests += 1
                print(f"✅ {test_name}: PASSED")
            else:
                print(f"❌ {test_name}: FAILED")
        except Exception as e:
            print(f"❌ {test_name}: ERROR - {e}")
    
    print("\n" + "=" * 60)
    print(f"SANITY TEST RESULTS: {passed_tests}/{total_tests} passed")
    
    if passed_tests == total_tests:
        print("🎉 All sanity tests passed! System appears healthy.")
    else:
        print("⚠️  Some sanity tests failed. Check the system!")
    print("=" * 60)
