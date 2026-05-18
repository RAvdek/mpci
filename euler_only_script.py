"""
compute_M2n.py
==============

Computes M_{2n} = GCD of c_{(n)}(X) over X in T_{2n} for n = 2, ..., N_MAX.

Also tests the conjectured formulas:
  - Even n: M_{2n} = denom(B_n / n!) where B_n is the standard Bernoulli number
  - Odd n >= 3: M_{2n} = (n-1)!
  - n = 1: M_2 = 2

Usage:
    python compute_M2n.py

Adjust N_MAX below to control how far the computation goes.
The computation for n=9 may take a few minutes; n=10 could take longer.
"""

import sys
import time
import logging
import sympy

# Use the integrated file (rename mpci_final.py to mpci.py or adjust import)
import mpci

# Adjust this to control how far to compute
N_MAX = 10

LOG_FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)

results = {}

print(f"Computing M_{{2n}} for n = 2 to {N_MAX}")
print("=" * 70)

for n in range(2, N_MAX + 1):
    t0 = time.time()
    val = mpci.get_euler_only_v2(n, log=True)
    elapsed = time.time() - t0
    results[n] = val
    print(f"  n={n:2d}: M_{{2n}} = {val}, time = {elapsed:.1f}s")
    print()

# Also record M_2 = 2 (n=1 case, trivial)
results[1] = 2

print()
print("=" * 70)
print("RESULTS AND CONJECTURE TEST")
print("=" * 70)
print()
print(f"{'n':>3} {'M_{2n}':>15} {'even formula':>15} {'odd formula':>15} {'match':>6}")
print("-" * 60)

for n in sorted(results.keys()):
    val = results[n]
    bn = sympy.bernoulli(n)

    if n == 1:
        # Special case
        even_pred = ""
        odd_pred = ""
        note = "(special)"
    elif n % 2 == 0:
        # Even: denom(B_n / n!)
        d = sympy.Rational(bn, sympy.factorial(n)).denominator
        even_pred = str(d)
        odd_pred = ""
        note = "✓" if d == val else "✗"
    else:
        # Odd: (n-1)!
        f = int(sympy.factorial(n - 1))
        even_pred = ""
        odd_pred = str(f)
        note = "✓" if f == val else "✗"

    print(f"{n:>3} {str(val):>15} {even_pred:>15} {odd_pred:>15} {note:>6}")

print()
print("Even formula: M_{2n} = denom(B_n / n!)")
print("Odd formula:  M_{2n} = (n-1)!")
print()

# Also print the polynomial generators used
print("=" * 70)
print("POLYNOMIAL GENERATORS")
print("=" * 70)
for n in range(1, N_MAX + 1):
    gen = mpci.get_polynomial_generator(n)
    s_val = sum(c * m.get_s_number() for c, m in gen)
    is_pp, p, k = mpci._is_prime_power(n + 1)
    if is_pp:
        target = f"±{p}"
    else:
        target = "±1"
    terms = " + ".join(f"{c}*[{m}]" for c, m in gen)
    print(f"  x_{n:2d} = {terms}")
    print(f"         s_{n} = {int(s_val)} (target {target}), #terms = {len(gen)}")