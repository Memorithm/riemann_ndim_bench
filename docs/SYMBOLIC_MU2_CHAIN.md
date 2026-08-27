# Exact post-perturbative mu² symbolic verifier

This checkpoint adds a deterministic symbolic layer after the exact second-order recurrence extraction already provided by `verify_math.py`.

It does **not** encode the final limiting constant and does not claim any implication for the Riemann hypothesis.

## Purpose

The blind dual-agent benchmark repeatedly reached the homogeneous irregular mode and Gamma tail constants but failed to derive the inhomogeneous µ² chain that connects the forcing to a generating function and then to a finite part.

`tools/symbolic_mu2.py` provides two exact, restricted operations for that post-perturbative stage.

## Hypergeometric / Pochhammer quotient

The `hypergeometric` mode takes rational numerator/denominator shifts and reconstructs the exact coefficient quotient

```text
t[k+1] / t[k].
```

If the caller supplies a candidate quotient, the result is fail-closed:

- exact equality -> `candidate_status=PROVED_EQUAL` and `exact_status=PROVED_BY_POCHHAMMER_QUOTIENT`;
- mismatch -> `candidate_status=MISMATCH` and `exact_status=REFUTED_CANDIDATE_RATIO`.

The µ² forcing subsequence

```text
a_k = binom(4k, 2k) / 16^k
```

is represented exactly as

```text
a_k = (1/4)_k (3/4)_k / ((1/2)_k (1)_k),
```

with quotient

```text
a_{k+1}/a_k
  = ((k+1/4)(k+3/4))/((k+1/2)(k+1)).
```

The regression suite checks the exact first terms

```text
1, 3/8, 35/128, 231/1024, 6435/32768, 46189/262144.
```

## Theta operators and Puiseux finite parts

The `finite-part` mode accepts:

- an exact algebraic generating expression in `z`;
- a polynomial in the sequence index `k`, interpreted as the same polynomial in `theta = z d/dz`;
- an optional exact algebraic extra term.

It then substitutes the fixed local coordinate

```text
x = sqrt(1-z),
z = 1-x^2,
```

computes an exact Puiseux/Laurent expansion and extracts the coefficient of `x^0`.

For

```text
A(z) = 1/2 * [
  (1-sqrt(z))^(-1/2)
  + (1+sqrt(z))^(-1/2)
],
```

the suite verifies exactly

```text
T_+(z) = (16/9) (8 theta^2 + theta) A(z)
```

has finite part

```text
c_+ = -sqrt(2)/9,
```

and

```text
B(z) = [sqrt(1+sqrt(z)) - sqrt(1-sqrt(z))]/sqrt(z),
```

```text
T_-(z)
  = (16/45) (16 theta^2 + 6 theta - 1/2) A(z)
    + (8/45) B(z)
```

has finite part

```text
c_- = 2 sqrt(2)/45.
```

These are exact symbolic statements about the supplied generating expressions. The verifier does not by itself prove that a caller selected the correct generating expression from the source recurrence; that derivation remains a separate evidence step and must not be silently skipped.

## Safety boundary

Model-provided expressions are parsed through a small Python AST whitelist. Arbitrary Python execution, attribute access and general SymPy parsing of untrusted text are not allowed. The only permitted function call in algebraic finite-part expressions is `sqrt(...)`.

## CI

The research CI installs pinned `sympy`, compiles the symbolic helper and runs `tools/test_symbolic_mu2.py` alongside the existing verifier/ledger suite.
