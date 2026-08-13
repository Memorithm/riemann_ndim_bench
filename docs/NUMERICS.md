# Numerical core: phase 1

## Purpose

This module provides an independent numerical layer for the research bench. It is deliberately separate from the experimental radial coordinate.

The goal is not to prove the Riemann hypothesis. The goal is to have reproducible numerical checks that future N-dimensional constructions must satisfy.

## Riemann zeta evaluation

The implementation uses an Euler-Maclaurin continuation of the zeta function. In the form used by the code,

```text
zeta(s) ~= sum_{n=1}^{N-1} n^(-s)
          + N^(1-s)/(s-1)
          + (1/2) N^(-s)
          + sum_{k=1}^{m} [B_(2k)/(2k)!] (s)_(2k-1) N^(1-s-2k).
```

This is the standard Bernoulli correction structure of the Euler-Maclaurin representation for the Hurwitz zeta function specialized to the Riemann zeta function.

Primary reference:

- NIST DLMF §25.11(iii), Euler-Maclaurin representations: https://dlmf.nist.gov/25.11.iii

The current `f64` implementation is a research probe, not an interval-certified evaluator. The code therefore compares two truncation resolutions and reports their difference as an **empirical cross-resolution discrepancy**. That difference is not claimed to be a rigorous error bound.

## Completed zeta function

The implementation forms

```text
xi(s) = (1/2) s (s-1) Gamma(s/2) pi^(-s/2) zeta(s)
```

and checks the exact identity

```text
xi(s) = xi(1-s).
```

Primary reference:

- NIST DLMF §25.4, equations 25.4.3 and 25.4.4: https://dlmf.nist.gov/25.4

The Gamma function is currently evaluated using a double-precision Lanczos approximation with the reflection formula. This is sufficient for the present regression checks but will not be treated as a certified high-precision backend.

## Independent regression checks

The test suite includes:

- `zeta(2) = pi^2/6`, an exact special value documented by NIST DLMF §25.6;
- `zeta(0) = -1/2`, an exact special value documented by NIST DLMF;
- `Gamma(1) = 1`;
- `Gamma(1/2) = sqrt(pi)`;
- conjugation symmetry of the numerical zeta evaluator;
- functional symmetry of `xi`;
- convergence consistency between two Euler-Maclaurin truncation resolutions.

Primary references:

- NIST DLMF §25.6: https://dlmf.nist.gov/25.6
- NIST DLMF §25.4: https://dlmf.nist.gov/25.4

## Current limitations

1. Arithmetic is IEEE-754 `f64` only.
2. The reported cross-resolution difference is not a proof of accuracy.
3. The Euler-Maclaurin parameter choice is conservative and not performance-tuned.
4. The current code is not intended for very large imaginary parts.
5. No conclusion about zero locations is drawn from this evaluator alone.

A later validation layer should add arbitrary precision and independent reference values before the N-dimensional spectral experiments are allowed to use numerically delicate conclusions.
