# Finite Weil sine truncation audit

This experiment introduces a non-polynomial **target** generator family without introducing a second implementation of the Riemann--Weil functional.

For the existing compact support coordinate

```text
t = (rho - a) / (b - a),  0 < t < 1,
```

the target generators are

```text
g_m(rho) = bump(rho) sin(m pi t),  m = 1, 2, ...
```

and the intended boundary-admissible directions are `h_m = Q g_m`, with the same differential operator `Q` used by the existing compact Weil experiments.

## Why a truncation audit

The current validated finite Riemann--Weil implementation evaluates pairings on the compact Legendre family. Rather than duplicating that source-sensitive implementation for sine functions, the sine enrichment is projected onto shifted Legendre polynomials:

```text
sin(m pi t) ~= sum_{j=0}^{N-1} c_{m,j} P_j(2t-1).
```

The coefficients are computed by Gauss--Legendre quadrature,

```text
c_{m,j} = (2j+1) integral_0^1 sin(m pi t) P_j(2t-1) dt.
```

For each declared parent dimension `N`, the coefficient vectors define a finite subspace in the already-computed parent pairing and Gram matrices. The transformed problem is evaluated as

```text
A_C = C^T A C,
G_C = C^T G C.
```

A dimension-matched leading-Legendre generalized spectrum is extracted from the same parent computation as a control.

## Recorded diagnostics

For every parent dimension the audit records:

- maximum sampled reconstruction residual of the sine enrichments;
- maximum coefficient L1 norm;
- raw minimum eigenvalue;
- Gram-normalized generalized minimum eigenvalue;
- same-dimension leading-Legendre generalized minimum;
- signed family delta between those two generalized minima;
- both Gram condition numbers;
- reconstructed critical boundary residual;
- parent directional pairing asymmetry;
- transformed whitening asymmetry.

Across parent dimensions it also reports observed spectral spans and last-step changes.

## Numerical boundary

The reconstruction residual is the maximum difference sampled on a fixed deterministic grid in the affine coordinate. It is **not** a certified uniform approximation bound.

It measures reconstruction of the enrichment `sin(m pi t)`. It is not a direct certified error estimate for `Q g_m`, for a Riemann--Weil pairing, or for the generalized eigenvalue.

Likewise, stabilization of the generalized minimum as the parent dimension increases is an empirical truncation diagnostic, not a proof of convergence of the infinite Legendre series under the Weil functional.

## Scientific boundary

At every finite parent dimension the computed sine approximation is still an element of the finite polynomial Legendre parent space. Therefore even a stable positive sequence does **not** establish:

- exact evaluation of the non-polynomial sine family under the Weil functional;
- convergence of the transformed quadratic form to a limiting sine-family form;
- density or completeness of the chosen test-function families;
- a uniform approximation theorem in the topology required by the Weil criterion;
- positivity of the Weil quadratic form on the complete admissible space;
- identification with semilocal `L^2(X_S)`;
- Conjecture 4.1;
- the Riemann hypothesis.

A stable negative finite direction would be experimentally important, but it would still require independent reproduction and stronger control of projection, quadrature, conditioning, and the action of `Q` before mathematical interpretation.

## Next research step

If the reconstruction and finite spectra stabilize coherently as the parent dimension grows, the next step is to derive or numerically audit the action of `Q` directly on the sine-enriched generators and compare that independent route against the Legendre-projected result. That comparison should precede any claim that the finite projection sequence represents the non-polynomial target family faithfully under the Weil functional.
