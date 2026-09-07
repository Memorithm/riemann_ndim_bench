# Finite Weil sine pairing refinement grid

This audit crosses two independent finite numerical axes for the non-polynomial sine target family:

1. shifted-Legendre parent dimension, which controls the finite approximation space;
2. coefficient and Riemann--Weil quadrature orders, which control numerical resolution.

For every parent dimension and quadrature level, the underlying direct-pairing audit compares

\[
\psi\!\left((h_i^{\mathrm{direct}})^* * h_j^{\mathrm{direct}}\right)
\]

with the pairing of the corresponding finite Legendre projection. It also checks the projected generic route against the historical parent matrix through `C^T A C`.

The refinement grid retains the complete per-level diagnostics and, for each parent dimension, reports:

- the observed interval and span of the maximum direct-vs-projected pairing residual;
- the observed interval and span of the normalized pairing residual;
- consecutive last-step changes for both residuals;
- the largest projected-generic versus parent-matrix disagreement;
- the largest forward/reverse pairing asymmetry;
- the largest critical Mellin boundary residual.

The per-level rows additionally retain the separate pole, archimedean and prime-power residuals from the full pairing decomposition.

## Interpretation

The two axes should not be conflated. A pairing residual that changes materially when quadrature is refined is not yet clean evidence about the Legendre truncation. Conversely, stability in quadrature while the residual changes with parent dimension is evidence about the finite approximation sequence rather than the integration resolution.

Observed spans and consecutive deltas are empirical diagnostics. They are not confidence intervals, certified error bounds, or a proof of convergence. No monotonicity or sign condition is imposed by the tests.

Even simultaneous stability in parent dimension and quadrature does not establish continuity of the Weil functional in a topology sufficient for the Weil criterion, density/completeness of the finite family, complete-space Weil positivity, Conjecture 4.1, or the Riemann hypothesis.

A later mathematical step must identify a topology/norm in which the admissible finite approximants are dense and in which the Weil functional is continuous or otherwise uniformly controlled. The present grid is intended to expose numerical behavior relevant to that theorem target, not replace it.
