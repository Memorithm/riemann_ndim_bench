# Finite Weil quadrature refinement

Status: **EMPIRICAL RESOLUTION DIAGNOSTIC / NOT CERTIFIED ERROR BOUNDS**

This layer re-evaluates one fixed compact support and one fixed basis dimension at a declared sequence of quadrature orders.

For each level it records the raw and Gram-normalized minimum eigenvalues, Gram condition number, critical-boundary residual, mixed-pairing asymmetry and whitening asymmetry.

The audit additionally reports the minimum, maximum and span of the values observed across the supplied levels, plus the absolute delta between the last two levels.

These quantities are **not** rigorous error bars, confidence intervals, or certified enclosures. They are empirical refinement diagnostics only. A small last-step delta can be evidence of numerical stabilization for the tested sequence, but by itself does not bound the remaining quadrature error.

## Release probe

```text
cargo run --release --example weil_refinement -- 4
```

The default probe uses support `(1/2, 7/2)`, dimension `N=4`, and the levels:

```text
(48, 48, 64, 64)
(72, 72, 96, 96)
(96, 96, 128, 128)
```

where each tuple is

```text
(correlation_order, archimedean_order, boundary_order, gram_order).
```

The example prints every sample and the observed spans/deltas.

## Interpretation

A candidate spectral feature should be compared simultaneously against:

- the observed refinement span;
- the last-step refinement delta;
- Gram conditioning;
- boundary residuals;
- pairing symmetry residuals;
- whitening symmetry residuals;
- changes under basis dimension and support window.

No sign is assumed by the API or tests.

## Scientific boundary

Refinement stability at finitely many orders does not certify the infinite-order integral, establish positivity on the complete admissible test-function space, prove density/completeness, justify the semilocal trace formula, establish Conjecture 4.1, or prove RH.
