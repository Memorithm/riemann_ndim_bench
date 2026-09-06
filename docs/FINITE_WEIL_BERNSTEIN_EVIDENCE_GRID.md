# Finite Weil Bernstein evidence grid

This experiment crosses the localized Bernstein-subspace axis with two numerical controls that were already audited independently:

- exact rational compact-support windows;
- declared quadrature-refinement levels.

For one fixed Bernstein index subspace, every `(support, quadrature level)` sample is computed by `audit_finite_weil_bernstein_subspace`. That routine builds one parent Legendre pairing/Gram pair `(A,G)`, transforms it to the selected Bernstein subspace, and extracts a dimension-matched leading-Legendre control from the same parent computation.

The primary finite quantities are therefore

```text
lambda_B = minimum generalized eigenvalue on the selected Bernstein subspace
lambda_L = minimum generalized eigenvalue on the leading Legendre control
delta    = lambda_B - lambda_L
```

The grid keeps `lambda_B`, `lambda_L`, and `delta` separate. It also records the Bernstein and control Gram condition numbers, the selected boundary residual, parent directional-pairing asymmetry, and Bernstein whitening asymmetry.

Across refinement levels, each support cell reports observed minima/maxima, spans, and last-step absolute deltas. These are empirical resolution diagnostics only. They are not confidence intervals, certified numerical error bounds, significance thresholds, or proof scores.

## Interpretation

This grid asks a narrow question: does the finite generalized spectral behavior attributed to the selected localized Bernstein subspace remain qualitatively stable when support and quadrature resolution are changed, relative to a same-dimension leading-Legendre control?

A stable family difference is stronger finite evidence that the observed behavior is not merely a single-resolution artifact. It is still a statement about declared finite subspaces.

## Scientific boundary

Even if every tested Bernstein and Legendre generalized minimum is positive and stable across all declared windows and refinement levels, this does **not** establish:

- positivity of the Weil quadratic form on the complete admissible function space;
- density or completeness of the tested families;
- a uniform approximation or convergence theorem;
- a certified numerical error enclosure for the infinite problem;
- identification with the semilocal Hilbert space `L^2(X_S)`;
- Conjecture 4.1;
- the Riemann hypothesis.

Likewise, a negative finite cell would be an experimentally important direction, but it would require independent replication, stronger numerical error control, and verification that conditioning and discretization are not responsible before any mathematical interpretation.

## Next research step

The Bernstein family used here still lives inside a finite polynomial parent space. The next family-level step should introduce compact generators outside that same finite polynomial span while preserving the existing boundary, Gram, support, and refinement diagnostics.
