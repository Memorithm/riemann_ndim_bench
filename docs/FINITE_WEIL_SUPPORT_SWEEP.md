# Finite Weil support-window sweep

Status: **FINITE PARAMETER SWEEP / NOT WEIL POSITIVITY**

This increment extends the Gram-normalized dimension sweep by varying the exact compact support window as an independent numerical parameter.

## Construction

For every declared rational open support `(a,b)`, the implementation constructs the existing compact bump, builds one maximum-dimensional finite Weil pairing matrix `A` and one multiplicative Gram matrix `G`, then extracts all leading-principal rows `N=1..N_max` without recomputing the expensive pairings inside that window.

Every row reports:

- `lambda_min(A_N)` in the raw coefficient basis;
- the generalized minimum eigenvalue of `A_N v = lambda G_N v`;
- minimum and maximum Gram eigenvalues;
- the Gram condition number;
- the maximum critical-boundary residual for the full window computation;
- the raw mixed-pairing symmetry residual;
- the whitening asymmetry residual.

The support endpoints remain exact positive rationals until the existing numerical bump/quadrature layer is entered.

## Release probe

Run the default manufactured width sweep with:

```text
cargo run --release --example weil_support_sweep -- 4 96 96 128 128
```

The positional parameters are:

```text
max_dimension correlation_order archimedean_order boundary_order gram_order
```

The default windows are:

```text
(3/4, 13/4)
(1/2, 7/2)
(1/4, 15/4)
```

They have the same arithmetic midpoint `rho=2` and progressively larger widths. They are manufactured numerical audit parameters only. No theoretical significance is attached to this center or to these widths.

The example emits CSV-compatible rows with exact endpoint numerators/denominators and all spectral/conditioning diagnostics.

## Interpretation

Changing the support can change both the finite Weil pairing and the conditioning of the chosen basis. Therefore an apparent movement of the generalized minimum eigenvalue must be read together with:

- dimension dependence;
- Gram conditioning;
- critical-boundary residuals;
- pairing and whitening symmetry residuals;
- quadrature refinement.

No sign is asserted by the support-sweep API or its structural tests.

## Scientific boundary

A positive generalized minimum on every tested dimension and every tested support window is still only finite numerical evidence for those declared spans and resolutions. It does not establish:

- positivity for all admissible compact test functions;
- completeness or density of the Legendre-bump family;
- uniform control in dimension or support size;
- identification of the numerical `L^2(d^*rho)` Gram norm with the semilocal `L^2(X_S)` Hilbert space;
- Conjecture 4.1;
- the Riemann hypothesis.

Likewise, a stable negative row would be a candidate finite negative direction requiring independent replication, quadrature refinement, and conditioning analysis before interpretation.
