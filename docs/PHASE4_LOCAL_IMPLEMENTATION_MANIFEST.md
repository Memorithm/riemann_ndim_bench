# Phase 4 local Rust implementation manifest

This file records the implementation that has been validated on Thor but has **not yet been upstreamed as source code**. It exists to prevent the research state from depending on one working tree.

Last consolidated: 2026-08-15.

## Local branch

The Phase 4 exploratory work has been performed locally in

```text
/root/riemann_ndim_bench
```

on branch

```text
phase4-pi-deformation
```

whose original base was the Phase 3 merge commit

```text
24e2743
```

The branch contains extensive uncommitted exploratory tests. It must not be reconstructed destructively from this manifest; the preferred upstream path is to inspect the Thor working tree and import validated pieces into focused modules/tests.

## Known local exploratory files

```text
tests/pi_radial_probe.rs
tests/pi_deformed_spectrum.rs
tests/radial_base_control.rs
tests/archimedean_radial_probe.rs
tests/archimedean_balance_control.rs
tests/semilocal_moment_probe.rs
tests/semilocal_p2_measure_probe.rs
tests/semilocal_p2_jacobi_probe.rs
tests/semilocal_p2_stieltjes_probe.rs
tests/semilocal_multi_prime_probe.rs
tests/semilocal_first_order_growth_fit.rs
```

The nested local `RiemannBench/` directory is an independent Git repository and is unrelated to the Phase 4 import. It must not be added or deleted as part of the upstreaming work.

## Negative controls already resolved

### Experimental pi deformation

The deformation

```text
u = (q/ln 2)^2
x_eff = x + u*x*(1-x/ln 2)
```

passed basic symmetry checks but failed specificity controls. Replacing `pi` by arbitrary bases showed that the effect was only a scalar deformation strength. This line is negative evidence and should not be revived as evidence for a distinguished role of `pi`.

### Archimedean balance alone

For

```text
A(s) = pi^(-s/2) Gamma(s/2),
q_inf(sigma,t)
 = (1/2) [log|A(s)| - log|A(1-s)|],
```

numerical evaluation at zeta zeros reproduced the functional-equation identity

```text
q_inf
 = (1/2) [log|zeta(1-s)| - log|zeta(s)|].
```

Thus this balance alone cannot characterize the critical line independently of zeta.

## Semilocal measure and Jacobi validation

The local implementation uses the semilocal weight

```text
|Gamma(1/4 + i t/2)|^2
* product_p 1/[1 - 2 p^(-1/2) cos(t log p) + p^(-1)].
```

The `p=2` moments, Stieltjes coefficients and low-dimensional Jacobi spectra were reproduced against the published rounded benchmarks.

Representative `p=2` Jacobi coefficients are

```text
a0 = 0.396868411277
a1 = 1.061680386162
a2 = 2.620426310397
a3 = 5.134945422201
a4 = 5.680705489972
a5 = 4.802674782590
a6 = 5.411988157671
a7 = 7.902019936502
```

The exact archimedean coefficients are

```text
a_n = sqrt((n+1/2)(n+1)).
```

## Multi-prime exploratory coefficients

Fine Stieltjes runs produced, for example,

```text
S={infinity,2,3}
0.280245981728
0.716587591679
2.765334984800
4.890780885101
4.489266293762
6.440770931015
6.189382751435
7.839356341293
```

and

```text
S={infinity,2,3,5}
0.225637784741
0.540678527503
2.650602730548
4.896849221684
4.801433880964
6.463423572310
6.067969664663
8.148271006043
```

These are exploratory finite-S numerical outputs; no direct published benchmark for these specific multi-prime sets has been claimed.

## Generalized prolate crossing implementation

The local implementation contains

```text
prolate_block_coefficients(
    jacobi,
    lambda,
    even_block,
    block_size,
)
```

using parity degree

```text
d = 2n     for W+
d = 2n+1   for W-.
```

The generalized crossing problem is transformed to

```text
K = B^(-1/2) (-A) B^(-1/2),
B_d = 2*pi*(4d+1).
```

Positive eigenvalues `mu` give

```text
lambda = sqrt(mu).
```

The generalized-EVD crossings were validated against direct bisection for block sizes `4, 6, 8`, both parities and several finite prime sets.

## Exact first-order q derivative

Known local functions include

```text
semilocal_alpha
archimedean_jacobi_squared
semilocal_jacobi_squared_q_derivative
semilocal_jacobi_relative_q_derivative
exact_first_order_prolate_coefficients_from_source
exact_first_order_prolate_derivative_closed_form
exact_first_order_prolate_crossing_derivatives
merged_exact_first_order_prolate_response
```

The source and closed-form derivatives agree to machine precision.

Validated tests include

```text
exact_first_order_q_closed_form_matches_source
exact_first_order_q_sign_structure_holds
exact_first_order_q_regression_matches_large_p_limit
```

The sign-corrected derivative was also checked to remain strictly diagonally dominant through block size `1024`.

## Tridiagonal EVD

The optimized path uses faer 0.24.4 tridiagonal self-adjoint EVD. The local helper

```text
exact_first_order_prolate_crossing_derivatives_tridiagonal(
    even_block,
    block_size,
)
```

was validated against the dense eigenvector solver for

```text
m = 8, 16, 32, 64, 128
```

with zero observed difference in the stored crossing values and derivatives.

## High-dimensional total-variation methods

### Global eigenvalues-only trace method

The first optimization computed

```text
Tr sqrt(K(+h)) - Tr sqrt(K(-h))
```

using tridiagonal eigenvalues only. It was accurate through moderate sizes but developed cancellation at `m=16384` when `h` became too small.

### Pairwise eigenvalue method

The safer local implementation computes each

```text
(lambda_j(+h) - lambda_j(-h)) / (2h)
```

before summation. Known helpers include

```text
first_order_sqrt_spectrum_tridiagonal_eigenvalues_only
pairwise_first_order_trace_derivative
pairwise_total_variation_derivative
```

This is the accepted high-dimensional method at present.

## Accepted total-variation checkpoints

```text
S(128)   = 3.970845543531
S(256)   = 4.640481894221
S(512)   = 5.359223651882
S(1024)  = 6.126883687871
S(2048)  = 6.943355708182
S(4096)  = 7.808580171428
S(8192)  = 8.72252476599664
S(16384) = 9.68517075774107
```

For `m=8192` and `m=16384`, independent interior-step tests validate the quadratic-in-`h^2` extrapolation. See `PHASE4_NUMERICAL_CHECKPOINTS_2026-08-15.md` for the detailed cancellation diagnostics and rejected extrapolants.

## Growth-fit test file

`tests/semilocal_first_order_growth_fit.rs` contains local Rust tests for

- quadratic-log growth versus a free `(log m)^p` fit;
- finite-size centered coefficients `A_m`;
- equal-parameter comparison between `(c log m+d)/m` and `c/m+d/m^2` corrections.

The five-point pre-`m=16384` fit favored the log-over-m correction model by a factor about `4.34` in RMSE. The homogeneous `m=8192` and `m=16384` pairwise checkpoints now require this fit to be rerun with the corrected six-point `A_m` sequence before any updated fitted limit is promoted.

## Current upstreaming priority

Do **not** upstream the exploratory file wholesale without review. Extract in this order:

1. exact `q=0` source/closed-form derivative helpers;
2. finite sign-lemma regression tests;
3. tridiagonal exact first-order EVD helper;
4. eigenvalues-only pairwise finite-difference helper;
5. compact high-dimensional regression checkpoints;
6. growth-analysis tests with only accepted homogeneous data.

Keep negative-control probes and exploratory multi-prime diagnostics separate from the stable library API.
