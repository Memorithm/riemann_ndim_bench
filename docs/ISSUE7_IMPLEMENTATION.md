# Issue #7 implementation checkpoint

This branch upstreams the exact finite `q=0` single-prime semilocal prolate derivative that had previously been validated only in the Thor exploratory workspace.

## Implemented source-derived formulas

The Rust module `src/semilocal.rs` implements

```text
alpha_0 = 1,
alpha_{n+1} = -((n+1/2)/(n+1)) alpha_n,
```

```text
a_n(0)^2 = (n+1/2)(n+1),
```

and

```text
(a_n^2)'(0)
  = -(1/sqrt(2))(2n+1)(4n+3) alpha_n.
```

It builds both `K(0)` and `K'(0)` for `W+` and `W-`. `K'(0)` is implemented twice: once from the unsimplified source-derived coefficient and once from the closed forms in `PHASE4_FIRST_ORDER_DERIVATION.md`, so the simplification can be regression-tested independently.

The crossing derivatives use `faer::linalg::solvers::SelfAdjointEigen` and the documented `U()` / `S()` decomposition factors:

```text
mu'_j(0) = u_j^T K'(0) u_j,
lambda'_j(0) = mu'_j(0)/(2 lambda_j).
```

## Regression coverage

`tests/semilocal_exact_derivative.rs` checks:

- exact initial terms of the stable `alpha_n` recurrence;
- the closed `(a_n^2)'(0)` formula against the unsimplified first-order coefficient;
- closed-form `K'(0)` against its independent source-derived construction;
- strict diagonal dominance of the sign-corrected perturbation through block size `1024`;
- negative `W+` and positive `W-` Rayleigh derivatives on normal test sizes;
- the independent large-prime aggregate targets at `m=16,24,32`;
- the documented `m=128` shape statistics;
- an ignored expensive regression reproducing total response through `m=1024`.

The expensive test is intentionally ignored in the default suite because this implementation uses the dense `SelfAdjointEigen` path required by issue #7. It remains available explicitly before promoting high-dimensional asymptotics.

## Independent reconstruction check

The formulas in this branch were independently reconstructed outside the Rust implementation. They reproduce, to floating-point precision, the historical targets recorded in the issue and Phase-4 notes, including

```text
m=16: m*mean_abs         = 0.28263208002937096
      m*trimmed_mean_abs = 0.22563342904799336

m=24: m*mean_abs         = 0.26113744999423477
      m*trimmed_mean_abs = 0.18744825444964030

m=32: m*mean_abs         = 0.24575569288899024
      m*trimmed_mean_abs = 0.16048091598744970

S(128) = 3.9708455435305754
S(256) = 4.640481894221456
```

## Validation boundary

GitHub Actions currently fails before repository steps are scheduled because of the pre-existing runner/billing state (`runner_id=0`, `steps=[]`). Thor is also currently unavailable. The branch has therefore been reviewed against the documented `faer 0.24.4` API and independently numerically reconstructed, but the repository Rust CI still needs to execute once a runner is available.

The issue branch has been refreshed onto the `main` state that already contains PR #10, so it is intended to merge without dropping the verified research harness.

This remains a finite-compression perturbation calculation. It does not identify these crossings with zeta zeros and does not imply the Riemann hypothesis.
