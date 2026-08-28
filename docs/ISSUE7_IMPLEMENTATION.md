# Issue #7 implementation checkpoint

Issue #7 upstreamed the exact finite `q=0` single-prime semilocal prolate derivative that had previously been validated only in the Thor exploratory workspace.

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

The dense reference crossing derivatives use `faer::linalg::solvers::SelfAdjointEigen` and the documented `U()` / `S()` decomposition factors:

```text
mu'_j(0) = u_j^T K'(0) u_j,
lambda'_j(0) = mu'_j(0)/(2 lambda_j).
```

The later `src/semilocal_tridiagonal.rs` path feeds the same exact `K(0)` directly to faer's self-adjoint tridiagonal EVD. It preserves the dense `SelfAdjointEigen` implementation as an independent differential oracle rather than replacing it.

## Regression coverage

`tests/semilocal_exact_derivative.rs` checks:

- exact initial terms of the stable `alpha_n` recurrence;
- the closed `(a_n^2)'(0)` formula against the unsimplified first-order coefficient;
- closed-form `K'(0)` against its independent source-derived construction;
- strict diagonal dominance of the sign-corrected perturbation through block size `1024`;
- negative `W+` and positive `W-` Rayleigh derivatives on normal test sizes;
- the independent large-prime aggregate targets at `m=16,24,32`;
- the documented `m=128` shape statistics;
- an ignored dense-reference regression reproducing total response through `m=1024`.

`tests/semilocal_tridiagonal_evd.rs` differentially validates dense and direct-tridiagonal eigensystems through block size `128` for both parity sectors, including `mu`, `lambda`, `mu'`, `lambda'`, ordering, and derivative signs.

`tests/semilocal_high_m_diagnostics.rs` is a no-harness scientific regression target executed by `cargo test --all-targets`. It evaluates the independently validated direct-tridiagonal route at `m=1024`, checks the documented high-block aggregate targets and all derivative signs, and prints runtime, Linux process peak RSS when available, minimum/maximum spectral values, minimum Rayleigh denominator, relative minimum `mu`, and finite-value counts.

The high-block diagnostics do not use the Stieltjes quadrature path.

## Independent reconstruction check

The exact implementation reproduces, to floating-point precision, the historical targets recorded in the issue and Phase-4 notes, including

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

The documented dense Thor reproduction also reaches `m=1024`; see `PHASE4_RUST_VALIDATION_2026-08-14.md`.

## Validation status

The original PR #11 was prepared while GitHub Actions could not obtain a runner because of the then-current account/billing state. That limitation is historical, not current repository status.

Subsequent PR #16 (`perf: add validated tridiagonal semilocal EVD path`) passed the complete repository CI on commit `452fdf33303d6b89b50ff6eacf324673171b7b6f`, including rustfmt, clippy with `-D warnings`, all Rust tests, the release smoke test, Python compilation, and the research-harness regression suite. GitHub Actions run: <https://github.com/Memorithm/riemann_ndim_bench/actions/runs/33123860169>.

For a reproducible high-block diagnostic on the current checkout, run:

```bash
cargo test --test semilocal_high_m_diagnostics
```

The normal complete scientific regression remains:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets
```

## Scientific boundary

This remains a finite-compression perturbation calculation. It does not identify these finite generalized-prolate crossings with zeta zeros and does not imply the Riemann hypothesis.
