# Exact q=0 semilocal m=1024 diagnostic checkpoint — 2026-08-28

This note records a reproducible high-block validation of the exact first-order single-prime semilocal prolate derivative already implemented for issue #7.

The calculation is a finite-compression perturbation calculation. It does not identify finite generalized-prolate crossings with zeta zeros and does not establish any implication for the Riemann hypothesis.

## Scope and solver separation

The source-derived dense reference implementation remains `src/semilocal.rs` and uses

```text
faer::linalg::solvers::SelfAdjointEigen
```

with its `U()` and `S()` factors to evaluate

```text
mu'_j(0) = u_j^T K'(0) u_j,
lambda'_j(0) = mu'_j(0)/(2 lambda_j).
```

The exact dense route and the source-derived formulas are documented in:

- `PHASE4_FIRST_ORDER_DERIVATION.md`;
- `PHASE4_FIRST_ORDER_SIGN_LEMMA.md`;
- `PHASE4_RUST_VALIDATION_2026-08-14.md`.

The current high-block diagnostic uses the later direct-tridiagonal faer path in `src/semilocal_tridiagonal.rs`. PR #16 validated that path differentially against the dense `SelfAdjointEigen` oracle through block size 128 for `mu`, `lambda`, `mu'`, `lambda'`, ordering, and derivative signs. The direct-tridiagonal route uses the same exact `K(0)` and `K'(0)` and does not use the Stieltjes quadrature path.

The historical dense Thor reproduction through `m=1024` remains recorded in `PHASE4_RUST_VALIDATION_2026-08-14.md`; the checkpoint below adds reproducible runtime, memory, denominator-stability, and sign diagnostics on the current GitHub Actions environment.

## Reproduction command

The no-harness test target is part of the normal all-target suite:

```bash
cargo test --test semilocal_high_m_diagnostics
```

and therefore also runs under:

```bash
cargo test --all-targets
```

## Exact CI environment

Successful GitHub Actions run:

- PR: `#17`;
- branch head commit: `c78ffefc5a178a60390c7ce299eb954e34863f44`;
- pull-request merge test commit: `0648fb3f6bff2b08df839a106c2242f4cccee28e`;
- workflow run: `33199569541`;
- runner OS: Ubuntu `24.04.4` LTS;
- runner image: `ubuntu-24.04`, image version `20260823.283.1`;
- Rust: `rustc 1.98.0 (88d9e12ae 2026-08-18)`;
- faer: `0.24.4`.

The complete CI passed:

```text
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets
cargo run --release
python research-harness compilation
python research-harness regression tests
```

## m=1024 runtime and memory

The diagnostic target reported:

```text
W+ elapsed              = 7.096762 s
W- elapsed              = 7.091576 s
total diagnostic time   = 14.188543 s
RSS before              = 4272 KiB
peak RSS after W+       = 33856 KiB
process peak RSS        = 35308 KiB
```

The RSS figure is Linux `VmHWM` for the diagnostic process. It is an observed process peak on this runner, not a portable analytical memory bound.

## Denominator and spectrum stability

For `W+`:

```text
min mu                   = 1.917241265476625e-4
max mu                   = 3.212770180442295e2
min mu / max mu          = 5.967564306802309e-7
min lambda               = 1.384644815639240e-2
max lambda               = 1.792420201973381e1
min 2*lambda denominator = 2.769289631278480e-2
```

For `W-`:

```text
min mu                   = 1.916305735494191e-4
max mu                   = 3.214354017767818e2
min mu / max mu          = 5.961713379738283e-7
min lambda               = 1.384306951327700e-2
max lambda               = 1.792861962831444e1
min 2*lambda denominator = 2.768613902655400e-2
```

No zero or non-finite Rayleigh denominator was encountered. The minimum observed denominators are comfortably separated from zero at this finite block size; no NaN or infinity was produced.

## Derivative sign diagnostics

For `W+`:

```text
min lambda'       = -2.469881901563398e-1
max lambda'       = -1.569217214465107e-4
negative count    = 1024
zero count        = 0
positive count    = 0
non-finite values = 0
```

For `W-`:

```text
min lambda'       =  1.568446486865127e-4
max lambda'       =  2.469871624562106e-1
negative count    = 0
zero count        = 0
positive count    = 1024
non-finite values = 0
```

Thus the finite-block sign lemma is reproduced at `m=1024`: every computed `W+` first-order crossing derivative is negative and every computed `W-` derivative is positive.

## Aggregate regression

The merged normalized vector is

```text
v_j = lambda'_j(0) / sqrt(m).
```

At `m=1024`, the CI diagnostic returned

```text
sum_j |lambda'_j|        = 6.126883687873787
m * mean_abs(v)          = 0.09573255762302792
m * trimmed_mean_abs(v)  = 0.02490830105318828
m * RMS(v)               = 0.4106039832167578
sqrt(m) * Linf(v)        = 0.2469881901563398
```

These reproduce the previously documented high-block targets within floating-point error. The diagnostic assertions use an absolute tolerance of `5e-8` for the independent tridiagonal high-block route; the actual observed discrepancies against the stored decimal checkpoints are many orders of magnitude smaller.

The small-block `m=16,24,32` regression targets remain covered by the dense exact implementation with the tighter existing `5e-12` tolerance.

## Independence from quadrature

`tests/semilocal_high_m_diagnostics.rs` calls only the exact semilocal tridiagonal eigensystem path. It does not import or call the Stieltjes quadrature support used by the independent large-prime validation route.

The two numerical routes therefore remain separable: quadrature can serve as an external regression, but it is not an input to the exact `q=0` derivative calculation.

## Interpretation boundary

All values in this note describe finite compressed generalized-prolate matrices and their first derivative with respect to `q=1/p` at `q=0`. No large-`m` theorem is inferred from this checkpoint alone.

This calculation does not identify the crossings with zeta zeros and does not establish an RH implication.
