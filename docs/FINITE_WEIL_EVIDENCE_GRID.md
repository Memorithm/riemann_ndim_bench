# Finite Weil evidence grid

Status: **CONSOLIDATED FINITE NUMERICAL DIAGNOSTICS / NOT WEIL POSITIVITY**

This layer combines the existing support-window, dimension and quadrature-refinement audits into one deterministic table.

## Reuse strategy

For each pair `(support window, quadrature level)`, the implementation computes the Gram-normalized Weil problem once at `N_max`. It then extracts the leading-principal subproblems `N=1..N_max` from the already-computed pairing and Gram matrices.

The results are regrouped by `(support window, dimension)` across quadrature levels. This avoids recomputing every dimension independently.

## Recorded data

For every quadrature sample the grid records:

- raw minimum eigenvalue;
- generalized minimum eigenvalue;
- Gram condition number;
- maximum critical-boundary residual;
- maximum independently evaluated pairing asymmetry;
- maximum whitening asymmetry.

For every `(support, dimension)` cell it also records:

- observed raw minimum/maximum and span across refinement levels;
- observed generalized minimum/maximum and span;
- the absolute raw and generalized delta between the last two supplied levels;
- the maximum Gram condition number seen across those levels.

No significance score, p-value, confidence level, or automatic proof-oriented classification is manufactured from these quantities.

## Release probe

```text
cargo run --release --example weil_evidence_grid -- 4
```

The default probe combines the three manufactured support windows

```text
(3/4, 13/4)
(1/2, 7/2)
(1/4, 15/4)
```

with the refinement levels

```text
(48, 48, 64, 64)
(72, 72, 96, 96)
(96, 96, 128, 128)
```

where each level is `(correlation, archimedean, boundary, Gram)` quadrature order.

The example emits CSV-compatible rows. It is intentionally a release research probe rather than a heavy default-CI workload.

## Interpretation

A spectral feature is more credible numerically when its scale can be compared directly with the observed quadrature drift, Gram conditioning and symmetry/boundary residuals. This comparison remains diagnostic: the observed refinement span is not a rigorous error enclosure, and no fixed ratio between eigenvalue magnitude and drift is declared sufficient by this module.

## Scientific boundary

Even a grid containing only positive generalized minima establishes positivity only on the tested finite spans, support windows and numerical resolutions. It does not establish:

- positivity for all admissible test functions;
- completeness or density of the chosen family;
- uniform convergence in dimension, support or quadrature order;
- equality of the numerical `L^2(d^*rho)` normalization with semilocal `L^2(X_S)`;
- the semilocal trace formula or Conjecture 4.1;
- the Riemann hypothesis.

A stable negative cell would likewise require independent replication and stronger numerical control before mathematical interpretation.
