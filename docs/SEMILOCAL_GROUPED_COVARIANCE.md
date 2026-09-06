# Finite grouped Q_S covariance audit

Status: `EXACT GROUPING / CERTIFIED MANUFACTURED REGRESSION`.

## Purpose

PR #47 certified a manufactured `Q_S^*` unit-covariance mechanism for one explicit non-zero `q in Q_S`: the local exponents are transported by the exact unit action while the archimedean scale compensates them.

PR #48 then made the finite arithmetic descent

```text
q = u m,
u in Q_S^*,
m in M_S
```

executable for finite sets of explicit `Q_S` samples, without erasing their local coordinates.

This increment combines those two layers. It compares:

1. an **ungrouped direct total**, evaluated in the original input order from the original factorizable lattice fixture;
2. a **grouped transported total**, independently re-evaluated after the exact unit transport and organized by the unique `m in M_S` representative.

A caller-supplied finite weight `w(m)` is evaluated once per represented `m` and applied to both sides. This lets regressions distinguish different orbit groups without silently changing the coefficient between the direct and grouped calculations.

## API

```text
compare_grouped_unit_covariance(
    samples,
    places,
    original_balls,
    max_abs_n,
    weight_for_m,
)
```

returns a `FiniteGroupedCovarianceAudit` containing:

- the input-order direct total;
- the grouped transported total;
- the total residual;
- the sum of the rigorous Gaussian truncation bounds from both independently evaluated sides;
- one `GroupedCovarianceGroupAudit` per represented `m`.

Floating-point roundoff is intentionally not folded into the certified Gaussian tail bound. Tests declare a separate roundoff allowance.

## Exact consistency checks

For every grouped term, the implementation independently recomputes the unit-covariance data and requires that:

```text
monoid representative == grouped m
```

and

```text
transported local balls == exact orbit-grouping transport
```

A disagreement fails closed as `InconsistentOrbitData` rather than being hidden by the final numerical sum.

## Regressions

### Dyadic case

For

```text
S_f = {2},
(k_2) = (1),
```

explicit samples are split into two finite orbit groups with representatives `m=3` and `m=5`. The weights are chosen as exact binary64 values `w(3)=1`, `w(5)=2`. All prime-power scales are dyadic, so the only tolerance beyond the certified Gaussian tails is a small floating-point accumulation allowance.

### Multi-prime case

For

```text
S_f = {2,3},
(k_2,k_3) = (1,-2),
```

samples with representatives `m=5` and `m=7` are regrouped while carrying the non-trivial transport of `45/8 = (2^-3 3^2)*5`. The manufactured weight is `w(m)=1/m`; the test therefore includes an explicit roundoff allowance.

## Scientific boundary

This is still a **finite manufactured regression**. It does not establish that an infinite sum over `Q_S` may be rearranged by analytic continuation, absolute convergence, or Hilbert-space descent. It does not construct `X_S`, prove the general Poisson map

```text
E(f)(x) = |x|^(1/2) sum_{m in M_S} f(mx),
```

establish Conjecture 4.1, Weil positivity, or RH.

The result is narrower: the exact finite orbit decomposition, local unit transport, and independently evaluated manufactured covariance fixture are mutually consistent when combined in one weighted finite regrouping.
