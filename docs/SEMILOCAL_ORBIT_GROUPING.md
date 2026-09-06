# Exact finite Q_S orbit grouping by M_S

Status: `SOURCE DECOMPOSITION / EXACT FINITE ALGEBRA`.

## Purpose

Immediately before equation (4.6) in *The Scaling Hamiltonian*, every non-zero

```text
q in Q_S
```

is written uniquely as

```text
q = u m,
u in Q_S^*,
m in M_S.
```

This is the arithmetic step that permits the source semilocal Poisson map to be written with an `M_S` index. PR #43 showed that keeping only the scalar representative `m` while deleting the finite local coordinates is not a faithful semilocal model.

This increment therefore groups **finite explicit samples** from `Q_S` by their unique `M_S` representative while retaining, for every term:

- the original exact rational `q`;
- the complete unit decomposition `u`;
- the sign of `u`;
- every finite valuation `v_p(u)`;
- the local ball factors transported by that unit.

No infinite sum is rearranged by this code.

## API

`group_qs_samples_by_m(samples, places, local_balls)` returns groups sorted by increasing `m`.

Each `QsOrbitTerm` carries an exact `QsUnitMonoidDecomposition`, so

```text
term.recomposed_sample() == term.sample()
```

is an exact integer/rational identity.

Terms inside a group preserve their input order. The grouping changes only the organization of the finite data, not its values.

## Representative regression

For

```text
S_f = {2,3},
```

the samples

```text
5,
10,
5/2,
15,
45/8,
-5
```

all have the same representative

```text
m = 5
```

with different unit factors in `Q_S^*`. Samples `7` and `14` form a separate `m=7` group.

With local ball exponents

```text
(k_2,k_3) = (1,-2),
```

the term `45/8 = (2^-3 3^2) * 5` retains the transported factors

```text
(k'_2,k'_3) = (4,-4),
```

whereas the term `5` retains the original local factors. Thus two terms with the same scalar representative `m` remain distinguishable by their finite local data.

This distinction is precisely what the scalar negative control of PR #43 erased.

## Zero

The source split `q=u m` is a statement for non-zero `q`. Zero is therefore rejected explicitly instead of being silently assigned to an artificial orbit group.

Boundary terms associated with `f(0)` and `Ff(0)` remain a separate part of the Poisson/Weil bookkeeping.

## Scientific boundary

This increment establishes only exact finite arithmetic grouping. It does **not** prove that an infinite sum over `Q_S` may be rearranged without analytic conditions, does not construct `X_S`, and does not establish the Hilbert-space map

```text
E(f)(x) = |x|^(1/2) sum_{m in M_S} f(mx).
```

It does not prove Conjecture 4.1, Weil positivity, or RH.

## Next Riemann-specific step

Combine this exact grouping with the certified covariance mechanism of PR #47. The next manufactured regression should assign an explicit contribution to every finite `Q_S` sample, group by `m`, transport the local coordinates with its unit factor, and verify that grouped and ungrouped finite totals coincide before taking any infinite or Hilbert-space limit.

That would make the algebraic descent underlying equation (4.6) executable without confusing finite regrouping with the general quotient Poisson theorem.
