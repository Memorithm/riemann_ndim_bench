# Manufactured Q_S^* unit covariance

Status: `EXACT BOOKKEEPING / CERTIFIED MANUFACTURED REGRESSION`.

## Purpose

The source discussion in *The Scaling Hamiltonian* states that additive summation over `Q_S` is invariant under multiplication by the unit group `Q_S^*`, and that every non-zero `q in Q_S` has a unique decomposition

```text
q = u m,
u in Q_S^*,
m in M_S.
```

PR #46 made the finite local action of `u` exact. This increment checks that the already-certified factorized Poisson fixture is numerically covariant under that action while retaining the finite local coordinates.

## Mechanism

For local ball exponents `k_p` and unit valuations `n_p`, PR #46 gives

```text
k_p -> k_p - n_p.
```

The transported finite factors therefore select the lattice scale

```text
A_u = prod_p p^(k_p-n_p).
```

The real coordinate is simultaneously rescaled by

```text
|u_infinity| = prod_p p^n_p.
```

Hence the exact exponent identity implies

```text
|u_infinity| A_u = prod_p p^k_p = A.
```

The implementation evaluates the original manufactured archimedean lattice sum at `A` and independently re-evaluates it at the transported effective step `|u_infinity| A_u`.

## Regressions

Two distinct cases are used.

### Dyadic exact-audit case

```text
S_f = {2},
q = 3/8 = 2^-3 * 3,
k_2 = 1.
```

The transport gives `k'_2=4`, so

```text
|u_infinity| = 1/8,
A_u = 16,
|u_infinity| A_u = 2 = A.
```

All these scale values are exactly representable in binary64. No roundoff allowance is required beyond the independent Gaussian tail bounds.

### Multi-prime case

```text
S_f = {2,3},
q = 45/8 = 2^-3 3^2 * 5,
(k_2,k_3) = (1,-2).
```

The transported exponents are `(4,-4)`. The exact exponent compensation is already proved by PR #46; this test only audits the corresponding floating-point lattice evaluation with an explicitly declared roundoff allowance.

## Scientific boundary

This is not yet the quotient Poisson map `E` on `X_S`. It verifies covariance of a narrow manufactured additive fixture under the exact unit action. It does not establish the Hilbert-space descent, Conjecture 4.1, Weil positivity, or RH.

## Next step

The next source-specific increment should group finite explicit `Q_S` samples by the unique representative `m in M_S`, carry the associated unit transport on the local coordinates, and compare the grouped result with the ungrouped additive calculation. That will test the concrete algebraic descent underlying equation (4.6) without reverting to the scalar prime-deletion shortcut rejected by PR #43.
