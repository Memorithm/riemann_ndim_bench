# Exact compact archimedean bound for the finite semilocal E bridge

Status: `EXACT SUPPORT BOOKKEEPING / MANUFACTURED ARCHIMEDEAN FIXTURE`.

Primary source:

- Alain Connes, Caterina Consani, *The Scaling Hamiltonian*, §4.1, especially equation (4.6), arXiv:1910.14368.

## Purpose

The finite bridge introduced after the exact `Q_S -> Q_S^* x M_S` descent still accepted a caller-supplied integer truncation `max_m` for

```text
E(f)(x) = |x|^(1/2) sum_{m in M_S} f(m x).
```

That parameter is harmless for a manufactured regression but is not an auditable consequence of the test function itself.

This increment removes that choice for one compactly supported archimedean fixture. The support endpoints and the positive archimedean coordinate are stored as exact reduced rationals, so the active integer `m` range is determined before any floating-point evaluation.

## Exact support arithmetic

Let

```text
support = (a,b),
x_infinity = r > 0,
```

with `a`, `b`, and `r` represented by `PositiveRational`.

An integer `m >= 1` is active exactly when

```text
a < m r < b.
```

The implementation computes the open integer interval from the exact rational quotients

```text
a/r,
b/r
```

using `u128` cross-products. Consequently:

- an integer lying exactly on either boundary is excluded;
- a non-integral upper quotient uses its exact floor;
- no `f64` rounding decision can add or remove an `m` at the support boundary.

For example,

```text
a = 1/2,
b = 7/2,
r = 1/2
```

gives exactly

```text
1 < m < 7,
```

hence

```text
m = 2,3,4,5,6
```

and `max_m = 6`.

## Numerical profile inside the exact support

After exact support membership is known, `CompactArchimedeanBump` evaluates the standard smooth bump

```text
exp(-1 / (t (1-t))),
t = (rho-a)/(b-a),
rho = m r.
```

The numerical value is evaluated in `f64`. This does not affect the truncation decision: outside the exact rational support the value is returned as exact zero without evaluating the bump.

The `MultiplicativeSupport` already used by the Weil-boundary layer is retained as the numerical support view; the rational endpoints are authoritative for discrete `m` inclusion.

## Finite E integration

`compare_compact_bruhat_e_bridge`:

1. derives the active integer `m` bounds exactly;
2. decomposes each explicit non-zero `q in Q_S` as `q = u m`;
3. removes samples whose representative lies outside the exact open archimedean support, because their manufactured archimedean contribution is exactly zero;
4. passes only active samples to the finite Bruhat descent bridge;
5. supplies the mechanically derived `max_m` to `SemilocalPoissonMonoid::finite_e_sum`;
6. uses the compact bump itself for the represented archimedean values.

Thus the finite enumeration bound is no longer a free caller parameter for this fixture.

## Scientific boundary

This increment does **not** establish that the manufactured bump is the full archimedean component required by the general semilocal Poisson theorem. In particular it does not:

1. represent a general Bruhat–Schwartz function on `A_S`;
2. construct `X_S = A_S / Q_S^*`;
3. prove Poisson summation on `Q_S`;
4. prove the Hilbert-space properties of `E`;
5. impose or prove the Weil boundary conditions on this bump;
6. prove Conjecture 4.1, Weil positivity, or RH.

The result is narrower: for this explicit finite compact-support regression, truncation is now derived from exact support data instead of being supplied externally.

## Next Riemann-specific step

The next useful increment is to connect a compact archimedean fixture to the existing source-locked Weil boundary operator

```text
Q = -(rho d/drho)^2 + 1/4
```

and numerically audit its two critical Mellin boundary moments. That should remain a separate test layer: satisfying those moments numerically is a validation of the manufactured test function, not a proof of Weil positivity.
