# Exact Q_S^* unit-orbit transport

Status: `SOURCE DECOMPOSITION / EXACT FINITE ALGEBRA`.

## Purpose

The source semilocal Poisson construction uses the quotient by

```text
Q_S^* = { +/- prod_p p^{n_p} }
```

and the unique decomposition of every non-zero localized rational

```text
q = u m,
u in Q_S^*,
m in M_S.
```

RiemannBench already implements this decomposition exactly. PR #43 also established a critical negative control: replacing the ordinary integer index set by `M_S` while discarding all finite local coordinates does **not** reproduce the semilocal Poisson/Fourier theory.

The missing bookkeeping is therefore the action of the unit `u` on the finite local factors themselves.

## Local action

Write

```text
u = sign * prod_p p^{n_p}.
```

For the elementary p-adic ball factor used by the current manufactured fixtures,

```text
1_{p^k Z_p}(u_p x_p)
```

is non-zero exactly when

```text
v_p(x_p) + n_p >= k,
```

or equivalently

```text
x_p in p^{k-n_p} Z_p.
```

Thus the exact unit transport is

```text
k_p -> k_p - n_p.
```

No numerical p-adic approximation is involved.

## Archimedean compensation

The same diagonal rational unit has archimedean absolute value

```text
|u_infinity| = prod_p p^{n_p}.
```

If the original local product selects the diagonal lattice scale

```text
A = prod_p p^{k_p},
```

then the transported finite local factors select

```text
A_u = prod_p p^{k_p-n_p}.
```

Multiplication of the archimedean coordinate by `u` restores the original effective scale place by place:

```text
n_p + (k_p - n_p) = k_p,
```

and globally

```text
|u_infinity| A_u = A.
```

The implementation keeps the exponent identity exact. A floating-point archimedean scale is exposed only as an audit hook.

## Concrete regression

For

```text
S_f = {2,3},
q = 45/8,
```

the exact source decomposition is

```text
q = u m,
u = 2^-3 3^2 = 9/8,
m = 5.
```

Starting from local exponents

```text
k_2 = 1,
k_3 = -2,
```

the transported product has

```text
k'_2 = 1 - (-3) = 4,
k'_3 = -2 - 2 = -4.
```

The compensated exponents are exactly

```text
-3 + 4 = 1,
 2 + (-4) = -2.
```

The external regression also covers a negative unit sign and three finite places.

## Relation to PR #45

PR #45 certified additive semilocal Poisson summation on a narrow factorized class by retaining explicit finite-place Fourier data.

This increment supplies the exact `Q_S^*` action needed to ask how that additive datum descends toward the source quotient map

```text
E(f)(x) = |x|^(1/2) sum_{m in M_S} f(mx).
```

It does **not** assert that the quotient descent is now proved. It only establishes the finite algebra that a faithful descent must preserve.

## Proof-status classification

- unique `q = u m` decomposition: `SOURCE FACT / EXACT IMPLEMENTATION` from the existing `semilocal_qs` layer;
- transport `k_p -> k_p-n_p`: `EXACT LOCAL ALGEBRA`;
- compensation `n_p + (k_p-n_p)=k_p`: `EXACT`;
- numerical evaluation of `|u_infinity|`: `AUDIT HOOK`, not needed for the exact identity;
- quotient map `E` as an operator on `X_S`: `OPEN IMPLEMENTATION`;
- semilocal Hilbert-space Fourier/Poisson theorem: `SOURCE THEOREM`, not independently re-proved;
- Conjecture 4.1: `SOURCE CONJECTURE`;
- Weil positivity / RH: `OPEN`.

## Next Riemann-specific step

Use the exact unit transport together with the certified factorizable additive Poisson fixture to build a manufactured quotient-orbit regression: group explicit diagonal `Q_S` samples by their unique `Q_S^*` orbit representative `m in M_S`, carry the local coordinates through the unit action, and verify that the grouped calculation agrees with the ungrouped additive fixture within certified truncation bounds.

That is the safe next step toward making the source map `E` faithful without reintroducing the scalar shortcut already refuted by PR #43.
