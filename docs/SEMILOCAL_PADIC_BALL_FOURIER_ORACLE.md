# Elementary p-adic ball Fourier oracle

Status: `SOURCE THEOREM / EXACT SYMBOLIC ORACLE` for one local finite place.

## Purpose

The semilocal Poisson bridge requires more than deleting prime multiples from a real lattice sum. PR #43 records that scalar shortcut as a negative control.

The next faithful step is therefore local: represent a class of p-adic Bruhat--Schwartz functions whose Fourier transform is known exactly under the same self-dual normalization already required by `SemilocalSpaceContract`.

## Source-normalized local identity

At a finite place `p`, choose the standard additive character with conductor `Z_p` and self-dual Haar measure normalized by

```text
vol(Z_p) = 1.
```

For every integer `k`, the annihilator of the additive subgroup

```text
p^k Z_p
```

is

```text
p^(-k) Z_p,
```

and its Haar volume is

```text
vol(p^k Z_p) = p^(-k).
```

Consequently the local Fourier transform satisfies

```text
F_p[1_{p^k Z_p}]
  = p^(-k) 1_{p^(-k) Z_p}.
```

In particular

```text
F_p[1_{Z_p}] = 1_{Z_p}.
```

This is the elementary finite-place self-duality needed by standard adelic Poisson summation.

## Exact representation

`PadicBall` stores the pair

```text
(p, k)
```

representing `p^k Z_p`.

Its Fourier image is not approximated numerically. It is represented symbolically as

```text
scale = p^(-k)
ball  = p^(-k) Z_p.
```

`PadicPowerScale` therefore stores an exact prime/exponent pair rather than a floating-point approximation to `p^e`.

Applying the transform twice gives exactly

```text
p^(-k) p^(k) 1_{p^k Z_p}
= 1_{p^k Z_p},
```

which is the expected Fourier-square identity on these even local balls.

## Diagonal Q_S audit hook

For a non-zero reduced `Q_S` rational `q`, the implementation evaluates

```text
q in p^k Z_p
iff
v_p(q) >= k.
```

This provides an exact bridge between the localized arithmetic already present in RiemannBench and the local p-adic Fourier fixture.

It is only an audit hook on the diagonal copy of `Q_S`; it is not a representation of arbitrary `Q_p` values.

## Deliberate boundary

This module does **not** implement:

- arbitrary p-adic numbers;
- arbitrary locally constant compactly supported functions on `Q_p`;
- numerical p-adic integration;
- a general local Fourier transform;
- the semilocal product space `A_S`;
- the quotient `X_S`;
- the global/semilocal Poisson identity;
- Conjecture 4.1;
- Weil positivity;
- RH.

It is a manufactured exact local oracle from which a later factorizable semilocal Fourier test can be built without guessing the finite-place transform.

## Proof status

- annihilator and volume formulas for p-adic balls: `SOURCE THEOREM / STANDARD LOCAL FOURIER FACT`;
- symbolic transform law in the implementation: `EXACT` representation of that formula;
- double-transform regression on the ball class: `EXACT` symbolic identity;
- diagonal `Q_S` membership tests: `EXACT` finite arithmetic;
- full semilocal Fourier/Poisson relation: still `OPEN IMPLEMENTATION` in RiemannBench.

## Next step

The next Riemann-specific increment should combine:

1. an archimedean Fourier pair already validated by the analytic Poisson oracle;
2. one or more elementary p-adic ball Fourier pairs from this module;
3. explicit independent local coordinates rather than a single real scalar;

into a tightly manufactured factorizable semilocal Fourier fixture.

That fixture must remain clearly separated from the full Hilbert-space theorem on `L^2(X_S)` and from Conjecture 4.1.
