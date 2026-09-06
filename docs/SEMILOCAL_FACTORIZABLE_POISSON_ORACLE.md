# Factorizable semilocal Poisson oracle

Status: `SOURCE THEOREM / CERTIFIED MANUFACTURED REGRESSION`.

## Purpose

The Riemann-specific semilocal bridge requires a faithful additive Fourier/Poisson control before the quotient map `E`, the Hilbert space `L^2(X_S)`, or Conjecture 4.1 can be treated numerically.

Two earlier controls delimit the problem:

- the archimedean Gaussian fixture verifies the source Poisson convention with certified tails;
- the scalar finite-place shortcut was refuted: merely deleting multiples of primes from a real sum is not semilocal Poisson theory.

The finite-place Fourier oracle then supplied the exact local identity

```text
F_p[1_{p^k Z_p}] = p^{-k} 1_{p^{-k} Z_p}
```

under the source-locked self-dual character and Haar normalization.

This increment combines those independently controlled factors into the first genuinely factorizable semilocal additive Poisson fixture in the bench.

## Manufactured class

Let

```text
S = {infinity} union S_f
```

and choose one integer `k_p` for each finite prime `p in S_f`.

The manufactured Bruhat--Schwartz function is

```text
f = f_infinity tensor prod_p 1_{p^{k_p} Z_p},
```

where `f_infinity` is the already-certified three-Gaussian archimedean fixture.

The diagonal lattice is

```text
Q_S = Z[S_f^{-1}].
```

For a diagonal rational `q in Q_S`, all finite indicators are non-zero exactly when

```text
v_p(q) >= k_p
```

for every finite place. Therefore the surviving diagonal lattice is

```text
A Z,
A = prod_p p^{k_p}.
```

This is an exact arithmetic reduction, not a fitted numerical observation.

## Fourier side

Place by place,

```text
F_p[1_{p^{k_p} Z_p}]
  = p^{-k_p} 1_{p^{-k_p} Z_p}.
```

Consequently the product finite Fourier scale is

```text
A^{-1} = prod_p p^{-k_p},
```

and the dual surviving diagonal lattice is

```text
A^{-1} Z.
```

The semilocal Poisson comparison therefore becomes

```text
sum_{n in Z} f_infinity(A n)

    ?=

A^{-1} sum_{n in Z} F_infinity f_infinity(A^{-1} n).
```

This is exactly the scaled archimedean Poisson identity expected from the product local Fourier theorem, but the implementation evaluates both finite prefixes independently rather than substituting one side for the other.

## Certification

`archimedean_poisson` now exposes certified bilateral lattice sums for both the source fixture and its closed-form Fourier transform.

For `|n| <= N`, the code evaluates the finite sum directly. Each omitted Gaussian tail is bounded by the standard monotone integral estimate, and the total two-sided tail bound is propagated through the exact finite-place Fourier scale.

The tests cover:

- the self-dual local unit ball `Z_p`;
- non-trivial exponents `k = -2,-1,1,2`;
- more than one finite place simultaneously;
- invalid/composite places and duplicate local specifications.

The asserted comparison budget is

```text
certified omitted-tail bound + explicitly declared floating-point allowance.
```

No zeta value or fitted spectral target appears in this oracle.

## Why this is not PR #43 again

PR #43 kept only a positive real scalar, replaced the integer index set by `M_S`, and discarded the finite local Fourier data. That route fails.

The present fixture instead starts from the additive diagonal lattice `Q_S` and includes explicit p-adic Bruhat factors together with their independently known local Fourier transforms. The finite local factors determine both the surviving lattice and the Fourier scaling.

Thus this regression tests an actual product Fourier datum rather than prime deletion alone.

## Deliberate boundary

This increment still does **not** implement:

- arbitrary elements of `Q_p`;
- arbitrary Bruhat--Schwartz functions on `A_S`;
- a general semilocal Fourier transform;
- the quotient `X_S = A_S / Q_S^*`;
- the quotient Poisson map `E` as a Hilbert-space operator;
- cutoff operators from Theorem 2.5;
- Conjecture 4.1;
- Weil positivity;
- RH.

The identity tested here is a manufactured additive semilocal Poisson regression on a factorized class whose local transforms are known exactly.

## Proof-status classification

- local p-adic ball Fourier law: `SOURCE THEOREM / STANDARD LOCAL FACT`;
- diagonal-lattice reduction to `A Z`: `EXACT ALGEBRA` for this manufactured class;
- archimedean Fourier pair: `EXACT CLOSED FORM`;
- finite lattice sums and tail certification: `NUMERICAL WITH RIGOROUS TRUNCATION BOUND`;
- agreement of the two semilocal sides: `CERTIFIED MANUFACTURED REGRESSION`;
- general semilocal Poisson theorem: `SOURCE THEOREM`, not re-proved here;
- quotient/Weil bridge: `OPEN IMPLEMENTATION / OPEN MATHEMATICAL BRIDGE`.

## Next Riemann-specific step

The next useful increment should move from this additive product fixture toward the quotient geometry rather than toward more generic Fourier machinery. A safe target is an explicit manufactured `Q_S^*`-orbit / `M_S` representative calculation showing how the faithful additive Poisson datum descends to the source map `E` without reverting to the scalar shortcut rejected by PR #43.
