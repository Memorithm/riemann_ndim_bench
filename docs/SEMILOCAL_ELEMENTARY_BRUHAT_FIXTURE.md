# Elementary semilocal Bruhat--Schwartz fixture

Status: source-locked finite-place fixture for the Riemann-specific semilocal Poisson bridge.

## Purpose

The repository already represents:

- the localized arithmetic ring `Q_S`;
- the unit/monoid split `Q_S^* x M_S`;
- the finite semilocal Poisson monoid `M_S`;
- the source contract requiring a `Q_S`-self-dual additive-character normalization.

Until this change, however, the finite local content of a test function was not represented explicitly. `SemilocalPoissonMonoid::finite_e_sum` therefore had to receive scalar values from a callback and could not say whether those values arose from an archimedean factor, a finite local factor, or a genuine semilocal test function.

This fixture adds one deliberately narrow class:

```text
phi_f = product_{p in S_f} 1_{Z_p}
```

restricted to the diagonal copy of `Q_S`.

## Exact finite-place statement

For a reduced non-zero rational `q=a/b` and `p in S_f`,

```text
v_p(q) = v_p(a) - v_p(b).
```

The elementary local factor is

```text
1_{Z_p}(q) = 1  iff  v_p(q) >= 0.
```

The product finite factor is therefore

```text
phi_f(q) = product_{p in S_f} 1_{Z_p}(q).
```

Zero belongs to every `Z_p` and receives value `1`.

These are exact finite arithmetic statements in the implementation.

## Compatibility with the Poisson monoid

Every `m in M_S` is prime to all finite places in `S`. Consequently

```text
v_p(mq) = v_p(q)
```

for every declared finite place `p`, and hence

```text
phi_f(mq) = phi_f(q).
```

The integration tests verify this directly against `SemilocalPoissonMonoid`.

This is the finite-place mechanism behind the source observation that the `Q_S` lattice sum can be reduced modulo the unit group to a sum indexed by `M_S`.

## Fourier normalization

For the standard self-dual local additive character and self-dual Haar normalization, the characteristic function of `Z_p` is the elementary local self-dual Bruhat--Schwartz factor used in the standard Poisson framework.

RiemannBench does **not** turn that source fact into a general semilocal Fourier transform in this PR. The module merely fixes the finite factor that a later faithful Fourier/Poisson implementation may use.

## Deliberate boundary

This fixture does **not** represent:

- a general element of the semilocal adele ring `A_S`;
- independent `Q_p` coordinates;
- the quotient `X_S = A_S / Q_S^*`;
- a general Bruhat--Schwartz function on `A_S`;
- the semilocal Fourier transform `F_alpha`;
- the Hilbert space `L^2(X_S)`;
- the cutoff projections of Theorem 2.5;
- Conjecture 4.1;
- Weil positivity or RH.

The method `evaluate_factorizable_diagonal` therefore accepts an externally supplied archimedean value and multiplies it by the exact finite factor only on a diagonal `Q_S` sample. It must not be interpreted as evaluation on arbitrary `A_S`.

## Proof status

- p-adic valuation arithmetic on `Q_S`: `EXACT` in the bench.
- invariance under multiplication by `M_S`: `EXACT` in the bench.
- use of the standard local unit-ball factor in the self-dual Poisson normalization: `SOURCE THEOREM / STANDARD LOCAL FOURIER FACT`.
- full semilocal Poisson/Fourier intertwining: `OPEN IMPLEMENTATION / SOURCE-FORMAL` in the current bench.
- Conjecture 4.1: `SOURCE CONJECTURE`.

## Next Riemann-specific step

The next faithful extension should introduce a narrowly defined semilocal Fourier/Poisson oracle on a manufactured factorizable class where every local Fourier factor is known independently. It must keep the distinction between a finite oracle and the full Hilbert-space statement on `X_S` explicit.
