# Semilocal Poisson map — source arithmetic contract

Status: source-locked arithmetic layer for the Poisson route; not a proof of the semilocal Poisson identity.

Primary source:

- Alain Connes, Caterina Consani, *The Scaling Hamiltonian*, §4.1, especially equations (4.2)–(4.6), arXiv:1910.14368.

## Why this is Riemann-specific

The source identifies the Poisson summation formula as a conceptual reason that ratios of local factors appear in the semilocal Fourier transform. It also states that the normalization of the additive characters is needed for the Weil inequality and should play a key role in Conjecture 4.1.

This is therefore part of the Riemann-specific bridge, unlike the generic Jacobi/resolvent machinery transferred to TDI-10.xx.

## The monoid M_S

For a finite place set `S` containing `infinity`, every element of `Q_S` is written uniquely as

```text
q = u m,
```

with

```text
u in Q_S^*,
m in M_S,
```

where `M_S` is the multiplicative monoid of positive integers prime to every finite prime in `S`.

For example,

```text
S = {infinity,2,3}
```

gives

```text
M_S = {1,5,7,11,13,17,...}.
```

`SemilocalPoissonMonoid` implements exactly this divisibility condition.

## Source E map

The source replaces the archimedean sum over all positive integers by

```text
E(f)(x) = |x|^(1/2) sum_{m in M_S} f(m x).
```

Formally, after multiplicative Fourier transform, this produces for the trivial character

```text
zeta(1/2-is) product_{p in S_finite} (1-p^(-1/2+is)).
```

This formula is the source motivation for the local-factor ratios appearing in the semilocal Fourier operator.

## What the Rust implementation does

RiemannBench does not yet have a full numerical representation of an adele `x in A_S` or a Bruhat-Schwartz function on `A_S`.

Therefore `SemilocalPoissonMonoid::finite_e_sum` deliberately receives a callback

```text
m -> f(m x)
```

rather than pretending that a scalar `x` represents the full semilocal object.

The caller also supplies `max_m`. When compact support proves that all terms above this bound vanish, the finite sum equals the source `E` sum for that datum.

The returned `FiniteESum` records:

- the half-density-scaled value;
- the unscaled raw sum;
- the number of terms;
- the truncation bound used.

## Relation to earlier contracts

The construction uses the finite place set from `semilocal_trace_contract`.

The required `Q_S`-self-dual character normalization remains encoded by `SemilocalSpaceContract`. This PR does not duplicate or weaken that requirement.

The compact-support -> finite-place reduction is handled independently by `weil_support`.

Together the currently executable chain is:

```text
Weil boundary conditions
    -> source Q transform
    -> compact support
    -> finite place set S
    -> self-dual character normalization contract
    -> monoid M_S
    -> finite arithmetic part of E(f).
```

## What remains open

This PR does not establish:

1. a representation of `A_S` or `X_S`;
2. Poisson summation on `Q_S` inside the bench;
3. the unitary semilocal Fourier transform `F_alpha`;
4. the equality connecting additive Fourier transform, inversion and local-factor ratios;
5. the semilocal trace formula numerically;
6. Conjecture 4.1;
7. RH.

## Proof status

- Definition of `M_S` and equation (4.6): `SOURCE DEFINITION`.
- Finite divisibility arithmetic in Rust: exact implementation.
- Callback-driven finite `E` sum: exact finite arithmetic when the supplied support bound is complete.
- Fourier/local-factor discussion above: source derivation/formal discussion, not newly proved here.
- Conjecture 4.1: `SOURCE CONJECTURE`.
