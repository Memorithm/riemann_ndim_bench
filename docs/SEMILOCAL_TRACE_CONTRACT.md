# Semilocal trace formula — source data contract

Status: Riemann-specific source contract; **not** a numerical implementation of the semilocal Hilbert space.

Primary source:

- Alain Connes, Caterina Consani, *The Scaling Hamiltonian*, §2.2–2.3, especially equations (2.12)–(2.30), arXiv:1910.14368.

## Why this contract exists

The current bench already contains finite Jacobi/prolate models. Those matrices must not be silently identified with the source cutoff operator `R_Lambda` acting on `L^2(X_S)`.

Before implementing any numerical semilocal trace object, the source-level data and normalizations are fixed explicitly.

## Finite set of places

The source assumes a finite set `S` containing the archimedean place.

It defines

```text
A_S = product_{v in S} Q_v,
Q_S = {q in Q : |q|_v <= 1 for all v not in S},
Q_S^* = {+/- p_1^n1 ... p_k^nk : p_j in S\{infinity}, n_j in Z},
X_S = A_S / Q_S^*,
C_{Q,S} = A_S^* / Q_S^*.
```

`FinitePlaceSet` stores only the finite primes and makes the inclusion of `infinity` implicit and mandatory.

This is metadata for the source construction. The Rust type does **not** pretend to realize the quotient spaces numerically.

## Basic additive characters and Haar measure

For each `v in S`, the source chooses a basic additive character `alpha_v` and normalizes additive Haar measure to be self-dual. Their product

```text
alpha = product_v alpha_v
```

defines the Fourier transform `F_alpha`, which extends unitarily to `L^2(X_S)`.

Section 4.1 stresses that, for the Poisson formula relevant to the Weil inequality, one uses choices for which `Q_S` is a self-dual lattice. Different such choices differ by scaling with `Q_S^*`, so the ambiguity disappears after passage to `X_S`.

RiemannBench therefore currently exposes exactly one accepted policy:

```text
BasicCharacterNormalization::QsSelfDual
```

No p-adic phase formula is invented here.

## Scaling and symmetric test normalization

The source scaling action satisfies

```text
theta_a(lambda)^* theta_a(lambda) = |lambda|_S.
```

Under the unitary identification with the idele-class representation, the test functions in equations (2.29) and (2.30) are related by

```text
h(lambda) = |lambda|^(1/2) f(lambda).
```

`symmetric_test_value` implements only this scalar normalization.

## Infrared and ultraviolet cutoffs

The source defines

```text
P_Lambda: support cutoff |x| <= Lambda,
P_hat_Lambda = F_alpha P_Lambda F_alpha^-1,
R_Lambda = P_hat_Lambda P_Lambda.
```

The present code stores only the shared positive scalar `Lambda`; it does not fabricate `P_Lambda`, `P_hat_Lambda`, `F_alpha`, or `R_Lambda`.

Lemma 2.4 contains the frequency endpoint

```text
2 log Lambda / (2 pi) = log Lambda / pi.
```

`SemilocalCutoff::quantized_band_endpoint` records this normalization.

## Theorem 2.5 normalization

For compactly supported `f`, the source theorem states asymptotically

```text
Tr(theta_a(f) R_Lambda)
  = 2 f(1) log Lambda
    + sum_{v in S} local_v(f)
    + o(1).
```

The code exposes only the exact leading normalization

```text
2 f(1) log Lambda.
```

It does **not** compute the trace, the principal-value local distributions, or the `o(1)` remainder.

## What is executable now

`src/semilocal_trace_contract.rs` verifies mechanically that later code uses:

- a finite set of actual prime places;
- the mandatory archimedean place;
- the `Q_S`-self-dual character convention;
- the source cutoff parameter;
- the exact band normalization from Lemma 2.4;
- the exact leading divergence from Theorem 2.5;
- the symmetric `h = |lambda|^(1/2) f` convention.

## What remains open

Still absent from RiemannBench:

1. an actual representation of Bruhat-Schwartz data on `A_S`;
2. the quotient Hilbert space `L^2(X_S)`;
3. the self-dual Fourier transform `F_alpha`;
4. the scaling representation `theta_a`;
5. the local factor-ratio unitary `u = product_v u_v o pi_v`;
6. the cutoff projections and `R_Lambda`;
7. a numerical or formal verification of Theorem 2.5;
8. the support-window positivity bridge of Conjecture 4.1.

Any future implementation must preserve this distinction.

## Scientific status

- Definitions and normalizations above: `SOURCE THEOREM` / `SOURCE DEFINITION`.
- Rust representation: implementation contract.
- Theorem 2.5 itself: source theorem, **not re-proved by these tests**.
- Semilocal support-window sufficiency: `SOURCE CONJECTURE`.
- RH: not proved.
