# Exact semilocal arithmetic for `Q_S`

Status: source-locked exact arithmetic layer.

Primary source:

- Alain Connes, Caterina Consani, *The Scaling Hamiltonian*, arXiv:1910.14368, equations (2.13)–(2.14) and the discussion immediately preceding equation (4.6).

## 1. Source ring

For a finite set of places `S` containing the archimedean place, the source defines

```text
Q_S = { q in Q : the denominator of q uses only primes p in S }.
```

Equivalently, equation (2.13) is

```text
Q_S = {q in Q : |q|_v <= 1 for every v not in S}.
```

RiemannBench represents the elementary rational form directly: reduced rational numerators and denominators, with denominator factors rejected unless they belong to the finite prime set.

## 2. Source unit group

Equation (2.14) gives

```text
Q_S^* = { +/- p_1^n1 ... p_k^nk : n_j in Z }.
```

The implementation represents a unit by:

- its sign;
- one integer valuation for each finite prime in `S`.

No floating-point logarithms or approximate factorization are involved.

## 3. Unit × monoid decomposition

Immediately before equation (4.6), the source uses the decomposition

```text
q = u m,
```

with

```text
u in Q_S^*,
m in M_S,
```

where `M_S` is the positive-integer monoid prime to every finite prime in `S`.

For a non-zero reduced rational, RiemannBench obtains this decomposition by removing every finite-place prime valuation from numerator and denominator. The remaining positive numerator is exactly the `M_S` component.

The implementation then recomposes the rational exactly as an audit check.

Zero is represented in `Q_S`, but it deliberately has no `unit × positive-monoid` decomposition.

## 4. Why this matters for the Riemann route

The source derives the semilocal Poisson map

```text
E(f)(x) = |x|^(1/2) sum_{m in M_S} f(m x)
```

from the fact that the lattice sum over `Q_S` is invariant under `Q_S^*` and that each non-zero lattice element has the `u m` decomposition above.

Before this PR, RiemannBench implemented `M_S` and the finite `E` sum but did not represent the source ring whose quotient produces that monoid. This layer closes that arithmetic bookkeeping gap.

## 5. Scientific boundary

This is `EXACT` finite arithmetic for the source definitions. It does not construct:

- the local fields `Q_v`;
- the semilocal adele ring `A_S`;
- the quotient `X_S = A_S / Q_S^*`;
- the self-dual additive character;
- the semilocal Fourier transform;
- Poisson summation as a Hilbert-space theorem;
- Conjecture 4.1 or Weil positivity.

No RH proof claim is made.
