# Phase 3 specification — faithful reproduction of Q-epsilon

Primary source: Alain Connes, Caterina Consani, *Weil positivity and Trace formula, the archimedean place*, arXiv:2006.13771.

This document is an implementation specification. It does **not** claim a new mathematical result.

## Source equations to reproduce

### Epsilon — equation (84)

For rho >= 1,

```text
epsilon(rho) = sum_n lambda(n) / sqrt(1-lambda(n)^2)
               * <xi_n | theta(rho^-1) zeta_n>.
```

The source also imposes the symmetry

```text
epsilon(rho^-1) = epsilon(rho).
```

### Q-epsilon — equation (99)

The source derives

```text
Q epsilon(rho) = sum_n lambda(n)^2 / (1-lambda(n)^2) * C_n(rho),
```

with

```text
C_n(rho)
 = rho^(1/2) integral_(rho^-1)^1
     x (xi_n^an)'(x) * rho x (xi_n^an)'(rho x) dx
 + rho^(-3/2) (xi_n^an)'(rho^-1) xi_n^an(1)
 - rho^(3/2) xi_n^an(1) (xi_n^an)'(rho).
```

The paper notes that this gives

```text
Q epsilon(1) = 0.
```

### Derivative at the boundary — equation following Lemma 4

```text
epsilon'(1+) = sum_n lambda(n)^2/(1-lambda(n)^2) * xi_n(1)^2.
```

The paper reports the first numerical contributions

```text
t(0) = 11.9719
t(1) = 8.77574
t(2) = 2.20528
t(3) = 0.0433983
t(4) = 0.000125459...
```

and a total value of approximately

```text
epsilon'(1+) ~= 22.9965.
```

These numbers are regression targets only; tolerances must reflect the precision and normalization of our independently reproduced prolate implementation.

## Truncation acceptance criterion

Section 6 states that Appendix F proves the sum of the first **11 terms** in the Q-epsilon series gives a uniform approximation up to `1e-11` on the relevant interval.

We will not mark the kernel implementation complete until we have reproduced this convergence behavior or can explain, from normalization differences verified against the source, why an equivalent bound is represented differently in code.

## Toeplitz reproduction targets — equations (105) and (106)

With

```text
omega = log(q),
N = floor(a/omega),
```

the discretized quadratic form uses

```text
Q_q(xi) = omega * sum_j sum_k conj(xi(j)) xi(j+k) Qepsilon(q^|k|),
```

and the corresponding Toeplitz first row is

```text
[Qepsilon(1), Qepsilon(q), ..., Qepsilon(q^N)].
```

After normalization by `2 epsilon'(1+)`, the finite operator is

```text
T_q = omega/(2 epsilon'(1+)) * Toeplitz(Qepsilon).
```

For the paper's numerical setting

```text
omega = 1e-3,
q = exp(1e-3),
a = log(2),
```

the paper reports approximately

```text
largest eigenvalue = 1.05177
second eigenvalue  = 0.687925
```

These values are **with the source Q-epsilon kernel**. They must not be used as evidence about RH until the whole upstream kernel has passed independent regression checks.

## Required implementation order

1. Reproduce the prolate eigenvalues/functions and the exact normalization used by the paper.
2. Reproduce the analytic continuation `xi_n^an` and its derivative needed in equation (99).
3. Reproduce the individual `t(n)` boundary contributions.
4. Reproduce `epsilon'(1+)`.
5. Reproduce `Qepsilon(1)=0`.
6. Demonstrate stable convergence of the 11-term truncation.
7. Feed the verified kernel into the existing `SymmetricToeplitz` infrastructure.
8. Only then compare the first two eigenvalues to `1.05177` and `0.687925`.

## Scientific firewall

A failure to reproduce a published number is a bug, normalization mismatch, or numerical issue until demonstrated otherwise. A successful reproduction is a validation of the implementation, not a proof of the Riemann hypothesis.

Source: https://arxiv.org/html/2006.13771
