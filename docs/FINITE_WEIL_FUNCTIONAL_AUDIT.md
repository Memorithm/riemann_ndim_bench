# Finite Weil functional audit

Status: **SOURCE-LOCKED FINITE NUMERICAL EXPERIMENT / NOT A POSITIVITY THEOREM**

This increment evaluates one compact convolution square with the Riemann--Weil source decomposition used by Connes--Consani.

Primary normalization reference:

- Alain Connes, Caterina Consani, *Spectral triples and zeta-cycles*, Enseign. Math. 69 (2023), §2.1, especially equations (2.7), (2.8), (2.31) and (2.32).

The compact test function `h1` is the `f = Qg` profile already validated by `semilocal_compact_weil`. In logarithmic coordinate

```text
phi(t) = h1(exp(t)).
```

For real `h1`, the multiplicative convolution square is represented by the ordinary log-coordinate autocorrelation

```text
theta(t) = integral phi(u+t) phi(u) du,
F(exp(t)) = theta(t).
```

The audit then evaluates

```text
psi(F)
  = Fhat(i/2) + Fhat(-i/2)
    - W_R(F)
    - sum_p W_p(F).
```

## Boundary / pole term

The two Mellin boundary moments of `h1` are computed by the existing compact Weil boundary audit. For real `h1`, multiplicative Fourier convolution gives

```text
Fhat(i/2) = Fhat(-i/2) = M_+(h1) M_-(h1),
```

so the reported pole term is

```text
2 M_+(h1) M_-(h1).
```

This keeps the residual boundary error visible instead of silently replacing it by exact zero.

## Non-archimedean term

Equation (2.31) is evaluated as the finite prime-power sum

```text
sum_{1 < m <= exp(L)} Lambda(m) m^(-1/2)
    [theta(log m) + theta(-log m)],
```

where `L` is the log-width of the support of `h1` and `Lambda(m)=log p` when `m=p^k`.

The upper integer bound is derived from the exact rational support ratio `upper/lower`, not from `exp(L)` in binary64. If that ratio is itself an integer prime power, the endpoint term is retained in the audit but set to exact zero because the compact autocorrelation vanishes at the support boundary.

## Archimedean term

The implementation uses the compact form of `W_R` in equation (2.32). The log-coordinate integral is evaluated by Gauss--Legendre quadrature. Near zero, the numerator is rearranged with `expm1` to reduce cancellation; this is only a numerical stabilization and does not alter the source formula.

## Numerical layers

The audit exposes three independently controlled quadrature orders:

1. inner log-autocorrelation;
2. archimedean `W_R` integral;
3. critical Mellin boundary moments.

Tests require convergence under refinement and verify the numerical symmetry `theta(t) ~= theta(-t)` at all represented prime powers.

## Manufactured fixture

For the current rational support `(1/2, 7/2)`, the exact support ratio is `7`. The represented prime powers are therefore

```text
2, 3, 4, 5, 7,
```

with `7` on the exact support boundary and hence contributing zero. The resulting scalar decomposition is stored only as a regression for this fixture.

## Scientific boundary

A positive value of `psi(F)` for this one manufactured convolution square proves only that this particular finite numerical experiment has that sign under the implemented source normalization. It does **not** prove:

- positivity for all admissible test functions on this support;
- a lower bound on the Weil quadratic form;
- Conjecture 4.1;
- the semilocal trace formula on `L^2(X_S)`;
- the Riemann Hypothesis.

The next mathematically relevant escalation is not to sample more arbitrary functions and infer a theorem. It is to move from scalar fixtures to a finite-dimensional basis of boundary-admissible test functions and construct the full Hermitian Weil quadratic-form matrix, with its smallest eigenvalue reported together with convergence and basis-size dependence.
