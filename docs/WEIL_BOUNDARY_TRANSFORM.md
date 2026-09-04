# Weil boundary transform — executable source contract

Status: source-locked Riemann-specific implementation note.

Primary source:

- Alain Connes, Caterina Consani, *Weil positivity and Trace formula, the archimedean place*, Lemma 3.3 and equations (55)–(56), arXiv:2006.13771.

## Source statement represented here

For compactly supported smooth functions on the multiplicative group `R_+^*`, the two boundary conditions are

```text
integral f(rho) rho^(+1/2) d*rho = 0,
integral f(rho) rho^(-1/2) d*rho = 0,
```

where `d*rho = d rho / rho`.

The source identifies the ideal defined by these conditions with the range of

```text
Q = -(rho d/drho)^2 + 1/4.
```

It also states that applying `Q` does not enlarge compact support.

In logarithmic coordinate

```text
x = log rho,
```

one has

```text
rho d/drho = d/dx,
Q = -d^2/dx^2 + 1/4.
```

Therefore the two critical multiplicative characters satisfy

```text
Q rho^(+1/2) = 0,
Q rho^(-1/2) = 0.
```

This algebraic annihilation is the mechanism behind the two moment conditions after integration by parts for a compactly supported smooth preimage.

## Repository API

`src/weil_boundary.rs` exposes:

- `MultiplicativeSupport` for finite positive support intervals;
- `character_multiplier(exponent) = 1/4 - exponent^2`;
- `q_from_log_second_derivative`;
- `q_from_rho_derivatives`;
- `q_on_support`, which returns exact zero outside the declared interval;
- `mellin_power_moment` with the source Haar measure `d*rho`;
- `critical_boundary_moments` for powers `+1/2` and `-1/2`.

## Independent manufactured regression

`tests/weil_boundary_transform.rs` uses the standard smooth logarithmic bump

```text
g(x) = exp(-1/(1-u^2)),  |u| < 1,
u = x / log(2),
```

with zero extension outside its support.

The test computes its second logarithmic derivative analytically, applies `Q`, and checks:

1. exact zero outside the declared support;
2. numerical vanishing of both critical Mellin moments;
3. the factor-two support convention used in the archimedean source window.

The critical-character annihilation itself is tested exactly through

```text
1/4 - (+1/2)^2 = 0,
1/4 - (-1/2)^2 = 0.
```

The quadrature residual is only a numerical regression of the integration-by-parts consequence; it is not presented as a proof of Lemma 3.3.

## What is not implemented here

This PR does not yet implement:

- the unique compactly supported inverse `g` for an arbitrary `f` in the ideal;
- the source convolution formula involving the distributions `Y` and `Y*`;
- the Sonin projection;
- the semilocal Hilbert space `L^2(X_S)`;
- the semilocal trace formula;
- any positivity theorem.

Those are separate proof dependencies.

## Proof status

- Differential identity for `Q`: `SOURCE THEOREM` / exact implementation.
- Support preservation contract: `SOURCE THEOREM`, represented directly by the local wrapper and manufactured test.
- Two critical-character zeros: exact algebra.
- Numerical moment residual of the manufactured bump: `BENCH REGRESSION`.
- Weil positivity: not claimed.
- RH: not claimed.
