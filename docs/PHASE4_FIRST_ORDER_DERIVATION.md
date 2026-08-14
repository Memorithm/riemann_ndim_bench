# Exact first-order single-prime semilocal prolate perturbation

This note isolates the finite-dimensional `q=1/p` derivative used in Phase 4. Source formulas and algebraic consequences are separated from large-`m` numerical observations.

## Source input

For `S={infinity,p}`, Proposition 7.2 of Connes–Consani–Moscovici, *On q-series and the moment problem associated to local factors*, arXiv:2403.01247v1, gives

```text
a_n(q)^2
  = (n + 1/2)(n + 1)
    * [1 + 2*sqrt(2)*(alpha_{n+1}-alpha_n) q + O(q^2)],
```

where

```text
alpha_n = (-4)^(-n) binom(2n,n),
q = 1/p.
```

The exact recurrence

```text
alpha_0 = 1,
alpha_{n+1} = -((n + 1/2)/(n + 1))*alpha_n
```

is preferable numerically.

Define

```text
A_n = (n + 1/2)(n + 1),
r_n = sqrt(2)*(alpha_{n+1}-alpha_n).
```

Then

```text
a_n(0) = sqrt(A_n),
(a_n^2)'(0) = 2 A_n r_n.
```

Using the recurrence for `alpha_n`, this derivative simplifies exactly to

```text
(a_n^2)'(0)
  = -(1/sqrt(2))*(2n+1)*(4n+3)*alpha_n.
```

## Generalized prolate matrix

The first cyclic-pair prolate construction of arXiv:2310.18423 is

```text
W_lambda = -s^2 + 2*pi*lambda^2*(4N+1) - 1/4.
```

For one parity block, degree `d` is

```text
d = 2i       for W+
d = 2i + 1   for W-.
```

At a zero crossing, with `mu=lambda^2`, the generalized problem is transformed to

```text
K(q) = B^(-1/2) (-A(q)) B^(-1/2),
B_d  = 2*pi*(4d+1).
```

Its non-zero entries are

```text
K_ii(q)
  = [a_{d-1}(q)^2 + a_d(q)^2 + 1/4] / B_d,

K_i,i+1(q)
  = a_d(q)*a_{d+1}(q)
    / sqrt(B_d * B_{d+2}),
```

with `a_{-1}=0`.

## Closed form for K'(0)

### Diagonal

Starting from

```text
K'_ii(0)
  = [(a_{d-1}^2)'(0) + (a_d^2)'(0)] / B_d,
```

and using

```text
alpha_{d-1} = -(2d)/(2d-1) * alpha_d,
```

the numerator is

```text
(1/sqrt(2))*[
    2d(4d-1) - (2d+1)(4d+3)
] * alpha_d

= -(3/sqrt(2))*(4d+1)*alpha_d.
```

Since `B_d=2*pi*(4d+1)`, the factor `(4d+1)` cancels exactly:

```text
K'_ii(0) = -3*alpha_d / (2*sqrt(2)*pi).
```

This also holds for `d=0` when `a_{-1}` is omitted.

### Off diagonal

The derivative is

```text
K'_{i,i+1}(0)
  = a_d(0)*a_{d+1}(0)*(r_d+r_{d+1})
    / sqrt(B_d*B_{d+2}).
```

Because

```text
r_d+r_{d+1}
  = sqrt(2)*(alpha_{d+2}-alpha_d),
```

and

```text
alpha_{d+2}-alpha_d
  = -(4d+5) alpha_d / [4(d+1)(d+2)],
```

we obtain

```text
K'_{i,i+1}(0)
 = -sqrt(2)*(4d+5)*alpha_d/(16*pi)
   * sqrt(
       (2d+1)(2d+3)
       / [(d+1)(d+2)(4d+1)(4d+9)]
     ).
```

## Immediate structural consequences

Since `alpha_d` has sign `(-1)^d`:

- every entry of the `W+` first-order tridiagonal perturbation is negative;
- every entry of the `W-` first-order tridiagonal perturbation is positive.

Using the central-binomial asymptotic

```text
|alpha_d| ~ 1/sqrt(pi*d),
```

we also obtain

```text
K'_ii(0) = O(d^(-1/2)),
K'_{i,i+1}(0) = O(d^(-1/2)).
```

More precisely, the ratio of the leading off-diagonal magnitude to the leading diagonal magnitude tends to `1/6`.

These are algebraic consequences of the source first-order coefficient; they do not depend on the Stieltjes quadrature.

## Eigenvalue derivative

For a fixed block size, let

```text
K(0) u_j = mu_j u_j,
||u_j||=1,
lambda_j=sqrt(mu_j).
```

The finite archimedean tridiagonal matrix has simple spectrum, so standard symmetric-matrix perturbation theory gives

```text
mu'_j(0) = u_j^T K'(0) u_j,

lambda'_j(0)
  = mu'_j(0)/(2*lambda_j).
```

Thus each fixed finite crossing has

```text
lambda_j(q)
  = lambda_j(0) + lambda'_j(0) q + O(q^2).
```

This is a finite-dimensional perturbation statement only.

## Independent large-p validation

The Thor Stieltjes experiment used

```text
v_j(p) = p*(lambda_j(p)-lambda_j(0))/sqrt(m)
```

at `p=1009,4001,16001,64007`. Its limiting aggregates agree closely with the direct derivative above:

| m | direct m*mean_abs | p=64007 | direct m*trimmed_mean_abs | p=64007 |
|---:|---:|---:|---:|---:|
| 16 | 0.282632080029 | 0.282631523985 | 0.225633429048 | 0.225632461293 |
| 24 | 0.261137449994 | 0.261136720091 | 0.187448254450 | 0.187447090630 |
| 32 | 0.245755692889 | 0.245754717957 | 0.160480915987 | 0.160479565847 |

The final two large-p increments of the scaled vector decay with effective order close to one in `1/p`, consistent with the unscaled crossing remainder being `O(p^-2)` at fixed `m`.

## Large-m study — reproduced in Rust

The exact derivative was independently evaluated beyond the sizes reached by the semilocal quadrature and has now been reproduced directly in Rust on Thor through block size `m=1024`. Full output and low-edge follow-up are recorded in `PHASE4_RUST_VALIDATION_2026-08-14.md`.

Let

```text
v_j = lambda'_j(0)/sqrt(m)
```

on the merged `W+ union W-` spectrum. With one eighth removed from each edge for the trimmed mean:

| m | m*mean_abs | m*trimmed_mean_abs | m*RMS | sqrt(m)*Linf |
|---:|---:|---:|---:|---:|
| 16 | 0.2826320800 | 0.2256334290 | 0.4132069031 | 0.2524050638 |
| 32 | 0.2457556929 | 0.1604809160 | 0.4141162673 | 0.2485066705 |
| 64 | 0.2094148973 | 0.1034553692 | 0.4138497346 | 0.2474067194 |
| 128 | 0.1754882382 | 0.0705903427 | 0.4130934071 | 0.2471011380 |
| 256 | 0.1450150592 | 0.0497878013 | 0.4121956353 | 0.2470172106 |
| 512 | 0.1184232308 | 0.0352144134 | 0.4113396754 | 0.2469943647 |
| 1024 | 0.0957325576 | 0.0249083011 | 0.4106039832 | 0.2469881902 |

The last-doubling effective exponents are approximately

```text
trimmed L1 : 1.499539430
RMS        : 1.002582610
Linf       : 0.500036066
```

and `m^(3/2)*trimmed_mean_abs(v)` is `0.7970656337` at `m=1024`. These are strong numerical asymptotic signals, but they are not proved `m -> infinity` theorems.

The untrimmed response is dominated by the low spectral edge. A simple harmonic `C/k` law on a moving window was tested and is not stable at the available sizes. A double-scaling description in `x=k/sqrt(m)` gives a better finite-size collapse, but the small-`x` limit of the scaling profile remains unresolved.

## Required follow-up

1. Isolate the exact derivative implementation from the exploratory test file into focused Rust code/tests.
2. Quantify the `W+`/`W-` pairing and the low-edge double-scaling profile.
3. Determine whether the untrimmed L1 growth has a genuine logarithmic asymptotic or a different slowly varying correction.
4. Keep finite-compression statements separate from any zeta-zero interpretation.

No statement in this note identifies these finite crossing parameters with Riemann-zeta zeros or proves RH.
