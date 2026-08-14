# Phase 4 Rust validation checkpoint — 2026-08-14

This note records results reproduced directly in Rust on Thor for the exact first-order single-prime semilocal prolate perturbation. It supersedes the earlier `pending Rust reproduction` status of the large-block targets in `PHASE4_FIRST_ORDER_DERIVATION.md`.

## Exact finite-block validation

The Rust implementation uses the source-derived recurrence

```text
alpha_0 = 1,
alpha_{n+1} = -((n+1/2)/(n+1))*alpha_n,
```

with

```text
(a_n^2)'(0) = -(1/sqrt(2)) (2n+1)(4n+3) alpha_n.
```

Three independent checks passed:

1. the unsimplified source-derived `K'(0)` and the closed form agree to about `1e-14`;
2. the sign-corrected `K'(0)` is strictly diagonally dominant through block size `m=1024`;
3. all computed `W+` crossing derivatives are negative and all `W-` crossing derivatives are positive.

The minimum diagonal-dominance margins at `m=1024` were

```text
W+  2.808620942441e-3
W-  2.807934238423e-3
```

## Agreement with the independent large-p Stieltjes experiment

For the merged normalized response

```text
v_j = lambda'_j(0)/sqrt(m),
```

Rust reproduced the previously independent large-prime limits:

| m | m*mean_abs | m*trimmed_mean_abs |
|---:|---:|---:|
| 16 | 0.2826320800294 | 0.2256334290480 |
| 24 | 0.2611374499942 | 0.1874482544496 |
| 32 | 0.2457556928890 | 0.1604809159874 |

The differences from the stored targets were about `1e-13`.

## Large-block Rust reproduction

The previously auxiliary large-block table was reproduced in Rust through `m=1024`:

| m | m*mean_abs | m*trimmed_mean_abs | m*RMS | sqrt(m)*Linf | m^(3/2)*trimmed_mean_abs |
|---:|---:|---:|---:|---:|---:|
| 16 | 0.2826320800294 | 0.2256334290480 | 0.4132069031489 | 0.2524050638280 | 0.9025337161920 |
| 32 | 0.2457556928890 | 0.1604809159874 | 0.4141162673372 | 0.2485066705391 | 0.9078171515660 |
| 64 | 0.2094148972661 | 0.1034553691615 | 0.4138497345639 | 0.2474067193632 | 0.8276429532919 |
| 128 | 0.1754882381797 | 0.07059034265272 | 0.4130934070885 | 0.2471011380024 | 0.7986385596163 |
| 256 | 0.1450150591944 | 0.04978780134720 | 0.4121956353477 | 0.2470172105957 | 0.7966048215551 |
| 512 | 0.1184232308169 | 0.03521441341300 | 0.4113396754123 | 0.2469943646897 | 0.7968112166348 |
| 1024 | 0.09573255762299 | 0.02490830105319 | 0.4106039832160 | 0.2469881901528 | 0.7970656337020 |

Effective exponents on the last doubling `512 -> 1024` were

```text
trimmed L1 : 1.499539430
RMS        : 1.002582610
Linf       : 0.500036066
```

This is strong numerical evidence, within this finite-compression model, for

```text
trimmed mean_abs(v) ~ const * m^(-3/2),
RMS(v)              ~ const * m^(-1),
Linf(v)             ~ const * m^(-1/2).
```

These remain numerical asymptotics, not proved `m -> infinity` theorems.

## Untrimmed low-edge behavior

The untrimmed total response

```text
sum_j |lambda'_j(0)|
```

was measured as

```text
m=128   3.970845543531
m=256   4.640481894221
m=512   5.359223651882
m=1024  6.126883687871
```

A global fit to `A log(m)+B` over these four sizes gives

```text
A    = 1.036844178588
B    = -1.084469067704
RMSE = 2.450593025561e-2
```

but the local logarithmic slopes are still drifting upward:

```text
128 -> 256   0.9660810423403
256 -> 512   1.036925169456
512 -> 1024  1.107499327011
```

The response is strongly concentrated at the low spectral edge. At `m=1024`, the first 128 merged modes contain about `69.7%` of the total absolute derivative, while the highest 128 contain only about `0.37%`.

A naive pair-smoothed harmonic law `|lambda'_k| ~ C/k` on a moving window was tested and rejected as a stable description at these sizes: the fitted exponents drifted from about `1.23` to `1.37`, and the inferred `2C` continued to overestimate the total logarithmic slope.

A double-scaling test with `x=k/sqrt(m)` gives a substantially better finite-size collapse for

```text
pair_abs(k,m) ~ m^(-1/2) F(k/sqrt(m)),
```

with relative RMS across `m=256,512,1024` of roughly `4%` to `8%` for sampled fixed `x`. This suggests a genuine low-edge scaling regime, but the current data are not sufficient to establish the small-`x` asymptotic of `F`.

## Scientific boundary

All statements here concern the exact first-order perturbation of finite generalized prolate compressions derived from the semilocal Jacobi `q=1/p` coefficient. They do not identify compression crossings with Riemann-zeta zeros and do not imply the Riemann hypothesis.
