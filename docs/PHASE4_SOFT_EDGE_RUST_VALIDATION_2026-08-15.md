# Phase 4 soft-edge Rust validation — 2026-08-15

This note records a direct Rust validation of the local asymptotic coefficients used in `PHASE4_LOG2_ASYMPTOTIC_HEURISTIC.md`.

The test uses the exact finite-dimensional coefficient generator already present in the local Phase 4 exploratory code on Thor. No eigendecomposition and no fit to the total-variation data is involved.

## Quantities tested

For row index `i`, after alternating conjugation of the positive off-diagonal generalized-prolate matrix, define

```text
b_i = K_{i,i+1},
V_i = K_{ii} - K_{i,i+1} - K_{i,i-1}.
```

For the sign-corrected first-order perturbation `H` (`H=-K'_+` for `W+`, `H=K'_-` for `W-`), define

```text
W_i = H_{ii} - H_{i,i+1} - H_{i,i-1}.
```

The formal soft-edge analysis predicts

```text
b_i ~ i/(4*pi),
V_i ~ 1/(64*pi*i),
W_i ~ 1/(2*pi^(3/2)*sqrt(i)),
pi*i*W_i/sqrt(b_i) -> 1.
```

The last limit is the combined coefficient entering the singular local trace density.

## Rust output

### W+

| i | `4*pi*b/i` | `64*pi*i*V` | `2*pi^(3/2)*sqrt(i)*W` | `pi*i*W/sqrt(b)` |
|---:|---:|---:|---:|---:|
| 64 | 1.009776959092385 | 0.997868601735470 | 0.998989743912227 | 0.994141723333574 |
| 256 | 1.002442119766469 | 0.999500529483240 | 0.999753746033940 | 0.998535218231022 |
| 1024 | 1.000610396238726 | 0.999877229608598 | 0.999938832637215 | 0.999633792825281 |
| 4096 | 1.000152590684167 | 0.999969437335134 | 0.999984732946113 | 0.999908447499060 |
| 16384 | 1.000038147147273 | 0.999992188637103 | 0.999996184790988 | 0.999977111835804 |

Final absolute errors from 1 at `i=16384`:

```text
4*pi*b/i                     3.8147147e-5
64*pi*i*V                    7.8113629e-6
2*pi^(3/2)*sqrt(i)*W         3.8152090e-6
pi*i*W/sqrt(b)               2.2888164e-5
```

### W-

| i | `4*pi*b/i` | `64*pi*i*V` | `2*pi^(3/2)*sqrt(i)*W` | `pi*i*W/sqrt(b)` |
|---:|---:|---:|---:|---:|
| 64 | 1.017589372064102 | 0.990150933875054 | 0.995118215453882 | 0.986480249900595 |
| 256 | 1.004395243378969 | 0.997553178699692 | 0.998779335513381 | 0.996591605451556 |
| 1024 | 1.001098677466935 | 0.999389308393292 | 0.999694826220265 | 0.999146107238087 |
| 4096 | 1.000274660996326 | 0.999847348043911 | 0.999923706173615 | 0.999786414433789 |
| 16384 | 1.000068664725392 | 0.999962228074839 | 0.999980926521432 | 0.999946596581505 |

Final absolute errors from 1 at `i=16384`:

```text
4*pi*b/i                     6.8664725e-5
64*pi*i*V                    3.7771925e-5
2*pi^(3/2)*sqrt(i)*W         1.9073479e-5
pi*i*W/sqrt(b)               5.3403418e-5
```

All test assertions with tolerance `1e-4` passed for both parity blocks.

## Interpretation

This validates, in Rust and independently of the global spectral fits, the four local limits used to obtain the formal soft-edge density

```text
h(i,theta)/sqrt(k(i,theta))
  ~ 1/(pi*i)
     / sqrt(theta^2 + 1/(16*i^2)).
```

Consequently the numerical candidate `1/(2*pi^2)` for the coefficient of `(log m)^2` is supported by the exact row-coefficient asymptotics, not merely by regression of `S(m)`.

This remains a validation of the local asymptotic input. A proof of the global trace asymptotic still requires uniform control of the finite-section resolvent / `K^(-1/2)` near the soft edge.

## Scientific boundary

These statements concern the asymptotics of the finite semilocal generalized-prolate perturbation only. They do not identify compression crossings with Riemann-zeta zeros and do not imply the Riemann hypothesis.
