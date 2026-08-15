# Phase 4 numerical checkpoints — 2026-08-15

This note records the current Rust-validated numerical state of the exact first-order single-prime semilocal prolate perturbation. It is intentionally separated from claims about the Riemann zeta zeros.

## Scientific boundary

All quantities below concern finite generalized prolate compressions derived from the first `q=1/p` coefficient of the semilocal Jacobi recurrence. No finite crossing is identified here with a Riemann-zeta zero, and nothing in this note proves the Riemann hypothesis.

## Exact first-order input

For the source coefficient

```text
alpha_0 = 1,
alpha_{n+1} = -((n+1/2)/(n+1))*alpha_n,
```

one has

```text
(a_n^2)'(0)
  = -(1/sqrt(2)) (2n+1)(4n+3) alpha_n.
```

For parity degree `d`, the generalized-prolate derivative is tridiagonal with

```text
K'_ii(0) = -3 alpha_d / (2 sqrt(2) pi)
```

and

```text
K'_{i,i+1}(0)
 = -sqrt(2)*(4d+5)*alpha_d/(16*pi)
   * sqrt(
       (2d+1)(2d+3)
       / [(d+1)(d+2)(4d+1)(4d+9)]
     ).
```

The finite-block sign lemma gives

```text
K'_+(0) < 0  as a quadratic form,
K'_-(0) > 0  as a quadratic form,
```

so every finite first-order `W+` crossing derivative is negative and every finite first-order `W-` derivative is positive.

## Tridiagonal eigensolver validation

The faer tridiagonal self-adjoint EVD implementation was checked against the earlier dense solver for both parity blocks at block sizes

```text
m = 8, 16, 32, 64, 128.
```

For the tested outputs the crossing values and first-order derivatives agreed bit-for-bit (`worst error = 0`). This enabled high-dimensional runs without storing a dense matrix.

## Direct first-order high-block reproduction

The exact derivative, using eigenvectors, was reproduced in Rust through block size `m=4096`.

Representative normalized statistics are

| m | m*mean_abs | m^(3/2)*trimmed_mean_abs | m*RMS | sqrt(m)*Linf |
|---:|---:|---:|---:|---:|
| 1024 | 9.57325576230e-2 | 7.97065633702e-1 | 4.10603983216e-1 | 2.46988190153e-1 |
| 2048 | 7.67139672726e-2 | 7.97227163828e-1 | 4.10008680259e-1 | 2.46986531198e-1 |
| 4096 | 6.10045325893e-2 | 7.97316385172e-1 | 4.09545298363e-1 | 2.46986087705e-1 |

These are numerical asymptotic signals only. In particular the trimmed response is consistent with an `m^(-3/2)` mean scale, while RMS and Linf exhibit different edge-sensitive scalings.

## Total first-order variation

Define

```text
S(m) = sum_j |lambda'_j(0)|
```

on the merged `W+ union W-` finite spectrum.

Because of the finite-block sign lemma,

```text
S(m)
 = sum_j lambda'_{-,j}(0)
   - sum_j lambda'_{+,j}(0).
```

The validated Rust checkpoints are

| m | S(m) | method |
|---:|---:|---|
| 128 | 3.970845543531 | exact first-order EVD |
| 256 | 4.640481894221 | exact first-order EVD |
| 512 | 5.359223651882 | exact first-order EVD |
| 1024 | 6.126883687871 | exact first-order EVD |
| 2048 | 6.943355708182 | exact first-order tridiagonal EVD |
| 4096 | 7.808580171428 | exact first-order tridiagonal EVD |
| 8192 | 8.72252476599664 | pairwise eigenvalue finite difference, validated window |
| 16384 | 9.68517075774107 | pairwise eigenvalue finite difference, validated window |

The earlier global-trace finite-difference value at `m=8192`,

```text
8.722522625208107
```

is retained only as a historical checkpoint. It differs from the homogeneous pairwise value by about `2.14e-6`, negligible relative to `S(8192)` but not negligible after taking second differences.

## Why pairwise finite differences were needed

The eigenvalues-only trace method first formed

```text
Tr sqrt(K(+h)) - Tr sqrt(K(-h))
```

before dividing by `2h`. At high dimension this subtracts two large nearly equal totals. At `m=16384`, reducing `h` exposed a cancellation floor: Richardson extrapolants ceased to follow the expected even-power hierarchy.

The safer method keeps the same eigenvalues-only tridiagonal solver but forms

```text
(lambda_j(+h) - lambda_j(-h)) / (2h)
```

mode by mode before summing. This substantially reduced cancellation.

### `m=16384` stable window

Pairwise centered differences gave

```text
h=5.00e-4   D=9.685679471063766
h=2.50e-4   D=9.685297949399001
h=1.25e-4   D=9.685202556488507
```

with raw difference ratios

```text
(D(5e-4)-D(2.5e-4)) / (D(2.5e-4)-D(1.25e-4))
  = 3.999476090960...
```

in excellent agreement with the expected factor `4` for a centered `O(h^2)` error.

The order-4 Richardson values were

```text
R12 = 9.685170775510747
R23 = 9.685170758851676
```

and the corresponding `h^2+h^4` extrapolate was

```text
Q123 = 9.685170757741071.
```

An independent interior point at

```text
h = 1.75e-4
```

gave

```text
observed  = 9.685232898533920
predicted = 9.685233082763617
relative prediction residual = -1.90217105068e-8.
```

The residual is about `0.296%` of the finite-`h` correction at that point. This is the current validation for the stored `S(16384)` checkpoint.

A smaller step `h=6.25e-5` was already affected by the roundoff floor and is explicitly rejected for extrapolation.

### `m=8192` stable window

The homogeneous pairwise values were

```text
h=5.00e-4   D=8.722779121464640
h=2.50e-4   D=8.722587591872799
h=1.25e-4   D=8.722540424778753
```

with

```text
Q123 = 8.722524765996642.
```

The independent interior validation at `h=1.75e-4` gave

```text
observed  = 8.722554927042479
predicted = 8.722555487118820
relative prediction residual = -6.42101248482e-8.
```

The residual is about `1.86%` of the finite-`h` correction. The extrapolated value from the three-point quadratic in `h^2` agrees with `Q123` to floating-point roundoff.

## Quadratic-log finite-size coefficient

For centered dyadic sizes define

```text
A_m = [S(2m) - 2 S(m) + S(m/2)] / [2 (log 2)^2].
```

If

```text
S(m) = A (log m)^2 + B log m + C + lower-order terms,
```

then `A_m -> A` when the lower-order corrections vanish.

Using the homogeneous checkpoints above gives

| center m | A_m | A_m - 1/(2 pi^2) | m*(gap)/log(m) |
|---:|---:|---:|---:|
| 256 | 5.11032354345484e-2 | +4.42643613380e-4 | 2.04351918689e-2 |
| 512 | 5.09084935580483e-2 | +2.47901736879e-4 | 2.03461180543e-2 |
| 1024 | 5.07978750345720e-2 | +1.37283213403e-4 | 2.02811198642e-2 |
| 2048 | 5.07359112365777e-2 | +7.53194154088e-5 | 2.02310650587e-2 |
| 4096 | 5.07022850427297e-2 | +4.16932215608e-5 | 2.05314061600e-2 |
| 8192 | 5.06827887066639e-2 | +2.21968854950e-5 | 2.01796287421e-2 |

The candidate

```text
1/(2 pi^2) = 0.05066059182116889...
```

is therefore numerically close to the observed sequence, but equality is **not established**.

The nearly constant quantity

```text
m * (A_m - 1/(2 pi^2)) / log(m)
```

is evidence that a correction of the form

```text
(c log m + d)/m
```

should be tested before interpreting naive `1/m` Richardson extrapolations. This remains a finite-size numerical observation; an analytic derivation is still required.

## Rejected numerical inferences

The following values must not be used as final scientific checkpoints:

- the global-trace order-8 extrapolation at `m=16384`;
- pairwise extrapolations involving `h=6.25e-5` at `m=16384`;
- any Richardson ratio mixing global-trace and pairwise summation methods;
- the earlier centered `A_4096` and `A_8192` values built from the obsolete global `S(8192)` checkpoint.

They are useful diagnostics of cancellation, not asymptotic data.

## Next mathematical target

For a fixed finite positive matrix `K(q)`, the trace derivative satisfies the exact finite-dimensional identity

```text
d/dq Tr sqrt(K(q)) |_{q=0}
  = (1/2) Tr(K(0)^(-1/2) K'(0)).
```

Together with the sign lemma this gives

```text
S(m)
 = (1/2) [
     Tr(K_-^(-1/2) K'_-)
     - Tr(K_+^(-1/2) K'_+)
   ].
```

The next substantive task is to exploit the explicit tridiagonal coefficients of `K(0)` and `K'(0)` to derive the large-`m` trace asymptotic analytically, rather than continue fitting larger finite matrices.
