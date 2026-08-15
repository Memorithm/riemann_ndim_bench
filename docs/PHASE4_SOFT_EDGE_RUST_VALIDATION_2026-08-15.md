# Phase 4 soft-edge Rust validation — 2026-08-15

This note records direct coefficient-level Rust checks for the soft-edge asymptotics used in the formal `(log m)^2` derivation.

The tests use only the exact finite coefficients already implemented in `tests/semilocal_multi_prime_probe.rs`. No eigendecomposition, no global fit, and no regression is involved.

## Leading row asymptotics

For parity degree

```text
d = 2i + epsilon,
epsilon = 0 for W+,
epsilon = 1 for W-,
```

the test `validates_soft_edge_row_asymptotic_coefficients` checks

```text
4*pi*b_i/i -> 1,
64*pi*i*V_i -> 1,
2*pi^(3/2)*sqrt(i)*W_i -> 1,
pi*i*W_i/sqrt(b_i) -> 1.
```

At `i=16384`:

```text
W+
4*pi*b/i                  1.000038147147273
64*pi*i*V                 0.9999921886371026
2*pi^(3/2)*sqrt(i)*W      0.9999961847909877
pi*i*W/sqrt(b)            0.9999771118358036

W-
4*pi*b/i                  1.000068664725392
64*pi*i*V                 0.9999622280748391
2*pi^(3/2)*sqrt(i)*W      0.9999809265214316
pi*i*W/sqrt(b)            0.9999465965815049
```

All four leading scaled quantities are within `7e-5` of one for both parities.

## First finite-size correction

The follow-up test `validates_soft_edge_first_finite_size_correction` uses the symmetric edge coefficient

```text
b_bar = (b_forward + b_backward)/2
```

and checks the first row correction predicted from the exact formulas:

```text
4*pi*b_bar/i
 = 1 + (4 epsilon + 1)/(8i) + O(i^-2),

64*pi*i*V
 = 1 - (4 epsilon + 1)/(8i) + O(i^-2),

2*pi^(3/2)*sqrt(i)*W
 = 1 - (4 epsilon + 1)/(16i) + O(i^-2),

pi*i*W/sqrt(b_bar)
 = 1 - (4 epsilon + 1)/(8i) + O(i^-2).
```

`V` is a small residual obtained by subtracting three quantities of order `i`. Its first correction becomes f64-cancellation limited before the other quantities; therefore the retained validation samples `V` at `i=4096`, while `b_bar`, `W`, and the trace prefactor are sampled at `i=16384`.

Retained errors against the predicted first-correction constants:

```text
W+
  b_bar            +2.8610011213e-6
  V                -1.8467529208e-4
  W                -8.3844570327e-6
  trace prefactor  -9.2188602139e-6

W-
  b_bar            +2.8609138099e-6
  V                -2.6241214255e-4
  W                +1.2713462638e-7
  trace prefactor  +1.3597491488e-5
```

In particular,

```text
W+ : i(C_i-1) -> -1/8,
W- : i(C_i-1) -> -5/8,
```

where

```text
C_i = pi*i*W_i/sqrt(b_bar(i)).
```

This is the row-level quantity that controls the logarithmic first finite-size correction in the formal singular trace calculation.

## Induced no-fit correction for the centered coefficient

The row analysis predicts

```text
D = 3/(8*pi^2)
```

for the coefficient of `log(m)/m` in the formal expansion of `S(m)`.

For

```text
A_m
 = [S(2m)-2S(m)+S(m/2)]/[2(log 2)^2],
```

this gives

```text
c_A
 = 3/[32*pi^2*(log 2)^2]
 = 0.01977063457049387...
```

for the coefficient in

```text
A_m
 = 1/(2*pi^2)
   + c_A log(m)/m
   + O(1/m).
```

The homogeneous Rust checkpoints give:

```text
m=256   0.02043519186892146
m=512   0.02034611805432550
m=1024  0.02028111986421427
m=2048  0.02023106505876094
m=4096  0.02053140616003269
m=8192  0.02017962874379530
```

for

```text
m*(A_m-1/(2*pi^2))/log(m).
```

At `m=8192` the relative gap to the no-fit prediction is

```text
2.068695224946514e-2
```

or about `2.07%`.

## Interpretation

The leading coefficient `1/(2*pi^2)` and the first finite-size coefficient above are now both supported by direct exact-row Rust asymptotics rather than by free regression against the global `S(m)` sequence.

The passage from row asymptotics to the global singular trace remains formal. A rigorous proof still requires uniform soft-edge control of `K^{-1/2}` or an equivalent resolvent representation.

## Scientific boundary

These are statements about finite semilocal prolate perturbations. They do not identify finite compression crossings with Riemann-zeta zeros and do not imply the Riemann hypothesis.
