# First finite-size correction to the soft-edge log-squared law

This note refines the formal soft-edge derivation of

```text
S(m) ~ [1/(2*pi^2)] (log m)^2 + O(log m)
```

by isolating the first row-level correction that can generate a `log(m)/m` term in the dyadic centered coefficient.

The result is still formal at the level of the global singular trace. The row asymptotics themselves are direct consequences of the exact finite coefficients and have been checked independently in Rust on Thor.

## 1. Symmetric local edge coefficient

Write the parity degree as

```text
d = 2 i + epsilon,
epsilon = 0 for W+,
epsilon = 1 for W-.
```

After alternating conjugation, use the symmetric frozen edge coefficient

```text
b_bar(i)
 = [K_{i,i+1} + K_{i,i-1}]/2.
```

Expanding the exact archimedean Jacobi coefficients gives

```text
b_bar(i)
 = i/(4*pi)
   * [1 + (4 epsilon + 1)/(8 i) + O(i^-2)].
```

The soft-edge residual satisfies

```text
V_i
 = 1/(64*pi*i)
   * [1 - (4 epsilon + 1)/(8 i) + O(i^-2)].
```

The Rust coefficient test confirms these constants. `V_i` is numerically delicate at very large `i` because it is obtained by cancellation of three quantities of order `i`; the `i=4096` checkpoint remains well resolved while the other row quantities remain stable through `i=16384`.

## 2. Sign-corrected perturbation

Let

```text
W_i
 = H_{ii} - H_{i,i+1} - H_{i,i-1}
```

for the sign-corrected first derivative (`H_+=-K'_+`, `H_-=K'_-`).

Using

```text
|alpha_d|
 = 4^(-d) binom(2d,d)
 ~ 1/sqrt(pi*d) * [1 - 1/(8d) + O(d^-2)]
```

and the exact forward/backward off-diagonal ratios, one finds that their `1/d` corrections cancel in the sum. Therefore

```text
W_i
 = 1/[2*pi^(3/2)*sqrt(i)]
   * [1 - (4 epsilon + 1)/(16 i) + O(i^-2)].
```

The corresponding frozen trace prefactor is

```text
C_i
 = pi*i*W_i/sqrt(b_bar(i))
 = 1 - (4 epsilon + 1)/(8 i) + O(i^-2).
```

Thus

```text
W+ : C_i = 1 - 1/(8i) + O(i^-2),
W- : C_i = 1 - 5/(8i) + O(i^-2).
```

The Rust test at `i=16384` gives

```text
W+  i(C_i-1) = -0.1250092188602139
W-  i(C_i-1) = -0.6249864025085117
```

in close agreement with `-1/8` and `-5/8`.

## 3. Why the infrared correction does not change the log/m coefficient

The local soft-edge symbol can be written formally as

```text
k(i,theta)
 = b_bar(i) theta^2 + V_i + ...
```

near `theta=0`. Its infrared scale is therefore

```text
a_i = sqrt(V_i/b_bar(i))
    = 1/(4i) * [1 + O(1/i)].
```

The singular integral is of the form

```text
integral dtheta / sqrt(theta^2 + a_i^2)
 = 2 log i + constant + O(1/i).
```

Changing `a_i` by a relative `O(1/i)` factor changes the logarithm only by `O(1/i)`. After multiplication by the leading `1/i` row prefactor this produces `O(1/i^2)`, whose tail is `O(1/m)` and contains no `log(m)/m` term.

Likewise, the `theta^2` part of the perturbation symbol removes the `1/|theta|` singularity and contributes no new logarithmic row term at this order.

Hence the coefficient of `(log i)/i^2` is controlled by the `1/i` correction to `C_i` alone.

## 4. Induced correction in S(m)

For one parity block the leading row density is

```text
(log i)/(2*pi^2*i).
```

Multiplying by

```text
1 + c_epsilon/i + O(i^-2)
```

with

```text
c_+ = -1/8,
c_- = -5/8
```

gives a logarithmic correction

```text
(c_epsilon/(2*pi^2)) * (log i)/i^2.
```

Using the tail asymptotic

```text
sum_{i<=m} (log i)/i^2
 = constant - (log m + 1)/m + O(log m/m^2),
```

the two parity corrections combine to

```text
S(m)
 = A log(m)^2 + B log(m) + C
   + D log(m)/m + O(1/m),
```

with the formal coefficient

```text
D
 = -(c_+ + c_-)/(2*pi^2)
 = 3/(8*pi^2)
 = 0.03799544386587666...
```

The sign is positive because `c_+ + c_- = -3/4` and the finite partial sum differs from its infinite limit by a negative tail.

## 5. Induced correction in the dyadic centered coefficient

Define

```text
A_m
 = [S(2m) - 2 S(m) + S(m/2)]
   / [2 (log 2)^2].
```

For

```text
f(m) = log(m)/m,
```

the centered dyadic combination is

```text
f(2m) - 2 f(m) + f(m/2)
 = [0.5 log(m) - 1.5 log(2)]/m.
```

Therefore the coefficient multiplying `log(m)/m` in `A_m` is

```text
c_A
 = D/[4 (log 2)^2]
 = 3/[32*pi^2*(log 2)^2]
 = 0.01977063457049387...
```

This is not fitted to the global `S(m)` data.

## 6. Comparison with homogeneous Rust checkpoints

For the retained centered coefficients,

| m | m(A_m-A)/log(m) | ratio to predicted c_A |
|---:|---:|---:|
| 256 | 0.02043519186892146 | 1.033613351967 |
| 512 | 0.02034611805432550 | 1.029107992552 |
| 1024 | 0.02028111986421427 | 1.025820379811 |
| 2048 | 0.02023106505876094 | 1.023288604451 |
| 4096 | 0.02053140616003269 | 1.038479877154 |
| 8192 | 0.02017962874379530 | 1.020686952249 |

where

```text
A = 1/(2*pi^2).
```

At `m=8192` the observed coefficient is about `2.07%` above the no-fit prediction.

The `m=4096` point is visibly less smooth than its neighbors; it mixes lower-dimensional direct checkpoints with the pairwise-refined `m=8192` value in the centered second difference and should not be overinterpreted in isolation.

After subtracting the predicted logarithmic correction, the residual

```text
m(A_m-A) - c_A log(m)
```

is of order a few `10^-3` over the retained range, consistent with an additional `d/m` term, but no coefficient for that term is claimed here.

## 7. Rust validation outcome

The dedicated Rust test `validates_soft_edge_first_finite_size_correction` passes after explicitly sampling `V_i` at `i=4096` to avoid the f64 cancellation floor while retaining `b_bar`, `W`, and `C_i` at `i=16384`.

Representative final errors against the predicted first corrections are:

```text
W+:
  b_bar            +2.86e-6
  V                -1.85e-4
  W                -8.38e-6
  trace prefactor  -9.22e-6

W-:
  b_bar            +2.86e-6
  V                -2.62e-4
  W                +1.27e-7
  trace prefactor  +1.36e-5
```

The test also checks that the `m=8192` observed finite-size coefficient is within `3%` of `c_A`; the measured relative gap is

```text
0.02068695224946514
```

that is, about `2.07%`.

## Scientific boundary

The row expansions and finite Rust checks are concrete. The passage from those row expansions to the global singular trace asymptotic remains formal until a uniform resolvent/soft-edge argument is supplied. In particular, this note is not a theorem about the Riemann zeta function and does not imply the Riemann hypothesis.
