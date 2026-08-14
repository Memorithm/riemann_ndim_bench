# Finite-block sign lemma for the first-order semilocal perturbation

This note proves a finite-dimensional consequence of the first `q=1/p` coefficient in Proposition 7.2 of arXiv:2403.01247. It is a derived lemma for the bench, not a theorem stated in the source paper.

## Setup

For parity degree `d=2i` (`W+`) or `d=2i+1` (`W-`), the exact first-order generalized-prolate matrix derivative has

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
     ),
```

where

```text
alpha_d = (-4)^(-d) binom(2d,d).
```

Thus `alpha_d>0` for even `d` and `alpha_d<0` for odd `d`.

Define the sign-corrected tridiagonal matrix

```text
H+ = -K'_+,
H- =  K'_-.
```

All diagonal and off-diagonal entries of each `H` are positive.

Let

```text
D_d = 3 |alpha_d| / (2 sqrt(2) pi)
```

be its diagonal magnitude and let `O_d` be the off-diagonal magnitude joining degree `d` to degree `d+2`.

## Forward off-diagonal bound

The exact ratio is

```text
R_d = O_d / D_d
```

with

```text
R_d^2
 = (2d+1)(2d+3)(4d+5)^2
   / [144(d+1)(d+2)(4d+1)(4d+9)].
```

For every integer `d>=0`,

```text
1/25 - R_d^2
 = [704 d^4 + 5472 d^3 + 11484 d^2 + 7408 d + 717]
   / [3600(d+1)(d+2)(4d+1)(4d+9)]
 > 0.
```

Hence

```text
R_d < 1/5.
```

## Backward off-diagonal bound

For `d>=2`, let

```text
L_d = O_{d-2} / D_d.
```

Using the exact recurrence for `alpha_d`,

```text
L_d^2
 = d(d-1)(4d-3)^2
   / [9(2d-3)(2d-1)(4d-7)(4d+1)].
```

Moreover

```text
1/4 - L_d^2
 = [512 d^4 - 1856 d^3 + 1776 d^2 - 108 d - 189]
   / [36(2d-3)(2d-1)(4d-7)(4d+1)].
```

Writing `d=x+2`, the numerator becomes

```text
512 x^4 + 2240 x^3 + 2928 x^2 + 1108 x + 43,
```

which is strictly positive for `x>=0`. Therefore

```text
L_d < 1/2.
```

## Strict diagonal dominance

For every interior row,

```text
(O_{d-2} + O_d) / D_d
  = L_d + R_d
  < 1/2 + 1/5
  = 7/10 < 1.
```

Boundary rows have only one off-diagonal and satisfy an even stronger bound.

Therefore both `H+` and `H-` are real symmetric, have positive diagonal, and are strictly diagonally dominant by absolute row sums. By Gershgorin, every eigenvalue of each `H` is strictly positive. Hence

```text
K'_+(0) is negative definite,
K'_-(0) is positive definite.
```

## Consequence for every finite crossing branch

For a normalized eigenvector `u_j` of the archimedean generalized matrix `K(0)`,

```text
mu'_j(0) = u_j^T K'(0) u_j.
```

Since `lambda_j=sqrt(mu_j)>0`,

```text
lambda'_j(0) = mu'_j(0)/(2 lambda_j).
```

Consequently, for **every finite block size and every simple crossing branch**,

```text
W+ : lambda'_j(0) < 0,
W- : lambda'_j(0) > 0.
```

So adjoining a sufficiently large prime (`q=1/p` small and positive) shifts every `W+` crossing downward to first order and every `W-` crossing upward to first order.

This statement is only about the finite compressed generalized-prolate crossings. It does not identify them with zeta zeros and has no direct implication for RH.

## Rust regression target

The eventual Rust implementation should test, for a range of block sizes:

1. the closed forms for `K'(0)` against the unsimplified source-derived formulas;
2. strict diagonal dominance of the sign-corrected matrix;
3. negative Rayleigh derivatives for all `W+` branches;
4. positive Rayleigh derivatives for all `W-` branches.
