# Phase 4 — Semilocal Jacobi/prolate findings

Status: research notes. These notes distinguish source-backed statements, finite-dimensional derivations, numerical observations, and open questions. Nothing here is a proof of the Riemann hypothesis, and finite prolate/Jacobi eigenvalues or zero-crossing parameters must not be identified with zeta zeros.

## Primary sources

1. Alain Connes, Caterina Consani, Henri Moscovici, *Zeta zeros and prolate wave operators*, arXiv:2310.18423v2 (2024).
2. Alain Connes, Caterina Consani, Henri Moscovici, *On q-series and the moment problem associated to local factors*, arXiv:2403.01247v1 (2024).

## 1. Source-backed framework

For a finite set of places `S` containing the archimedean place, arXiv:2403.01247 associates a determinate moment problem whose measure is governed by the squared modulus of the product of the local factors restricted to the critical line. For `S={infinity,p}`, the moments, orthogonal polynomials, and Jacobi matrix are power series in

```text
q = 1/p.
```

For the even measure used here the Jacobi diagonal is zero and multiplication by the spectral variable is represented by off-diagonal coefficients `a_n`.

The prolate construction used by the bench is the first cyclic-pair construction in arXiv:2310.18423. In the orthogonal-polynomial representation it has the form

```text
W_lambda = -s^2 + 2*pi*lambda^2*(4*N + 1) - 1/4,
```

where `s` is multiplication by the spectral variable and `N` is the polynomial grading. The paper states that this formulation is meaningful in the general orthogonal-polynomial setting and extends to the semilocal setting.

## 2. Finite parity blocks and generalized zero crossings

For a parity block of size `m`, with polynomial degree

```text
d = 2*n       for W+
d = 2*n + 1   for W-
```

write

```text
W(lambda) = A + lambda^2 B.
```

The current bench uses

```text
A_nn     = -(a_{d-1}^2 + a_d^2) - 1/4
A_n,n+1  = -a_d*a_{d+1}
B_nn     = 2*pi*(4*d + 1).
```

`B` is positive diagonal. A zero of one finite compressed eigenvalue is therefore the generalized symmetric eigenproblem

```text
(-A) v = mu B v,       mu = lambda^2,
```

or equivalently the ordinary symmetric problem

```text
K = B^(-1/2) (-A) B^(-1/2).
```

The crossing parameters are `lambda_j = sqrt(mu_j)`.

### Numerical controls already passed

- Generalized-eigenvalue crossings reproduce direct bisection of the compressed prolate eigenvalues at approximately machine precision (`~1e-15`).
- Leading principal compressions satisfy the expected Cauchy interlacing.
- The archimedean case reproduces the explicit published parity-block coefficients.

These controls validate the finite matrix algebra. They do **not** turn the crossings into zeta zeros.

## 3. Archimedean scaling control

The exact archimedean Jacobi coefficients are

```text
a_n(0) = sqrt((n + 1/2)*(n + 1)).
```

Using these exact coefficients, the generalized crossings were pushed to block size `m=128`. The normalized bulk quantiles stabilize under division by `sqrt(m)`, while the smallest crossing scales approximately as `1/sqrt(m)`.

For example, for `W+`:

| m | min*sqrt(m) | q50/sqrt(m) | max/sqrt(m) |
|---:|---:|---:|---:|
| 16 | 0.440992833766 | 0.208763200161 | 0.500363471537 |
| 32 | 0.442156684250 | 0.218135385827 | 0.523812390475 |
| 64 | 0.442659044960 | 0.222980141973 | 0.538664862586 |
| 128 | 0.442891884827 | 0.225437588276 | 0.548061951426 |

The upper edge converges much more slowly; no limiting constant is claimed for it.

## 4. Multi-prime finite-S experiments

Direct Stieltjes quadrature and independent moment/Gram-Schmidt reconstruction agree for the tested low orders. The cases explored include

```text
{infinity}
{infinity,2}
{infinity,2,3}
{infinity,2,3,5}.
```

At block sizes through `m=32`, the normalized semilocal bulk approaches the archimedean benchmark in merged `W+ union W-` Wasserstein-type distances. On the finite grid `m=12,16,20,24,28,32`, a `1/m`-like decay was visible for several trims, but the fitted exponent changes under deep trimming. This must therefore be treated as a finite-size/pre-asymptotic observation, not as an asymptotic theorem.

## 5. Prime identity at fixed cardinality

At fixed number of finite primes the spectral correction depends strongly on the identity of the prime(s), not only on cardinality.

For singleton sets, with the untrimmed mean coefficient `C = mean_m(m D_m)` over `m=12,16,20,24,28,32`, high-resolution `170/190` quadrature gives approximately

| finite prime set | C |
|---|---:|
| {2} | 0.108897469278 |
| {3} | 0.069878275078 |
| {5} | 0.035980958657 |
| {7} | 0.022398815288 |

The same ordering survives the tested trims.

For pairs, untrimmed values are approximately

| finite prime set | C |
|---|---:|
| {2,3} | 0.136458550358 |
| {2,5} | 0.117142870775 |
| {2,7} | 0.114660392834 |
| {3,5} | 0.082430758712 |
| {3,7} | 0.079392823720 |
| {5,7} | 0.052194099817 |

The pair response is not a simple sum of singleton responses.

## 6. Exact first-order q coefficient from the source

Proposition 7.2 of arXiv:2403.01247 gives, for both parities,

```text
a_n(q)^2
  = (n + 1/2)(n + 1)
    * [1 + 2*sqrt(2)*(alpha_{n+1} - alpha_n)*q + O(q^2)],
```

with

```text
alpha_n = (-4)^(-n) * binom(2*n, n).
```

Hence, defining

```text
r_n = sqrt(2)*(alpha_{n+1} - alpha_n),
```

we have

```text
a_n(q) = a_n(0) * [1 + r_n*q + O(q^2)].
```

For numerical work `alpha_n` should be generated by the stable recurrence

```text
alpha_0 = 1,
alpha_{n+1} = -((n + 1/2)/(n + 1))*alpha_n,
```

rather than by factorial/binomial evaluation.

## 7. Finite-dimensional first-order prolate perturbation

For fixed block size `m`, define

```text
K(q) = B^(-1/2) (-A(q)) B^(-1/2).
```

At `q=0`, the derivative follows directly from the source coefficient above. Let

```text
A2_n = (n + 1/2)(n + 1).
```

Then

```text
(d/dq a_n(q)^2)|_{q=0} = 2*A2_n*r_n.
```

For degree `d` in the chosen parity block,

```text
K'_nn(0)
  = [ (a_{d-1}^2)'(0) + (a_d^2)'(0) ] / B_nn,

K'_{n,n+1}(0)
  = a_d(0)*a_{d+1}(0)*(r_d + r_{d+1})
    / sqrt(B_nn*B_{n+1,n+1}).
```

Because the finite symmetric tridiagonal archimedean matrix has simple eigenvalues, if

```text
K(0) u_j = mu_j u_j,
||u_j|| = 1,
lambda_j(0) = sqrt(mu_j),
```

standard finite-dimensional perturbation theory gives

```text
mu'_j(0) = u_j^T K'(0) u_j,

lambda'_j(0)
  = [u_j^T K'(0) u_j] / [2*lambda_j(0)].
```

Thus, at **fixed finite m**, the source q-series plus ordinary symmetric-matrix perturbation imply

```text
lambda_j(q)
  = lambda_j(0) + lambda'_j(0) q + O(q^2)
```

for each simple branch near `q=0`.

This is a finite-dimensional statement. It says nothing by itself about the limit `m -> infinity` and nothing about RH.

## 8. Independent numerical validation of the first-order derivative

The Thor run compared the merged spectral vector

```text
v_j(p) = p*(lambda_j(p)-lambda_j(infinity))/sqrt(m)
```

at `p=1009,4001,16001,64007` using independent Stieltjes quadrature. The vector converges strongly; the successive changes decrease by about a factor four when `p` is multiplied by about four, consistent with an `O(1/p)` remainder in `v_j(p)` and therefore an `O(1/p^2)` remainder in the unscaled crossing.

A separate evaluation of the exact derivative above predicts the following aggregate first-order limits:

| m | exact m*mean(|v|) | Thor p=64007 | exact m*trimmed-mean(|v|) | Thor p=64007 |
|---:|---:|---:|---:|---:|
| 16 | 0.282632080029 | 0.282631523985 | 0.225633429048 | 0.225632461293 |
| 24 | 0.261137449994 | 0.261136720091 | 0.187448254450 | 0.187447090630 |
| 32 | 0.245755692889 | 0.245754717957 | 0.160480915987 | 0.160479565847 |

The agreement is strong and independently checks the large-p quadrature against the source-derived linear coefficient.

## 9. Large-m behavior of the exact first-order response: exploratory

Once the first derivative is computed directly, no large-p quadrature is needed. A preliminary high-m study of the exact finite-dimensional derivative suggests different scalings for different norms:

- the trimmed interior mean absolute response appears compatible with `m^(-3/2)`;
- the RMS response appears compatible with `m^(-1)`;
- the largest fixed-rank low-edge response appears compatible with `m^(-1/2)`;
- the high edge is much less singular;
- the full untrimmed L1 response may contain a logarithmic low-edge contribution.

These statements are **not yet accepted findings**. They were obtained in an auxiliary calculation and must be reproduced in Rust with explicit convergence tests before being promoted. In particular, they warn that the earlier finite-prime `~1/m` bulk fit at `m<=32` may be pre-asymptotic and that the limits `q -> 0` and `m -> infinity` may not be uniform.

## 10. Next implementation steps

1. Add a reusable Rust implementation of `alpha_n`, `r_n`, and the exact first-order Jacobi derivative at `q=0`.
2. Add the exact matrix derivative `K'(0)` and Rayleigh-quotient derivative of each finite crossing.
3. Add regression tests comparing the exact derivative against the independent large-prime Stieltjes result at `m=16,24,32`.
4. Use the exact derivative to push block sizes to at least `m=1024` in Rust and test the candidate norm scalings above.
5. Resolve the low-edge profile explicitly before drawing any infinite-dimensional conclusion.
6. Refactor the large exploratory `semilocal_multi_prime_probe.rs` into reusable source modules and focused tests before merge.
7. Run `cargo fmt --all -- --check`, `cargo test --all-targets`, and `cargo clippy --all-targets -- -D warnings` before any Phase 4 merge.

## 11. Interpretation boundary

The most defensible current result is:

> For fixed finite prolate compression size `m` in the single-prime semilocal model `S={infinity,p}`, the Jacobi q-series of arXiv:2403.01247 induces an ordinary analytic first-order perturbation of the finite generalized prolate matrix in `q=1/p`. The resulting source-derived first derivative agrees closely with independent large-p Stieltjes computations.

What is **not** established:

- an infinite-dimensional limiting operator theorem for this perturbation;
- a universal `1/m` correction law;
- an identification of finite crossings with zeta zeros;
- any implication proving RH.
