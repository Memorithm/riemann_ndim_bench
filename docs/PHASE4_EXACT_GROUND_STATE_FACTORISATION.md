# Exact ground-state factorisation of the archimedean parity blocks

This note records an exact finite-dimensional structure underlying the soft-edge asymptotics. Unlike the local-symbol argument, the identities below are algebraic and hold before taking any large-block limit.

## 1. Archimedean Jacobi recurrence

At `q=0` the orthonormal Jacobi coefficients are

```text
a_n^2 = (n + 1/2)(n + 1).
```

Define the positive central-binomial sequence

```text
A_n = 4^(-n) binom(2n,n).
```

It satisfies the exact recurrence

```text
A_{n+1}/A_n = (n + 1/2)/(n + 1).
```

Set

```text
r_n = sqrt(A_n),
y_n = i^n r_n.
```

Then

```text
a_n r_{n+1} = (n + 1/2) r_n,
a_{n-1} r_{n-1} = n r_n.
```

For the zero-diagonal Jacobi matrix `J` with off-diagonal entries `a_n`, this gives exactly

```text
(J y)_n
 = a_{n-1} y_{n-1} + a_n y_{n+1}
 = (i/2) y_n.
```

Hence

```text
(J^2 + 1/4) y = 0.
```

No asymptotic approximation is used here.

## 2. Parity restriction and generalized-prolate normalization

For one parity block write

```text
d = 2i + epsilon,
epsilon in {0,1}.
```

The generalized-prolate normalization is

```text
B_d = 2*pi*(4d+1).
```

At `lambda=0`, the positive generalized crossing matrix is

```text
K = B^(-1/2) (J^2 + 1/4) B^(-1/2)
```

restricted to the chosen parity.

Since

```text
y_d = phase_epsilon * (-1)^i sqrt(A_d),
```

the alternating conjugation `(Ux)_i=(-1)^i x_i` turns the zero solution into

```text
u_i = sqrt(B_d A_d).
```

Thus, for the infinite parity Jacobi matrix `T = U K U`, whose off-diagonal entries are negative,

```text
T u = 0
```

exactly.

## 3. Exact conductance factorisation

Let `b_i>0` denote the magnitude of the `T` off-diagonal joining parity degrees `d` and `d+2`:

```text
b_i = a_d a_{d+1} / sqrt(B_d B_{d+2}).
```

Define

```text
c_i = b_i u_i u_{i+1}.
```

Using the recurrence for `A_n`, this simplifies exactly to

```text
c_i = (d + 1/2)(d + 3/2) A_d.
```

Let `T_m` be the principal block on indices `0,...,m-1`. For any vector `f`, set

```text
g_i = f_i/u_i,
g_m = 0.
```

The ground-state identity gives the exact finite quadratic form

```text
<f,T_m f>
 = sum_{i=0}^{m-1} c_i |g_i-g_{i+1}|^2.
```

Equivalently, if `D` is the upper first-difference matrix `(Dg)_i=g_i-g_{i+1}` with `g_m=0`, `U0=diag(u_i)`, and `C=diag(c_i)`, then

```text
T_m = U0^(-1) D^T C D U0^(-1).
```

This is an exact weighted discrete Sturm-Liouville / ground-state factorisation.

## 4. Exact inverse kernel

The triangular difference matrix is explicitly invertible:

```text
(D^(-1)h)_i = sum_{k=i}^{m-1} h_k.
```

Therefore

```text
T_m^(-1)
 = U0 D^(-1) C^(-1) (D^T)^(-1) U0.
```

Its entries are

```text
(T_m^(-1))_{ij}
 = u_i u_j * sum_{k=max(i,j)}^{m-1} 1/c_k.
```

The inverse of the original positive-off-diagonal block is obtained by restoring the alternating signs:

```text
K_m^(-1) = U T_m^(-1) U.
```

Thus the absolute inverse kernel is known exactly at every finite block size.

## 5. Rust validation of the exact identities

The Rust test `validates_exact_ground_state_factorization_and_inverse_kernel` checks the finite identities directly.

For a `64 x 64` principal block it gives

```text
W+ worst zero-mode residual                  2.13e-14
W+ worst conductance relative error          4.64e-16
W+ worst |T_m G - I|                         2.49e-14

W- worst zero-mode residual                  2.84e-14
W- worst conductance relative error          4.10e-16
W- worst |T_m G - I|                         1.78e-14
```

These values are consistent with the identities above up to floating-point roundoff.

## 6. Large-index weights and intrinsic length

The central-binomial asymptotic gives

```text
A_d ~ 1/sqrt(pi*d).
```

Hence

```text
u_i^2 = B_d A_d ~ 8*sqrt(pi)*sqrt(d),
c_i ~ d^(3/2)/sqrt(pi).
```

The exact ratio relevant to the Liouville coordinate is independent of `A_d`:

```text
u_i^2/c_i = B_d / [(d+1/2)(d+3/2)].
```

Therefore

```text
sqrt(u_i^2/c_i) ~ sqrt(8*pi/d) ~ 2*sqrt(pi/i),
```

and

```text
L_m = sum_{i<m} sqrt(u_i^2/c_i)
    ~ 4*sqrt(pi*m).
```

At `m=16384` the Rust test gives

```text
W+ : L_m/[4 sqrt(pi m)] = 0.9939072509788144
W- : L_m/[4 sqrt(pi m)] = 0.9924381195814174.
```

The remaining relative deficit is of order `m^(-1/2)`, compatible with an additive boundary correction to `L_m`.

## 7. Continuum critical-Bessel interpretation

Keeping only the leading row coefficients of the alternating block gives

```text
K_cont
 = -(1/(4*pi)) d/dx [x d/dx]
   + 1/(64*pi*x).
```

With `x=r^2`, this becomes

```text
K_cont
 = (1/(16*pi))
   [-d^2/dr^2 - (1/r)d/dr + 1/(4r^2)].
```

On the radial measure `r dr`, the standard unitary multiplication by `sqrt(r)` converts

```text
-d^2/dr^2 - (1/r)d/dr + nu^2/r^2
```

to

```text
-d^2/dr^2 + (nu^2-1/4)/r^2.
```

Here `nu=1/2`, so the inverse-square term cancels exactly at leading order and the Liouville model is the free one-dimensional Laplacian, up to the overall factor `1/(16*pi)`.

The continuum discussion is asymptotic; the ground-state factorisation and finite inverse formula above are exact.

## 8. Exact Hardy/Copson representation

Define

```text
R_m = U0 D^(-1) C^(-1/2).
```

Then

```text
T_m^(-1) = R_m R_m^T.
```

The entries are explicitly

```text
(R_m)_{ik} = u_i / sqrt(c_k),  k>=i,
```

and zero otherwise. Asymptotically this is a weighted Hardy/Copson kernel

```text
i^(1/4) k^(-3/4) 1_{k>=i}.
```

In the continuum variable `r=sqrt(i)`, this becomes formally

```text
(R f)(r)
 ~ 4 sqrt(pi) r integral_r^R f(s)/s ds.
```

Removing the scalar factor gives the limiting Hardy operator

```text
(C f)(r) = r integral_r^R f(s)/s ds,
```

with formal inverse

```text
C^(-1) = -d/dr + 1/r.
```

The associated principal second-order factor satisfies

```text
(C C^*)^(-1)
 = (d/dr + 1/r)(-d/dr + 1/r)
 = -d^2/dr^2.
```

Thus the critical free-Laplacian cancellation is visible directly from the exact inverse factorisation.

## 9. Spectral validation and correction of the boundary phase

A dedicated Rust test checks the low square-root spectrum through `m=16384`.

The scaled spacings

```text
(lambda_{j+1}-lambda_j) L_m/pi
```

converge to 1 for both parity blocks. At `m=16384`, the mean over the first eight gaps is

```text
W+ : 0.9939026771501716
W- : 0.9924185037937292.
```

The initially guessed phase `j+1/2` is **not** supported. Instead the data give

```text
lambda_j L_m/pi -> j+1.
```

At `m=16384`:

```text
W+ : lambda_0 L_m/pi = 0.9939034484480231
W- : lambda_0 L_m/pi = 0.9924191935620239.
```

Therefore one must not assign a regular Neumann condition at the transformed left endpoint merely from the absence of a left finite-difference term in the `g`-form. The left endpoint is singular before the Liouville transform, and the admissible transformed phase is selected by that singular structure.

A sharper finite-size observation is that dividing the measured mean spacing by `L_m/[4 sqrt(pi m)]` gives

```text
W+ : 0.9999953981333386
W- : 0.9999802347498538
```

at `m=16384`. Likewise the first scaled eigenvalue divided by the same length ratio gives

```text
W+ : 0.9999961741593217
W- : 0.9999809297738368.
```

Thus the cleaner leading law is currently

```text
lambda_j * 4 sqrt(pi m)/pi -> j+1.
```

The difference between `L_m` and `4 sqrt(pi m)` is an additive boundary correction and is largely compensated by the spectral endpoint phase.

Full numerical output is recorded in `PHASE4_FREE_LAPLACIAN_SPECTRAL_VALIDATION_2026-08-15.md`.

## 10. Relation to soft-edge Jacobi theory

The exact Hardy representation places the remaining singular trace problem close to established unbounded-Jacobi methods. A relevant primary reference is Grzegorz Świderski, *Periodic perturbations of unbounded Jacobi matrices III: The soft edge regime*, arXiv:1707.06486. Its theorems are not claimed to prove the weighted `K^(-1/2)` trace asymptotic studied here.

## 11. Next rigorous target

The next target is to combine the exact representation

```text
T_m^(-1)=R_m R_m^T
```

with the validated free-Laplacian low-spectrum scale to control

```text
Tr(T_m^(-1/2) H_m)
```

uniformly at the soft edge.

The exact identities remove the need to treat `T_m` as an opaque dense matrix. What remains is a singular weighted trace estimate for a concrete Hardy operator with explicit finite-section weights.

## Scientific boundary

All exact identities in this note concern finite archimedean generalized-prolate blocks. The continuum and spectral-limit interpretations are asymptotic statements about this model. They do not identify finite compression eigenvalues with Riemann-zeta zeros and do not imply the Riemann hypothesis.
