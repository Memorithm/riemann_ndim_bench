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
A_{n+1}/A_n
 = (n + 1/2)/(n + 1).
```

Set

```text
r_n = sqrt(A_n),
y_n = i^n r_n.
```

Then

```text
a_n r_{n+1}
 = (n + 1/2) r_n,

a_{n-1} r_{n-1}
 = n r_n.
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
K
 = B^(-1/2) (J^2 + 1/4) B^(-1/2)
```

restricted to the chosen parity.

Since

```text
y_d
 = phase_epsilon * (-1)^i sqrt(A_d),
```

the alternating conjugation

```text
(Ux)_i = (-1)^i x_i
```

turns the zero solution into the positive sequence

```text
u_i
 = sqrt(B_d A_d).
```

Thus, for the infinite parity Jacobi matrix

```text
T = U K U,
```

whose off-diagonal entries are negative,

```text
T u = 0
```

exactly.

## 3. Exact conductance factorisation

Let `b_i>0` denote the magnitude of the `T` off-diagonal joining parity degrees `d` and `d+2`:

```text
b_i
 = a_d a_{d+1}
   / sqrt(B_d B_{d+2}).
```

Define the conductance

```text
c_i = b_i u_i u_{i+1}.
```

Using the recurrence for `A_n`, this simplifies exactly to

```text
c_i
 = (d + 1/2)(d + 3/2) A_d.
```

Let `T_m` be the principal block on indices `0,...,m-1`. For any vector `f`, set

```text
g_i = f_i/u_i,
g_m = 0.
```

The ground-state identity `T u=0` gives the exact finite quadratic form

```text
<f,T_m f>
 = sum_{i=0}^{m-1}
     c_i |g_i-g_{i+1}|^2.
```

Equivalently, if `D` is the upper first-difference matrix

```text
(Dg)_i = g_i-g_{i+1},
g_m=0,
```

and `U0=diag(u_i)`, `C=diag(c_i)`, then

```text
T_m
 = U0^(-1) D^T C D U0^(-1).
```

This is an exact weighted discrete Sturm-Liouville / ground-state factorisation.

## 4. Exact inverse kernel

The triangular difference matrix is explicitly invertible:

```text
(D^(-1)h)_i
 = sum_{k=i}^{m-1} h_k.
```

Therefore

```text
T_m^(-1)
 = U0 D^(-1) C^(-1) (D^T)^(-1) U0.
```

Its entries are

```text
(T_m^(-1))_{ij}
 = u_i u_j
   * sum_{k=max(i,j)}^{m-1} 1/c_k.
```

The inverse of the original positive-off-diagonal block is obtained by restoring the alternating signs:

```text
K_m^(-1)
 = U T_m^(-1) U.
```

Thus the absolute inverse kernel is known exactly at every finite block size.

## 5. Large-index weights

The central-binomial asymptotic gives

```text
A_d
 ~ 1/sqrt(pi*d).
```

Hence

```text
u_i^2
 = B_d A_d
 ~ 8*sqrt(pi)*sqrt(d),

c_i
 = (d+1/2)(d+3/2) A_d
 ~ d^(3/2)/sqrt(pi).
```

The exact ratio relevant to the Liouville coordinate is independent of `A_d`:

```text
u_i^2/c_i
 = B_d / [(d+1/2)(d+3/2)].
```

Therefore

```text
sqrt(u_i^2/c_i)
 ~ sqrt(8*pi/d)
 ~ 2*sqrt(pi/i).
```

The accumulated intrinsic length through row `m` is consequently

```text
L_m
 = sum_{i<m} sqrt(u_i^2/c_i)
 ~ 4*sqrt(pi*m).
```

This is the discrete counterpart of the `r=sqrt(i)` soft-edge coordinate.

## 6. Continuum critical-Bessel interpretation

Keeping only the leading row coefficients of the alternating block gives

```text
K_cont
 = -(1/(4*pi)) d/dx [x d/dx]
   + 1/(64*pi*x).
```

With

```text
x = r^2,
```

this becomes

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

This explains structurally why the soft edge is critical rather than gapped and why the natural finite-section length scales like `sqrt(m)`.

The continuum discussion is asymptotic; the ground-state factorisation and finite inverse formula above are exact.

## 7. Why this matters for the proof route

The exact factorisation rewrites the inverse square-root problem in terms of a weighted first-difference / Hardy-type operator rather than an opaque dense spectral function.

Define

```text
R_m
 = U0 D^(-1) C^(-1/2).
```

Then

```text
T_m^(-1) = R_m R_m^T.
```

The entries of `R_m` are explicitly

```text
(R_m)_{ik}
 = u_i / sqrt(c_k)
```

for `k>=i`, and zero otherwise.

Asymptotically,

```text
u_i ~ const * i^(1/4),
sqrt(c_k) ~ const * k^(3/4),
```

so `R_m` is a concrete weighted Hardy/Copson triangular operator with leading kernel

```text
i^(1/4) k^(-3/4)  1_{k>=i}.
```

This places the remaining singular trace problem much closer to the established theory of Hardy operators and inverses of soft-edge Jacobi matrices.

A relevant primary reference is Grzegorz Świderski, *Periodic perturbations of unbounded Jacobi matrices III: The soft edge regime*, arXiv:1707.06486. That work derives explicit inverse/resolvent formulas for soft-edge Jacobi matrices and reduces boundedness/compactness questions to Hardy-type estimates. Its theorems do not directly prove the weighted `K^(-1/2)` trace asymptotic studied here.

## 8. Next rigorous target

The next target is to exploit

```text
T_m^(-1) = R_m R_m^T
```

and

```text
T_m^(-1/2) = |R_m^T|
```

in a suitable polar-decomposition sense, then control the weighted trace

```text
Tr(T_m^(-1/2) H_m)
```

through the explicit Hardy operator `R_m` and the exact row asymptotics of `H_m`.

The immediate finite checks are:

1. verify `T u=0` from the Rust coefficient generator for both parities;
2. verify the conductance formula `c_i=(d+1/2)(d+3/2)A_d`;
3. compare the closed inverse kernel above with a dense inverse for small blocks;
4. measure the intrinsic length `L_m/(4*sqrt(pi*m)) -> 1`.

## Scientific boundary

These exact identities concern the finite archimedean generalized-prolate blocks and their semilocal first-order analysis. They do not identify finite compression crossings with zeta zeros and do not imply the Riemann hypothesis.
