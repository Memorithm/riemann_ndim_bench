# Formal origin of the `(log m)^2` law

This note gives a **formal asymptotic derivation** of the numerical candidate

```text
S(m) ~ (1/(2*pi^2)) (log m)^2 + O(log m)
```

for the total first-order semilocal prolate response.

It is not yet a theorem. The exact finite-dimensional formulas used below are established in the bench; the passage from a slowly varying tridiagonal matrix to a local trace symbol, including uniform control at the soft spectral edge, still needs a proof.

## 1. Exact finite-dimensional trace identity

For either parity block, let

```text
K(q) = K(0) + q K'(0) + O(q^2)
```

with `K(0)>0` at finite block size. Then

```text
d/dq Tr sqrt(K(q)) |_{q=0}
  = (1/2) Tr(K(0)^(-1/2) K'(0)).
```

The finite-block sign lemma gives

```text
K'_+(0) < 0,
K'_-(0) > 0
```

as quadratic forms. Therefore, with

```text
H_+ = -K'_+(0),
H_- =  K'_-(0),
```

both `H_+` and `H_-` are positive and

```text
S(m)
 = (1/2) Tr(K_+^(-1/2) H_+)
   + (1/2) Tr(K_-^(-1/2) H_-).
```

The two parity blocks have the same leading large-index coefficients, so it is enough to derive one block and double the result.

## 2. Archimedean Jacobi coefficients

At `q=0`,

```text
a_n^2 = (n + 1/2)(n + 1).
```

For parity degree `d`,

```text
B_d = 2*pi*(4d+1),
```

and

```text
K_ii
  = [a_{d-1}^2 + a_d^2 + 1/4] / B_d,

K_i,i+1
  = a_d a_{d+1} / sqrt(B_d B_{d+2}).
```

The parity degree is

```text
d = 2i + epsilon,
epsilon in {0,1}.
```

The large-`d` expansions are

```text
K_ii
 = d/(4*pi)
   + 1/(16*pi)
   + 5/(64*pi*d)
   + O(d^-2),

K_i,i+1
 = d/(8*pi)
   + 5/(32*pi)
   + 3/(128*pi*d)
   + O(d^-2).
```

The preceding edge coefficient satisfies

```text
K_i,i-1
 = d/(8*pi)
   - 3/(32*pi)
   + 3/(128*pi*d)
   + O(d^-2).
```

Hence the soft-edge residual is

```text
K_ii - K_i,i+1 - K_i,i-1
 = 1/(32*pi*d) + O(d^-2)
 = 1/(64*pi*i) + O(i^-2).
```

This `1/i` term is crucial: it supplies the natural infrared scale that regularizes the local square-root singularity.

## 3. Alternating conjugation and soft-edge form

Let

```text
(Ux)_i = (-1)^i x_i.
```

The positive off-diagonal entries of `K` become negative after conjugation. To leading order the resulting operator has the discrete Sturm-Liouville form

```text
U K U^*
  = - nabla^* (b_i nabla) + V_i + lower-order terms,
```

with

```text
b_i = i/(4*pi) + O(1),
V_i = 1/(64*pi*i) + O(i^-2).
```

Freezing the coefficients near a large index `i` gives the local soft-edge symbol

```text
k(i,theta)
  = 2 b_i (1 - cos theta) + V_i + ...
```

and therefore, for small `theta`,

```text
k(i,theta)
  = i*theta^2/(4*pi)
    + 1/(64*pi*i)
    + lower-order terms

  = i/(4*pi)
    * [theta^2 + 1/(16 i^2)]
    + lower-order terms.
```

Thus the soft edge is quadratic in `theta`, with a built-in scale

```text
theta_IR ~ 1/(4i).
```

No ad hoc numerical cutoff is required at leading order.

## 4. First-order perturbation at the same edge

The exact closed form is

```text
K'_ii(0)
  = -3 alpha_d / (2 sqrt(2) pi),
```

with

```text
|alpha_d| ~ 1/sqrt(pi*d).
```

For the sign-corrected matrices `H_+` and `H_-`, let `D_i` be the diagonal magnitude and `O_i` the off-diagonal magnitude. The exact ratio satisfies

```text
O_i / D_i -> 1/6.
```

After the same alternating conjugation, the soft-edge multiplication part is

```text
D_i - O_i - O_{i-1}.
```

Using `d=2i+epsilon`, its leading term is

```text
h(i,0)
 = 1/(2*pi^(3/2)*sqrt(i))
   + O(i^(-3/2)).
```

The gradient part of `H` is lower order in the soft-edge singular integral, so formally

```text
h(i,theta)
 = 1/(2*pi^(3/2)*sqrt(i))
   + O(i^(-1/2)*theta^2)
   + O(i^(-3/2)).
```

## 5. Local trace density

For one parity block,

```text
T_m
 = (1/2) Tr(K^(-1/2) H).
```

The formal local-symbol trace density is

```text
(1/2) * (1/(2*pi))
* h(i,theta) / sqrt(k(i,theta)).
```

Substituting the leading soft-edge forms gives

```text
h(i,theta) / sqrt(k(i,theta))
 ~ 1/(pi*i)
    * 1/sqrt(theta^2 + 1/(16 i^2)).
```

Therefore

```text
local density
 ~ 1/(4*pi^2*i)
    * integral dtheta
      / sqrt(theta^2 + 1/(16 i^2)).
```

The symmetric soft-edge integral has the leading behavior

```text
integral_{-theta_0}^{theta_0}
  dtheta / sqrt(theta^2 + 1/(16 i^2))
 = 2 log i + O(1),
```

for any fixed small `theta_0>0`.

Hence one parity contributes

```text
T_m
 ~ sum_{i<=m}
      [log i]/[2*pi^2*i]

 ~ 1/(4*pi^2) (log m)^2
    + O(log m).
```

There are two parity blocks with the same leading coefficient. Consequently the formal prediction is

```text
S(m)
 ~ 1/(2*pi^2) (log m)^2
   + O(log m).
```

This is exactly the numerical candidate approached by the centered dyadic coefficient `A_m`.

## 6. Comparison with the Rust checkpoints

The homogeneous centered coefficients are

```text
A_256  = 0.0511032354345484
A_512  = 0.0509084935580483
A_1024 = 0.0507978750345720
A_2048 = 0.0507359112365777
A_4096 = 0.0507022850427297
A_8192 = 0.0506827887066639
```

while

```text
1/(2*pi^2)
 = 0.0506605918211689...
```

The finite-size gap remains resolved numerically and should not be suppressed by fitting. The quantity

```text
m * [A_m - 1/(2*pi^2)] / log(m)
```

stays near `0.0202` over the validated range, suggesting a next correction compatible with

```text
A_m
 = 1/(2*pi^2)
   + [c log m + d]/m
   + ...
```

but this next term is not derived here.

## 7. What remains to prove

The coefficient `1/(2*pi^2)` above follows formally from exact leading coefficients, but a proof requires uniform control of the singular trace near `theta=0`. In particular one needs to justify:

1. replacement of the finite slowly varying tridiagonal matrices by the local symbol in a trace involving the singular function `x -> x^(-1/2)`;
2. uniform treatment of the soft-edge region where `k(i,theta)` is of order `1/i`;
3. summability of the errors after the `i`-sum;
4. parity-boundary contributions and the first finitely many rows;
5. the next-order term responsible for the observed `log(m)/m` correction in `A_m`.

A generalized locally Toeplitz / discrete Sturm-Liouville approach is a plausible framework, but ordinary bounded continuous test-function spectral distribution results are not by themselves sufficient because `x^(-1/2)` is singular at the soft edge.

## Scientific boundary

This derivation concerns only the large-block asymptotics of a finite semilocal first-order prolate perturbation. Even if the asymptotic coefficient is proved, that would not identify the finite crossings with zeta zeros and would not prove the Riemann hypothesis.
