# Certified archimedean Poisson oracle

Status: manufactured regression for a source theorem used by the Riemann-specific bridge.

Primary source:

- Alain Connes, Caterina Consani, *The Scaling Hamiltonian*, arXiv:1910.14368, §4.1, equations (4.2)–(4.3).

## Source identity

For an even Schwartz function satisfying

```text
f(0) = 0,
Ff(0) = 0,
```

and the Fourier convention

```text
Ff(xi) = integral_R f(x) exp(-2 pi i x xi) dx,
```

the source writes

```text
E(f)(x) = x^(1/2) sum_{n>=1} f(n x),
```

and Poisson summation gives

```text
E(Ff)(x) = E(f)(x^-1),  x>0.
```

This identity is the archimedean prototype for the semilocal map introduced later in equation (4.6).

## Manufactured analytic fixture

The bench uses

```text
f(x)
 = (1/3) exp(-pi x^2)
 - (4/3) exp(-4 pi x^2)
 +       exp(-9 pi x^2).
```

The two source boundary conditions hold algebraically:

```text
f(0) = 1/3 - 4/3 + 1 = 0.
```

For

```text
F[exp(-pi a x^2)](xi)
 = a^(-1/2) exp(-pi xi^2/a),
```

one also has

```text
Ff(0)
 = 1/3 - (4/3)/2 + 1/3
 = 0.
```

Thus no numerical Fourier transform is used in the oracle.

## Certified truncation

Each E-sum is truncated at a finite integer `N`. For a Gaussian component

```text
c exp(-pi a (n x)^2),
```

monotonicity plus the elementary Gaussian integral estimate gives

```text
sum_{n>N} |c| exp(-pi a x^2 n^2)
 <= |c| exp(-pi a x^2 N^2) / (2 pi a x^2 N).
```

After multiplication by `sqrt(x)`, this supplies an explicit absolute tail bound for each finite E-sum.

The integration tests compare

```text
E(Ff)(x)
```

and

```text
E(f)(1/x)
```

for several reciprocal scales and require the residual to lie inside the two certified truncation bounds plus a declared floating-point roundoff allowance.

## Proof status

- Poisson summation / equation (4.3): `SOURCE THEOREM`.
- Fourier transform of each Gaussian: exact analytic identity used as a manufactured oracle.
- Finite E-sums: numerical values with explicit rigorous truncation bounds.
- Floating arithmetic equality: regression under a declared roundoff allowance.

The test does not re-prove Poisson summation in general. Its purpose is to verify that the bench's conventions, half-density factor, reciprocal scaling and Fourier normalization agree with the source on an independently controlled function.

## Riemann-specific role

This closes a normalization checkpoint between the archimedean E-map and the semilocal Poisson program. It does not establish the semilocal Poisson identity on `A_S` or `X_S`, does not prove Conjecture 4.1, and does not prove RH.
