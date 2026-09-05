# Scalar semilocal Poisson negative control

Status: `REFUTED ROUTE` for a tempting implementation shortcut.

## Question

After quotienting the source lattice sum by `Q_S^*`, equation (4.6) of Connes--Consani, *The Scaling Hamiltonian*, is indexed by the monoid `M_S`:

```text
E(f)(x) = |x|^(1/2) sum_{m in M_S} f(m x).
```

RiemannBench already represents `M_S`. A tempting shortcut is therefore:

1. discard the finite local coordinates entirely;
2. keep only a positive real scalar `x`;
3. replace the ordinary positive-integer sum by `M_S`;
4. expect the archimedean reciprocal Poisson identity to survive.

This document records an explicit negative control showing that this shortcut is false.

It does **not** refute the source semilocal Poisson formula, whose variable lies in the genuine semilocal quotient and whose Fourier transform uses the self-dual additive character on the full semilocal space.

## Analytic manufactured fixture

Use the Fourier convention

```text
F[e^(-pi a x^2)](xi) = a^(-1/2) e^(-pi xi^2/a).
```

Define

```text
f(x)
 = (1/3) e^(-pi x^2)
   -(4/3) e^(-4 pi x^2)
   +       e^(-9 pi x^2).
```

Then

```text
Ff(x)
 = (1/3) e^(-pi x^2)
   -(2/3) e^(-pi x^2/4)
   +(1/3) e^(-pi x^2/9).
```

The two source boundary conditions hold exactly at the coefficient level:

```text
f(0)  = 1/3 - 4/3 + 1 = 0,
Ff(0) = 1/3 - 2/3 + 1/3 = 0.
```

Thus failure of the scalar surrogate cannot be blamed on violating the archimedean Poisson boundary conditions.

## Negative control

Take

```text
S = {infinity, 2},
M_S = {positive odd integers},
x = 3/4.
```

The test compares the two scalar expressions

```text
sqrt(x)     sum_{m in M_S} Ff(m x)
```

and

```text
sqrt(x^-1)  sum_{m in M_S} f(m x^-1).
```

Both sums are truncated at `m <= 31`. Gaussian integral tails give an explicit omitted-term bound far below floating-point roundoff.

The residual is of order

```text
4.5e-2,
```

whereas the certified truncation tail is negligible and the declared roundoff allowance is `2e-14`.

The mismatch is therefore structural, not a truncation artifact.

## Interpretation

The result is:

```text
REFUTED ROUTE:
real scalar x + deletion of multiples of finite-place primes
is not the semilocal Poisson/Fourier theory.
```

This is expected mathematically. The source construction also needs:

- the finite local coordinates in `A_S`;
- a basic additive character normalized so that `Q_S` is self-dual;
- the resulting semilocal Fourier transform;
- the quotient geometry `X_S = A_S / Q_S^*`.

The index set `M_S` by itself does not encode those structures.

## Relation to existing positive control

The repository already contains an independent archimedean Poisson oracle verifying the reciprocal identity for an analytic Schwartz fixture when there are no finite places. This negative control complements that regression:

- archimedean Poisson convention: positive control;
- naive scalar finite-place deletion: negative control.

Together they protect the next implementation from accidentally interpreting prime deletion as a semilocal Fourier transform.

## Next implementation obligation

The next safe step is local rather than global: represent an exact finite-place Fourier oracle for elementary p-adic balls under the declared self-dual normalization.

Only after independent local Fourier factors exist should RiemannBench attempt a factorizable semilocal Fourier/Poisson oracle.

This test makes no claim about Conjecture 4.1, Weil positivity, or RH.
