# Finite Bruhat–Schwartz bridge to the semilocal E sum

Status: `SOURCE-GUIDED / EXACT FINITE LOCAL ALGEBRA / MANUFACTURED FINITE E REGRESSION`.

Primary source:

- Alain Connes, Caterina Consani, *The Scaling Hamiltonian*, §4.1, especially equation (4.6), arXiv:1910.14368.

## Source target

For the semilocal place set `S`, the source uses the unique non-zero decomposition

```text
q = u m,
u in Q_S^*,
m in M_S,
```

and writes the semilocal analogue of the Poisson map as

```text
E(f)(x) = |x|^(1/2) sum_{m in M_S} f(m x).
```

The repository already had each required finite arithmetic component separately:

- exact `Q_S` arithmetic and the unique `q = u m` decomposition;
- exact unit transport of local balls `p^k Z_p`;
- exact membership of diagonal `Q_S` samples in `p^k Z_p`;
- exact finite grouping by the representative `m in M_S`;
- the callback-driven finite `SemilocalPoissonMonoid::finite_e_sum`.

This increment connects those pieces without claiming the infinite or Hilbert-space theorem.

## Complete elementary local product

`FiniteLocalBallProduct` represents

```text
prod_{p in S_f} 1_{p^{k_p} Z_p}
```

with exactly one ball for every declared finite place.

For an explicit diagonal `q in Q_S`, membership is evaluated using the existing exact p-adic oracle:

```text
q in p^k Z_p  <=>  v_p(q) >= k.
```

For a monoid representative `m in M_S`, every finite valuation is zero by definition, so

```text
m in p^k Z_p  <=>  0 >= k.
```

No floating-point arithmetic enters these local indicators.

## Exact term-by-term descent check

For each explicit non-zero sample

```text
q = u m,
```

and each original local exponent `k_p`, the exact unit transport already implemented in the repository produces

```text
k'_p = k_p - v_p(u).
```

The bridge checks independently that

```text
prod_p 1_{p^{k_p} Z_p}(q)
==
prod_p 1_{p^{k'_p} Z_p}(m).
```

A disagreement is a hard error. This is the finite local algebra that must survive quotient descent; it is not inferred from equality of final floating-point totals.

## Manufactured finite archimedean profile

`compare_finite_bruhat_e_bridge` accepts a callback

```text
m -> a(m)
```

for the represented monoid elements. The callback is evaluated exactly once per represented `m`, and non-finite values are rejected.

For this finite fixture, unrepresented `m <= max_m` are defined to have zero contribution. This is an explicit manufactured finite support, not a proof that a general Bruhat–Schwartz function has been represented or that an infinite tail vanishes.

The direct side preserves the caller's original `Q_S` sample order:

```text
sum_q a(m(q)) * local_original(q).
```

The grouped side uses the exact orbit grouping and transported local coordinates:

```text
sum_m a(m) * sum_{q : m(q)=m} local_transported(q -> m).
```

Finally those grouped contributions are passed through

```text
SemilocalPoissonMonoid::finite_e_sum(modulus, max_m, ...),
```

which applies the source half-density factor `|x|^(1/2)` and enumerates `M_S` independently.

## Representative regression

For

```text
S_f = {2},
k_2 = 0,
```

the samples

```text
3, 6, 3/2, -12
```

all have representative `m=3` but different unit exponents. Their direct local indicators are

```text
1, 1, 0, 1.
```

After unit transport the exponents are evaluated on the same representative `m=3`; the transported indicators remain exactly

```text
1, 1, 0, 1.
```

Thus the finite local coordinates change under descent while membership is preserved term by term.

A second group at `m=5` and a nonconstant manufactured archimedean profile make the final `E` regression depend on both the monoid representative and the transported finite coordinates.

## Scientific boundary

This increment does **not**:

1. represent a general adele `x in A_S`;
2. represent a general Bruhat–Schwartz function on `A_S`;
3. construct the quotient `X_S = A_S / Q_S^*`;
4. justify rearranging an infinite `Q_S` sum;
5. prove the semilocal Poisson identity;
6. establish the general Hilbert-space map `E`;
7. prove Conjecture 4.1, Weil positivity, or RH.

It establishes only the finite exact local descent and a manufactured finite compatibility check with the existing `M_S` enumerator.

## Next Riemann-specific step

The next useful increment is to replace the arbitrary finite archimedean callback by a certified compactly supported archimedean fixture whose support bound mechanically determines `max_m`. That would remove one remaining user-supplied truncation choice while still staying strictly on the finite, auditable side of equation (4.6).
