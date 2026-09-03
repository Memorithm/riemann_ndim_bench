# Phase 4 — frozen-cavity soft-edge contraction asymptotic

Status: **formal local asymptotic derived from exact finite coefficients**

This note isolates the large-degree behavior of the contraction factor of the
row-frozen cavity map introduced in the finite-to-frozen resolvent programme.
It does **not** prove a uniform finite-section estimate, summability of the full
variable-coefficient transport error, or the global singular trace asymptotic.

## 1. Exact local coefficients

For parity degree

```text
d = 2 i + epsilon,   epsilon in {0,1},
```

the exact `K(0)` diagonal at an interior row is

```text
a_d
 = [2 d^2 + d + 3/4] / [2 pi (4 d + 1)].
```

The outgoing and incoming edges are respectively

```text
e_d^+
 = (1/(2 pi))
   sqrt(
     (d+1/2)(d+1)(d+3/2)(d+2)
     / [(4d+1)(4d+9)]
   ),
```

and

```text
e_d^-
 = (1/(2 pi))
   sqrt(
     (d-3/2)(d-1)(d-1/2)d
     / [(4d-7)(4d+1)]
   ).
```

The frozen symmetric edge used by `frozen_row_cavity_fixed_point` is

```text
b_d = (e_d^- + e_d^+) / 2.
```

These formulas come directly from `build_k0`; no fitted parameter enters the
derivation.

## 2. Large-degree expansion

Direct expansion in inverse powers of `d` gives

```text
a_d
 = d/(4 pi)
   + 1/(16 pi)
   + 5/(64 pi d)
   - 5/(256 pi d^2)
   + O(d^-3),
```

```text
e_d^+
 = d/(8 pi)
   + 5/(32 pi)
   + 3/(128 pi d)
   - 15/(512 pi d^2)
   + O(d^-3),
```

and

```text
e_d^-
 = d/(8 pi)
   - 3/(32 pi)
   + 3/(128 pi d)
   + 9/(512 pi d^2)
   + O(d^-3).
```

Therefore

```text
b_d
 = d/(8 pi)
   + 1/(32 pi)
   + 3/(128 pi d)
   - 3/(512 pi d^2)
   + O(d^-3),
```

and the local soft-edge gap is

```text
a_d - 2 b_d
 = 1/(32 pi d)
   - 1/(128 pi d^2)
   + O(d^-3).
```

The leading `1/d` gap is the same infrared scale that appears in the earlier
local-symbol derivation.

## 3. Positive cavity fixed point

The frozen half-line Schur map is

```text
F_d(x) = a_d - b_d^2/x.
```

Its positive stable fixed point is

```text
q_d = [a_d + sqrt(a_d^2 - 4 b_d^2)] / 2.
```

Using the expansions above,

```text
a_d^2 - 4 b_d^2
 = 1/(64 pi^2)
   - 11/(256 pi^2 d^2)
   + O(d^-3),
```

so

```text
sqrt(a_d^2 - 4 b_d^2)
 = 1/(8 pi)
   - 11/(64 pi d^2)
   + O(d^-3).
```

Consequently

```text
q_d
 = d/(8 pi)
   + 3/(32 pi)
   + 5/(128 pi d)
   - 49/(512 pi d^2)
   + O(d^-3).
```

## 4. Local contraction

The derivative of the frozen cavity map at its positive fixed point is

```text
kappa_d = F_d'(q_d) = b_d^2/q_d^2.
```

The quotient first has the expansion

```text
b_d/q_d
 = 1
   - 1/(2d)
   + 1/(4d^2)
   + O(d^-3),
```

hence

```text
kappa_d
 = 1
   - 1/d
   + 3/(4d^2)
   + O(d^-3).
```

In particular,

```text
d (1-kappa_d) -> 1
```

and

```text
d^2 (kappa_d - 1 + 1/d) -> 3/4.
```

This rules out a degree-independent contraction gap of the form
`kappa_d <= rho < 1` at the soft edge.

## 5. Formal cumulative consequence

Along a fixed parity block, `d_i = 2 i + epsilon`. Therefore

```text
kappa_{d_i}
 = 1 - 1/(2i+epsilon) + O(i^-2),
```

and formally

```text
log kappa_{d_i}
 = -1/(2i) + O(i^-2).
```

Thus a product of *frozen* contraction factors over rows has polynomial rather
than geometric decay:

```text
prod_{j=r}^s kappa_{d_j}
 = Theta((r/s)^(1/2))
```

at the formal asymptotic level, up to an endpoint-dependent constant.

This statement applies only to the frozen factors. The exact finite transport
multipliers from the cavity-error identity also contain the edge-ratio and
finite-cavity-ratio corrections exposed by `semilocal_cavity_drift`. Controlling
those corrections remains an open step.

## 6. Rust validation

`semilocal_cavity_asymptotic` exposes the displayed second-order approximants
and a diagnostic comparing

```text
1 - 1/d + 3/(4 d^2)
```

to the exact frozen contraction computed by `frozen_row_cavity_fixed_point`.
The regression tests compare the local coefficient series directly with
`build_k0` for both parity sectors and check the scaled limits above at finite
large degrees.

Those tests are finite numerical validation of the implementation and of the
algebraic expansion; they are not a substitute for a uniform remainder proof.

## 7. Remaining proof step

The next required estimate is not a uniform contraction theorem. It is a
weighted cumulative transport estimate combining:

1. the polynomial frozen contraction above;
2. the exact incoming/outgoing edge-ratio corrections;
3. the exact finite-cavity-ratio corrections;
4. the edge-asymmetry and fixed-point-variation pieces of the drift forcing;
5. explicit treatment of the finitely many initial rows where the frozen local
   symbol is not positive.

Only after these terms are controlled uniformly in the soft-edge scaling can
the finite slowly varying tridiagonal resolvent be replaced rigorously by its
local frozen model in the singular trace.

## Scientific boundary

This note concerns a local asymptotic of a finite semilocal prolate model. It
does not identify finite crossings with zeta zeros and does not prove the
Riemann hypothesis.
