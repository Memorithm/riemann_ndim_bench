# Phase 4 — exact first shift response of the left cavity

Status: **exact finite first-derivative identity at zero shift**

This note differentiates the exact left Schur cavity of the semilocal `K(0)`
parity blocks with respect to the additive resolvent shift. The result is exact
at `t=0` for every finite row and both parity sectors.

It is a first-order statement only. In particular, it does not assert that the
same multiplicative law holds at finite `t`, nor that the Taylor remainder is
uniform in the row index.

## 1. Shifted left cavity

For one finite parity block, define

```text
L_0(t) = K_00 + t,
```

and for `i>0`

```text
L_i(t)
 = K_ii + t - K_i,i-1^2/L_{i-1}(t).
```

At zero shift the preceding exact result gives, for parity degree

```text
d = 2 i + epsilon,
```

```text
L_i(0)
 = (d+1/2)(d+3/2) / [2 pi (4d+1)].
```

The incoming Schur correction is also exact:

```text
C_i
 := K_i,i-1^2/L_{i-1}(0)
 = d(d-1) / [2 pi (4d+1)].
```

For the left endpoint, `C_0=0` in the corresponding parity degree.

## 2. Constant Riccati gap

Subtracting the incoming correction from the closed left denominator gives

```text
L_i(0) - C_i
 = [(d+1/2)(d+3/2)-d(d-1)] / [2 pi (4d+1)]
 = [3d+3/4] / [2 pi (4d+1)]
 = 3/(8 pi).
```

Thus the combination that closes the differentiated recurrence is exactly
independent of row and parity:

```text
L_i(0) - C_i = 3/(8 pi).
```

## 3. Differentiated Schur recurrence

Let

```text
D_i = L_i'(0).
```

Differentiating the finite recurrence gives

```text
D_0 = 1,
```

and

```text
D_i
 = 1
   + [K_i,i-1^2/L_{i-1}(0)^2] D_{i-1}.
```

Set

```text
c = 8 pi/3.
```

At either parity left endpoint,

```text
L_0(0)=3/(8 pi),
```

so

```text
D_0=1=c L_0(0).
```

Assume `D_{i-1}=c L_{i-1}(0)`. Then

```text
D_i
 = 1 + c K_i,i-1^2/L_{i-1}(0)
 = 1 + c C_i.
```

The constant-gap identity gives

```text
1
 = c [L_i(0)-C_i],
```

hence

```text
D_i = c L_i(0).
```

This closes the induction.

Therefore, exactly,

```text
L_i'(0) = (8 pi/3) L_i(0)
```

for every row and both parity sectors.

Equivalently,

```text
d/dt log L_i(t) |_(t=0) = 8 pi/3.
```

Using the degree form of `L_i(0)`, the derivative itself is

```text
L_i'(0)
 = (2d+1)(2d+3) / [3(4d+1)].
```

## 4. Meaning for the soft-edge programme

This identity strengthens the zero-shift anchor in two ways:

1. the left cavity itself is exact at `t=0`;
2. its **relative first shift response is also exact and row-independent**.

Thus, for each fixed row,

```text
L_i(t)
 = L_i(0) [1 + (8 pi/3)t] + O_i(t^2).
```

The subscript on `O_i(t^2)` is essential. The present result does not show that
the second-order remainder is bounded uniformly as `i -> infinity`.
Consequently it does not by itself justify replacing the full positive-shift
left cavity by a common multiplicative deformation over the whole soft-edge
region.

The next question is the size and row dependence of the second and higher shift
derivatives, or an alternative monotone comparison valid on the natural
soft-edge shift scale.

## 5. Rust validation

`semilocal_zero_shift_response` exposes:

- the exact row response;
- the exact shift derivative;
- the universal logarithmic derivative `8 pi/3`;
- the derivative prefix over a finite parity block.

The regression suite differentiates the Schur recurrence directly using the
exact `K(0)` edges and exact zero-shift left denominators, then compares the
result row by row with the closed derivative. It also validates the constant
Riccati gap and the direct degree formula.

These tests validate the implementation. The mathematical identity is the
finite induction displayed above.

## Scientific boundary

This result controls one first derivative of one finite semilocal cavity. It
does not prove a uniform positive-shift resolvent estimate, the global singular
trace asymptotic, an identification with zeta zeros, or the Riemann hypothesis.
