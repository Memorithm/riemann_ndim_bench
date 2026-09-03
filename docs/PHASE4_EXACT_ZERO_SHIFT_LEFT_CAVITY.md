# Phase 4 — exact zero-shift left cavity

Status: **exact finite algebraic identity**

This note records a closed form for the left-to-right Schur cavity of the
finite semilocal `K(0)` parity blocks. Unlike the frozen-cavity asymptotic, this
identity is exact at shift zero and does not depend on the right endpoint of the
finite block.

## 1. Exact recurrence

For one parity sector write

```text
d = 2 i + epsilon,
epsilon in {0,1},
B_d = 2 pi (4d+1).
```

Let `L_i` be the left Schur denominator of the finite `K(0)` block:

```text
L_0 = K_00,
L_i = K_ii - K_i,i-1^2/L_{i-1}.
```

The exact diagonal coefficient is

```text
K_ii
 = [a_{d-1}^2 + a_d^2 + 1/4] / B_d,
```

with

```text
a_n^2 = (n+1/2)(n+1).
```

For `i>0`, the incoming edge is

```text
K_i,i-1^2
 = a_{d-2}^2 a_{d-1}^2 / (B_{d-2} B_d).
```

## 2. Closed form

The left cavity is exactly

```text
L_i
 = (d+1/2)(d+3/2) / B_d
 = (d+1/2)(d+3/2) / [2 pi (4d+1)].
```

This holds in both parity sectors, including the left endpoint.

At `d=0` (`W+`) and `d=1` (`W-`) the incoming term is absent and the formula
reduces directly to `K_00`.

## 3. Induction

Assume the formula at the preceding parity row, whose degree is `d-2`:

```text
L_{i-1}
 = (d-3/2)(d-1/2) / B_{d-2}.
```

Using

```text
a_{d-2}^2 = (d-3/2)(d-1),
a_{d-1}^2 = (d-1/2)d,
```

the incoming Schur correction simplifies exactly to

```text
K_i,i-1^2 / L_{i-1}
 = d(d-1) / B_d.
```

The exact diagonal numerator is

```text
a_{d-1}^2 + a_d^2 + 1/4
 = 2d^2 + d + 3/4.
```

Therefore

```text
L_i
 = [2d^2 + d + 3/4 - d(d-1)] / B_d
 = [d^2 + 2d + 3/4] / B_d
 = (d+1/2)(d+3/2) / B_d.
```

This closes the induction without an asymptotic approximation.

## 4. Ground-state interpretation

The exact ground-state factorisation writes the alternating-conjugated block as

```text
T_m = U0^(-1) D^T C D U0^(-1).
```

With

```text
u_i^2 = B_d A_d,
c_i = (d+1/2)(d+3/2) A_d,
```

the same cavity denominator is

```text
L_i = c_i/u_i^2.
```

Thus the closed Schur formula is the Riccati form of the exact weighted
ground-state factorisation.

## 5. Consequence for the finite-to-frozen programme

At shift zero, the left boundary contribution does not require a cumulative
Schur estimate: the entire left cavity is already known in closed form and is
independent of the right endpoint.

This does **not** make the full singular resolvent problem trivial:

- the right cavity retains the finite right boundary;
- positive shifts `t>0` no longer obey this zero-shift closed form;
- the square-root trace representation requires control over a continuum of
  shifts;
- comparison with the local frozen point `q_i` still has a nonzero local
  mismatch.

The exact identity nevertheless removes one source of uncertainty at the soft
edge and provides an anchor for perturbative-in-shift estimates.

## 6. Rust validation

`semilocal_zero_shift_cavity` exposes the exact row formula and the exact
incoming Schur correction.

The regression suite checks:

1. the closed form against `cavity_green_bands(0)` over complete finite blocks;
2. the induction identity `edge^2/L_{i-1} = d(d-1)/B_d`;
3. reconstruction of the exact `K(0)` diagonal;
4. independence of the left-cavity prefix from the right endpoint.

These numerical checks validate the implementation; the identity itself is the
finite algebraic induction displayed above.

## Scientific boundary

This result is an exact identity for one finite semilocal cavity at zero shift.
It does not prove the global singular trace asymptotic, does not identify finite
crossings with zeta zeros, and does not prove the Riemann hypothesis.
