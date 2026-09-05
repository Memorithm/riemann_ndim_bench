# Semilocal Poisson Fourier multiplier

Status: source-locked Riemann-specific bridge component.

Primary source:

- Alain Connes, Caterina Consani, *The Scaling Hamiltonian*, arXiv:1910.14368, §4.1, especially equations (4.5)–(4.6) and the discussion immediately following them.

## 1. Source map

For a finite semilocal place set `S`, let `M_S` be the monoid of positive integers prime to every finite prime in `S`. The source introduces

```text
E(f)(x) = |x|^(1/2) sum_{m in M_S} f(m x).
```

RiemannBench represents the exact arithmetic monoid and finite support-truncated sum in `src/semilocal_poisson.rs`.

## 2. Fourier reading in the source

With multiplicative Fourier parameter `t`, the source writes formally

```text
sum_{m in M_S} m^(-1/2 + i t)
  = zeta(1/2 - i t)
    product_{p in S_f} (1 - p^(-1/2 + i t)).
```

Equivalently, if

```text
s = 1/2 - i t,
```

then the finite-place multiplier is

```text
zeta(s) product_{p in S_f} (1 - p^(-s)).
```

The finite product deletes the Euler factors belonging to the finite places already included in the semilocal space.

## 3. Critical scientific boundary

The cited paper explicitly presents this Fourier discussion at a formal/heuristic level on the critical line. The Dirichlet series

```text
sum m^(-1/2 + i t)
```

is not absolutely convergent there.

Therefore the implementation deliberately separates:

- `EXACT FINITE ALGEBRA`: `product_{p in S_f}(1-p^(-s))` for finite `S`;
- `EXACT CONVERGENT DOMAIN`: direct Dirichlet-series checks for real `s=sigma>1`;
- `SOURCE FORMAL`: substitution `s=1/2-it` on the critical line;
- `NOT CLAIMED`: convergence of the critical-line Dirichlet series or a proof of the Poisson/Fourier intertwining on the full semilocal Hilbert space.

No code path silently upgrades the source's formal critical-line manipulation to a theorem.

## 4. Executable convergent-domain oracle

For real `sigma>1`, absolute convergence gives

```text
sum_{m in M_S} m^(-sigma)
  = zeta(sigma) product_{p in S_f} (1-p^(-sigma)).
```

The bench checks the `sigma=2` instances

```text
S={infinity}:       pi^2/6,
S={infinity,2}:     pi^2/8,
S={infinity,2,3}:   pi^2/9.
```

A finite prefix through `N` carries the elementary certified tail bound

```text
sum_{m in M_S, m>N} m^(-sigma)
 <= sum_{n>N} n^(-sigma)
 <= N^(1-sigma)/(sigma-1).
```

Thus these regressions do not depend on a fitted constant or on a critical-line continuation.

## 5. Role in the RH dependency graph

This module closes only the arithmetic/Fourier-multiplier bookkeeping around `E`:

```text
finite place set S
  -> M_S
  -> E map
  -> deleted Euler factors in multiplicative Fourier variables.
```

It does not yet construct:

- a general adele `x in A_S`;
- the quotient Hilbert space `L^2(X_S)`;
- the semilocal Fourier transform `F_alpha`;
- the source Poisson identity as an operator equality;
- the Sonin projection or the semilocal Weil positivity estimate;
- Conjecture 4.1.

Those remain separate proof dependencies.

Nothing in this document or module constitutes a proof of the Riemann hypothesis.
