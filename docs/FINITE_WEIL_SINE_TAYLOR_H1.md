# Finite Weil sine Taylor H1 envelope

This note records the analytic bridge between the compact sine family and the finite-support H1 continuity envelope.

## Target and approximant

On the exact compact support `(a,b)`, use

`p(t) = sin(m*pi*t)`, `t=(rho-a)/(b-a)`,

and let `p_n` be the degree-`n` Taylor polynomial of `p` about `t=1/2`.

For derivative order `r=0,1,2,3`, Taylor's theorem on `|t-1/2| <= 1/2` gives

`sup |p^(r)-p_n^(r)| <= omega^(n+1) (1/2)^(n+1-r) / (n+1-r)!`,

with `omega=m*pi` and `n>=2`.

These are analytic remainder bounds. They are not fitted to sampled data.

## Standard bump derivative majorants

Write

`beta(t)=exp(-1/(t(1-t)))`

and `x=1/(t(1-t)) >= 4`. Each of `beta`, `beta'`, `beta''`, `beta'''` is `exp(-x)` times a polynomial in `x` and `1-2t`. Bounding `|1-2t| <= 1` and maximizing each monomial `x^k exp(-x)` on `x>=4` produces explicit finite derivative envelopes.

The implementation uses

- `B0 = sup exp(-x)`,
- `B1 <= sup x^2 exp(-x)`,
- `B2 <= 2 M2 + 2 M3 + M4`,
- `B3 <= 12 M3 + 12 M4 + 6 M5 + M6`,

where `Mk = sup_{x>=4} x^k exp(-x)`.

## Propagation through Q

Let `g=beta p` and

`D = rho d/drho`, `Q=-D^2+1/4`.

With `A=rho/(b-a)`, one has exactly

`Qg = -A g' - A^2 g'' + g/4`,

and

`D(Qg) = -(3/4) A g' - 3 A^2 g'' - A^3 g'''`.

Therefore uniform bounds on `g,g',g'',g'''` imply uniform bounds on `Qg` and `D(Qg)`. Since `d*rho=d rho/rho=d log rho`, a support of logarithmic width `L=log(b/a)` gives

`||f||_2 <= sqrt(L) sup |f|`.

The same propagation is applied to the Taylor remainder `beta(p-p_n)`. This yields upper envelopes for

- `||h||_2`,
- `||h_n||_2`,
- `||h-h_n||_2`,
- `||Dh||_2`,
- `||Dh_n||_2`,
- `||D(h-h_n)||_2`.

The approximant norms are bounded safely by target norm plus error norm.

## Connection to the finite Weil pairing

The resulting six quantities are passed directly to `H1ApproximationUpperBounds` from `weil_h1_continuity`. For a pair of sine modes this gives the deterministic finite-support perturbation envelope for the pole, archimedean, prime-power, and total compact Weil pairing pieces.

This is stronger than the preceding quadrature-only continuity probe because the norm inputs are derived from closed inequalities rather than inferred from sampled integrals.

## Numerical-arithmetic boundary

The inequalities above are mathematical. The Rust implementation evaluates their constants in binary64. It does not use outward-rounded interval arithmetic, so the returned floating-point number is not a formal machine certificate for the exact real constant. The code must therefore not label it as a proof-assistant certificate or a certified interval bound.

No result here establishes density or completeness of the Taylor family in the topology required by the full Weil criterion. No finite pairing sign, finite-dimensional positivity result, or shrinking perturbation envelope proves complete-space Weil positivity, Conjecture 4.1, or the Riemann hypothesis.
