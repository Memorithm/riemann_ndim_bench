# Compact Weil boundary audit

Status: **SOURCE-LOCKED TEST-FUNCTION REGRESSION / NOT A POSITIVITY RESULT**

This increment connects the compact archimedean fixture from the finite semilocal `E` bridge to the source boundary operator

\[
Q = -(\rho\,\partial_\rho)^2 + \frac14.
\]

For the smooth compact generator `g`, the implementation forms `f = Qg` using an analytic second derivative in the logarithmic coordinate `x = log rho`. No finite-difference derivative is used.

The existing `weil_boundary` module then evaluates the two critical Mellin moments

\[
\int_0^\infty f(\rho)\rho^{1/2}\,d^*\rho,
\qquad
\int_0^\infty f(\rho)\rho^{-1/2}\,d^*\rho.
\]

For a compactly supported smooth function in the image of `Q`, these are the two boundary conditions encoded by the source-level contract already present in the repository. The regression checks that the numerical quadrature residuals tend toward zero as the quadrature order is refined.

## What is exact

- the compact support envelope inherited from `CompactArchimedeanBump`;
- zero outside that declared support;
- the symbolic differential operator `Q` used by the repository;
- the analytic formula for the second logarithmic derivative of the manufactured bump.

## What is numerical

- evaluation of the smooth bump in binary64;
- evaluation of `Qg` in binary64;
- Gauss--Legendre quadrature of the two critical Mellin moments;
- the final residual tolerance.

The implementation guards exponential underflow near the support boundary. This is a numerical stabilization of the manufactured `C^infinity` bump, not a new analytic statement.

## Scientific boundary

Passing this audit establishes only that the manufactured compact test function behaves numerically like an element of the source boundary-condition class. It does **not** establish:

- Weil positivity;
- the semilocal trace formula;
- the Hilbert-space properties of the Poisson map `E`;
- Conjecture 4.1;
- the Riemann Hypothesis.

A subsequent increment may use this validated test-function class as input to a finite Weil-functional experiment, but positivity must remain a separate quantity to be tested rather than inferred from the vanishing boundary moments.
