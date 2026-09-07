# Finite Weil sine direct-Q audit

This audit cross-checks the Legendre-truncation route for the non-polynomial target generators

\[
g_m(\rho)=b(\rho)\sin(m\pi t),\qquad t=\frac{\rho-a}{b-a}.
\]

The existing truncation route first projects the sine profile into shifted Legendre polynomials and then applies the already validated images

\[
h_j=Q\bigl(b(\rho)P_j(2t-1)\bigr).
\]

The direct route instead evaluates the target under

\[
Q=-\left(\rho\frac{d}{d\rho}\right)^2+\frac14
\]

using the product identity

\[
Q(bs)=sQb-2(Db)(Ds)-bD^2s,
\qquad D=\rho\frac{d}{d\rho}.
\]

`Qb` and the compact bump value are reused from the existing compact Weil implementation. The sine derivatives are analytic.

For every declared parent dimension the audit reports, separately:

- sampled generator reconstruction residual before applying `Q`;
- maximum direct `Qg_m` amplitude;
- maximum projected `Qg_m` amplitude;
- maximum absolute direct-vs-projected `Q` residual;
- that residual normalized by the maximum direct `Q` amplitude;
- maximum critical Mellin boundary residual for the direct route;
- maximum reconstructed boundary residual for the projected route.

The sampling grid, coefficient quadrature order and boundary quadrature order are explicit experiment parameters and are retained in `DirectSineQAuditConfig`.

## Interpretation

A decreasing generator residual does not imply a decreasing `Q` residual: differentiation can amplify approximation error. The direct-Q audit makes this amplification visible rather than silently attributing the projected result to the non-polynomial target family.

The reported maxima are deterministic sampled diagnostics. They are not certified sup-norm bounds. The normalized `Q` residual is a numerical comparison scale, not a significance score or proof score.

Small direct boundary moments numerically check the expected range-of-`Q` condition for these manufactured compact generators. They do not establish the full Weil criterion.

Even if the direct and projected routes converge closely across increasing parent dimensions, this remains finite numerical evidence. It does not establish convergence of the Weil functional in a required topology, density/completeness of the target family, complete-space Weil positivity, semilocal `L^2(X_S)`, Conjecture 4.1, or the Riemann hypothesis.

A later step may use the direct compact functions inside an independently generalized mixed Weil-pairing evaluator. That stronger step should only be interpreted after this local `Q` consistency check is stable.
