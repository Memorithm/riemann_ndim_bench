# Finite Weil direct sine pairing audit

This audit extends the local direct-`Q` check to the complete finite compact Riemann--Weil pairing used by the repository.

For

\[
g_m(\rho)=b(\rho)\sin(m\pi t),\qquad t=\frac{\rho-a}{b-a},
\]

the direct route evaluates

\[
h_m^{\mathrm{direct}}=Qg_m
\]

with the analytic product-rule implementation from `weil_sine_direct_q`, then evaluates mixed pairings

\[
\psi\!\left((h_i^{\mathrm{direct}})^* * h_j^{\mathrm{direct}}\right)
\]

through `weil_compact_pairing`.

The projected route computes the shifted-Legendre coefficients of the sine profile at a declared parent dimension, evaluates the resulting compact linear combination through the same generic pairing engine, and independently checks its total against

\[
C^T A C,
\]

where `A` is the existing finite Legendre Weil pairing matrix.

For every parent dimension the audit retains:

- maximum pairing amplitude;
- maximum direct-vs-projected pairing residual;
- the same residual normalized by the observed pairing amplitude;
- maximum disagreement between the projected generic route and `C^T A C`;
- separate direct-vs-projected residuals for the critical pole term, the archimedean term and the finite prime-power total;
- raw forward/reverse pairing asymmetry for both routes;
- maximum critical Mellin boundary residual for both routes.

The generic compact pairing implementation is cross-regressed against the existing Legendre matrix normalization. The matrix comparison is therefore a compatibility check on the new abstraction as well as a control on the sine projection.

## Interpretation boundary

The shifted-Legendre approximation remains finite. The reported maxima are numerical diagnostics at declared quadrature orders and parent dimensions; they are not certified error bounds. Agreement of the direct and projected routes does not establish convergence in a topology sufficient for the Weil criterion.

No finite sign, small residual, stable parent-dimension sequence, or agreement between implementations establishes density/completeness, complete-space Weil positivity, the semilocal Hilbert-space statement, Conjecture 4.1, or the Riemann hypothesis.

The useful next mathematical question after stable direct-pairing agreement is no longer implementation identity. It is to identify a topology and approximation theorem strong enough to control the Weil functional as the non-polynomial family is approximated by finite admissible subspaces.
