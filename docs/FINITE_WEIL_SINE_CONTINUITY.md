# Finite sine Weil continuity probe

This audit connects two quantities that were previously reported separately:

1. the approximation error between the direct compact sine function
   `h = Q(bump * sin(m*pi*t))` and its finite shifted-Legendre approximation
   `h_N`;
2. the residual between their complete finite compact Riemann--Weil pairings.

For each mode the approximation error is measured numerically in

`||f||_2^2 = integral |f(rho)|^2 d rho / rho`,

using declared Gauss--Legendre quadrature in logarithmic coordinate. The audit
also records the maximum pointwise error observed at those quadrature nodes.
Neither quantity is a certified uniform error bound.

For a pair of modes `(i,j)` the diagnostic perturbation scale is

`e_i max(||h_j||_2, ||h_j,N||_2) + e_j max(||h_i||_2, ||h_i,N||_2)`,

where `e_k = ||h_k-h_k,N||_2`.

The reported `observed_l2_continuity_quotient` is the absolute direct/projected
Weil pairing residual divided by this scale when the scale is nonzero. It is a
finite numerical response quotient. It is not asserted to be an operator norm
or a bound valid outside the sampled family.

The pairing decomposition remains the same source-locked implementation used by
`weil_compact_pairing`: critical pole term minus the real-place term minus the
finite prime-power total. The probe reports residuals of those three pieces
separately so that apparent stability cannot hide cancellation between source
terms.

## Interpretation boundary

A stable or bounded observed quotient across the tested parent dimensions would
be useful evidence that the chosen finite approximants behave continuously in
this declared numerical `L2(d*rho)` metric. It would not prove that the Weil
functional is continuous in this norm on a complete admissible function space.
In particular this module does not establish:

- a certified quadrature error bound;
- uniform convergence of the Legendre series or its Q-images;
- continuity of the complete Weil functional in an analytically specified
  topology;
- density or completeness of the sine or polynomial test families;
- complete-space Weil positivity;
- the semilocal `L2(X_S)` identification;
- Conjecture 4.1;
- the Riemann hypothesis.

The intended next mathematical step, if the numerical quotient remains stable,
is to identify a topology strong enough to control the archimedean singular
kernel and prime-power evaluations and then prove an explicit continuity
estimate there. The present probe supplies diagnostics for choosing and testing
that topology; it is not that proof.
