# Finite Weil subspace sensitivity audit

This layer probes a different numerical axis from support width, finite dimension,
and quadrature refinement: the choice of the finite test-function subspace.

The parent family remains

\[
g_j(\rho)=\operatorname{bump}(\rho)P_j(2t-1),\qquad h_j=Qg_j,
\]

with the same compact support, boundary operator, Riemann--Weil decomposition,
and multiplicative Gram normalization used by the preceding finite Weil audits.

Instead of always taking the leading span

\[
\operatorname{span}\{h_0,\ldots,h_{N-1}\},
\]

this module accepts explicit strictly increasing degree sets such as

\[
\{0,1,2\}\quad\text{and}\quad\{0,2,4\}.
\]

The second set is not a diagonal rescaling of the first and is not the same
finite subspace.  It therefore provides a controlled sensitivity test while
holding the support, quadrature levels, and finite Weil functional fixed.

## Reuse of the validated parent matrices

For all requested subspaces the implementation computes one parent pairing
matrix `A` and one parent Gram matrix `G`, extending only to the largest requested
Legendre degree.  Each subspace is extracted as the corresponding principal or
non-principal indexed restriction of these same matrices.

This matters because a comparison between subspaces is then not confounded by a
second implementation of the Riemann--Weil decomposition.

For every selected subspace the audit records:

- the exact Legendre degree indices used;
- the raw finite Weil eigenvalues;
- the multiplicative-Gram eigenvalues;
- the generalized eigenvalues of `A v = lambda G v`;
- the Gram condition number;
- the maximum selected boundary residual;
- the whitening asymmetry residual.

It also carries the maximum raw directional pairing asymmetry of the full parent
matrix as a conservative shared diagnostic.

## Numerical interpretation

The generalized spectrum removes arbitrary diagonal basis scaling inside a
fixed selected span, but it does not make spectra from different subspaces
identical or directly comparable to a complete-space operator spectrum.
Differences between selected subspaces are empirical finite-dimensional
sensitivity diagnostics.

No automatic significance or proof threshold is defined.

## Scientific boundary

Agreement of positive generalized minima across several selected finite
subspaces would strengthen finite numerical evidence only.  It would not prove:

- positivity of the Weil quadratic form on the complete admissible space;
- density or completeness of the tested subspaces;
- uniform convergence of the finite restrictions;
- identification of the numerical Gram norm with semilocal `L^2(X_S)`;
- Conjecture 4.1;
- the Riemann hypothesis.

Likewise, a negative selected finite eigenvalue would require independent
replication, refinement, and conditioning checks before any mathematical
interpretation.

## Reproduction

```bash
cargo run --release --example weil_subspace_audit
```

The example compares the leading degree set `0:1:2:3` with the even degree set
`0:2:4:6` at the same support and quadrature levels.

## Next research step

If this subspace audit is numerically stable, the next stronger test is to add a
genuinely different compact generator family rather than another selection from
the Legendre parent family.  That step should preserve the same provenance,
boundary, Gram-conditioning, and refinement diagnostics and should not be mixed
with a simultaneous change of support or source normalization.
