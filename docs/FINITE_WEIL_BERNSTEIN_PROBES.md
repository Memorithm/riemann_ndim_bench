# Localized finite Weil Bernstein probes

This layer introduces a different family of compact polynomial enrichments while
reusing the validated finite Weil pairing matrix `A` and multiplicative Gram
matrix `G`.

For a fixed degree `D`, define

\[
B_{D,k}(t)=\binom{D}{k}t^k(1-t)^{D-k},\qquad 0\le k\le D,
\]

and

\[
g_{D,k}(\rho)=\operatorname{bump}(\rho)B_{D,k}(t),\qquad
h_{D,k}=Qg_{D,k}.
\]

Bernstein functions are spatially localized in the affine support coordinate
`t`, unlike individual global Legendre modes. They therefore provide a useful
new set of finite directions for sensitivity testing.

## Reuse of the parent Weil computation

Every degree-`D` Bernstein polynomial is expanded as

\[
B_{D,k}(t)=\sum_{j=0}^{D} c_{j,k}P_j(2t-1).
\]

Since `Q` is linear,

\[
h_{D,k}=\sum_{j=0}^{D} c_{j,k}h_j.
\]

The implementation therefore computes the existing Legendre parent matrices
only once at dimension `D+1`, then evaluates

\[
A_k=c_k^T A c_k,
\qquad
G_k=c_k^T G c_k,
\qquad
R_k=\frac{A_k}{G_k}.
\]

No second implementation of the Riemann--Weil prime/archimedean decomposition is
introduced.

## Coefficient construction

The shifted Legendre identity

\[
P_j(2t-1)=(-1)^j\sum_{m=0}^{j}(-1)^m
\binom{j}{m}\binom{j+m}{m}t^m
\]

is combined with the beta moment

\[
\int_0^1 B_{D,k}(t)t^m\,dt
=
\frac{D!(k+m)!}{k!(D+m+1)!}.
\]

The resulting `f64` coefficients are audited by reconstructing each Bernstein
polynomial on a deterministic grid. The reconstruction residual is a numerical
implementation check, not a proof error bound.

## Recorded evidence

For every selected Bernstein index the audit records:

- Legendre expansion coefficients;
- coefficient L1 norm;
- raw finite Weil quadratic value `c^T A c`;
- multiplicative Gram norm squared `c^T G c`;
- generalized one-dimensional Rayleigh quotient;
- both critical boundary moments and their maximum residual;
- maximum sampled Bernstein reconstruction residual.

The shared parent audit additionally reports Gram conditioning, raw pairing
asymmetry, and whitening asymmetry.

## Scientific boundary

The full degree-`D` Bernstein family spans exactly the same polynomial space as
Legendre degrees `0..D`. Consequently, this experiment probes different
localized **directions** inside that finite parent space. It does not by itself
supply a new dense test-function class.

Positive Rayleigh quotients for many Bernstein directions do not prove:

- positivity of the finite parent matrix in untested directions;
- positivity on the complete admissible Weil space;
- density or completeness;
- a uniform approximation theorem;
- semilocal `L^2(X_S)` identification;
- Conjecture 4.1;
- the Riemann hypothesis.

A negative direction would likewise require independent replication,
quadrature refinement, and conditioning analysis before interpretation.

## Reproduction

```bash
cargo run --release --example weil_bernstein_probe
```

The example probes indices `0,2,4,6` in the degree-6 Bernstein family on the
same manufactured support used by the existing finite Weil audits.

## Next step

The next stronger experiment is to combine selected localized Bernstein
functions into a mixed finite subspace, form `C^T A C` and `C^T G C`, and audit
its full generalized spectrum. Selecting fewer than all `D+1` Bernstein
functions then produces a genuinely different finite subspace from the leading
Legendre span of the same dimension.
