# Mixed localized Bernstein finite Weil subspace

This experiment strengthens the one-direction Bernstein probes by forming a
full generalized finite Weil problem on several selected localized Bernstein
functions simultaneously.

For degree `D`, let `c_k` be the shifted-Legendre coefficient vector of

\[
B_{D,k}(t)=\binom{D}{k}t^k(1-t)^{D-k}.
\]

For selected indices `k_1,...,k_N`, assemble the coefficient matrix

\[
C=[c_{k_1}\;\cdots\;c_{k_N}].
\]

The existing validated parent matrices are restricted by congruence:

\[
A_B=C^TAC,\qquad G_B=C^TGC,
\]

and the audit solves

\[
A_Bv=\lambda G_Bv.
\]

No second Riemann--Weil decomposition is introduced.

## Why this is stronger than a basis rescaling

A complete degree-`D` Bernstein family spans the same polynomial space as
Legendre degrees `0..D`, so changing the complete basis alone cannot change the
generalized spectrum. This fact is used as a regression: the complete degree-2
Bernstein basis must reproduce the generalized spectrum of the three-dimensional
Legendre parent space.

The experiment of interest selects fewer than `D+1` Bernstein functions. For
example,

\[
D=6,\qquad k\in\{0,2,4,6\}
\]

defines a four-dimensional localized subspace inside the seven-dimensional
Legendre parent space. It is therefore generally different from the
dimension-matched leading Legendre subspace

\[
\operatorname{span}\{h_0,h_1,h_2,h_3\}.
\]

Both are evaluated from the same parent `(A,G)` computation.

## Recorded evidence

The audit records:

- Bernstein degree and selected indices;
- parent Legendre dimension;
- raw finite Weil eigenvalues;
- multiplicative Gram eigenvalues;
- generalized eigenvalues;
- Gram condition number;
- reconstructed critical boundary residual;
- whitening asymmetry;
- parent raw directional-pairing asymmetry;
- dimension-matched leading-Legendre generalized minimum and Gram condition.

The leading-Legendre comparison is a control, not a claim that either finite
subspace is canonical.

## Scientific boundary

A positive minimum generalized eigenvalue on this selected Bernstein subspace
would establish only positivity of the numerically represented restriction,
subject to the reported numerical diagnostics. It would not prove:

- positivity on all degree-`D` directions;
- positivity on the complete admissible Weil space;
- density or completeness of the tested family;
- uniform convergence of finite restrictions;
- identification of the numerical Gram norm with semilocal `L^2(X_S)`;
- Conjecture 4.1;
- the Riemann hypothesis.

Likewise, a negative finite result would require independent replication,
quadrature refinement, conditioning analysis, and preferably a second
implementation before mathematical interpretation.

## Reproduction

```bash
cargo run --release --example weil_bernstein_subspace
```

The default example compares degree-6 Bernstein indices `0:2:4:6` with the
leading four-dimensional Legendre control, at the same compact support and
quadrature settings.

## Next research step

If the mixed localized subspace audit is stable, cross it with the existing
quadrature-refinement and support-window axes. That will distinguish a real
subspace effect from conditioning or resolution drift. Only after that should
we introduce a compact generator family not contained in the same finite
polynomial parent space.
