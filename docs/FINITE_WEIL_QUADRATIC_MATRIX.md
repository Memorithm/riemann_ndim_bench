# Finite Weil quadratic-form matrix

Status: **FINITE-BASIS NUMERICAL RESTRICTION / NOT WEIL POSITIVITY**

This increment lifts the scalar convolution-square audit to a finite-dimensional real-symmetric matrix of mixed Riemann--Weil pairings.

## Basis

Fix the same compact rational support `(a,b)` used by the compact boundary audit. Let

```text
t = (rho-a)/(b-a),
y = 2t-1,
g_j(rho) = bump(rho) P_j(y),
h_j = Q g_j,
Q = -(rho d/drho)^2 + 1/4,
```

where `P_j` is the Legendre polynomial of degree `j`.

The bump is `C^infinity` with compact support, and multiplying it by a polynomial does not enlarge support. The implementation differentiates the Legendre recurrence analytically and evaluates `Q g_j` without finite differences.

Each `h_j` is independently audited against the two critical Mellin boundary moments.

## Matrix entries

For basis vectors `h_i,h_j`, the raw mixed entry is

```text
A_ij = psi(h_i^* * h_j).
```

The mixed log correlation is evaluated directly. Both orientations `(i,j)` and `(j,i)` are computed for off-diagonal entries; their difference is reported as a numerical symmetry residual. The stored matrix entry is the average of the two independently evaluated orientations.

The source decomposition remains separated into:

- the residual critical pole term;
- the real-place distribution `W_R`;
- the finite prime-power sum;
- the resulting `psi` pairing.

The support ratio is derived from the exact rational endpoints, so the non-archimedean prime-power window has the same exact endpoint treatment as the scalar audit.

## Spectral audit

The symmetric matrix is diagonalized with the repository's existing `faer::linalg::solvers::SelfAdjointEigen` machinery. The eigenvalues are sorted, and the smallest eigenvalue is exposed.

Tests do **not** assume that this eigenvalue is positive. Instead they check:

1. degree-zero agreement with the scalar `h=Qg` implementation;
2. one-dimensional agreement with the scalar Weil functional audit;
3. independently evaluated off-diagonal symmetry;
4. critical-boundary residuals for every basis vector;
5. Cauchy interlacing consistency of the minimum eigenvalue for leading principal basis expansions;
6. convergence of the minimum eigenvalue under quadrature refinement.

The ordinary eigenvalues of `A` use the Euclidean norm of the coefficient vector and therefore depend on the normalization of the chosen basis. Their magnitude is not an intrinsic spectrum of the admissible function space. The sign/inertia of the quadratic form on the declared finite span is the basis-invariant information relevant at this stage.

A subsequent normalization layer should build an independent positive Gram matrix, for example in the declared multiplicative `L^2(d*rho)` norm, and solve the generalized problem

```text
A v = lambda G v.
```

Only that normalized generalized spectrum should be used to compare eigenvalue magnitudes across differently scaled bases or growing dimensions.

## Scientific boundary

If every eigenvalue of a declared finite matrix is non-negative, this means only that the Weil quadratic form is numerically non-negative on that finite span at the declared resolution. It does not establish positivity on the full infinite-dimensional admissible test-function space.

Conversely, a stable negative eigenvalue would be scientifically important for this chosen finite span because it would exhibit a numerical negative direction for the implemented source functional; it would require independent verification before interpretation.

Neither sign proves or disproves RH by itself. Extending a finite-basis observation to the Weil criterion would require a justified density/completeness argument and uniform control of numerical/operator errors. Conjecture 4.1 is not assumed by this module.
