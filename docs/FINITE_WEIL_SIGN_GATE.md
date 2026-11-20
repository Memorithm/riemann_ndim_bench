# Conditional finite Weil sign gate

Status: implemented finite-dimensional diagnostic; **not an RH proof**.

Implementation: [src/weil_sign_gate.rs](../src/weil_sign_gate.rs).

## Mathematical contract

The existing compact Riemann--Weil matrix builder provides a real-symmetric
floating-point matrix B indexed by a finite family of compact,
boundary-admissible Weil test functions. Let A denote the corresponding exact
finite pairing matrix in the **same** basis.

The caller supplies a full n-by-n array of **justified absolute entrywise
error bounds** E, assumed to satisfy

    |A_ij - B_ij| <= E_ij.

These E values must account for the source formula, every quadrature stage,
all approximations, symmetrization, and floating-point roundoff. The numerical
difference between two quadrature refinements is **not** such a bound.

For a nonzero real coefficient vector v, the finite quadratic form satisfies

    |v^T A v - v^T B v|
      <= sum_(i,j) |v_i| max(E_ij,E_ji) |v_j|.

Dividing by v^T v produces the conditional Rayleigh interval reported by
the gate. When its upper endpoint is strictly negative, v is a negative
direction for A, provided all error assumptions hold.

For all real-symmetric matrices A satisfying the entrywise enclosure,
the Gershgorin theorem also yields

    lambda_min(A)
      >= min_i [
           B_ii - E_ii
           - sum_(j != i) (|B_ij| + max(E_ij,E_ji))
         ].

A **strictly positive** right-hand side implies positive definiteness of A
on the declared finite span. This is a sufficient but not necessary condition.
The gate intentionally reports Undetermined when neither strict margin wins.

## Scope / non-implications

- The matrix and the test-function family are finite-dimensional.
- The norm in the Rayleigh quotient is the Euclidean coefficient norm, not a
  source semilocal Hilbert norm or a Gram-normalized generalized spectrum.
- Positivity on one finite span does not establish positivity of the complete
  Weil form, convergence as n grows, semilocal sufficiency, or RH.
- A negative direction, **if genuinely certified** for the correct full Weil
  pairing, would be a serious mathematical obstruction requiring independent
  verification; the present floating-point result alone is not a refutation.
- This implementation uses binary64 arithmetic without outward rounding.
  The mathematical implications are conditional on error assumptions and on
  reliable arithmetic evaluation, not machine-verifiable certificates.
- Supplying E=0 for quadrature outputs is only meaningful when exactness has
  been independently demonstrated; no such assumption is made by this module.

## Reproducible tests

    cargo fmt --all -- --check
    cargo clippy --locked --all-targets --all-features -- -D warnings
    cargo test --locked weil_sign_gate
    cargo test --locked --all-targets

The included tests exercise positive, negative, and inconclusive synthetic
finite matrices, invalid inputs, and the adapter to the existing compact Weil
matrix builder. They do not substitute for validated error enclosures.

## Next proof-obligation milestone

Introduce genuine, independently validated error envelopes for the compact
pairing matrix, ideally interval/outward-rounded quadrature bounds with
controlled input roundoff. The current H1 continuity module can supply
conditional analytic estimates for projection errors; it does not, by itself,
supply rigorous entrywise integration errors.

The next Riemann-specific step remains the open semilocal-to-Weil bridge in
[RH_PROOF_DEPENDENCY_AUDIT.md](RH_PROOF_DEPENDENCY_AUDIT.md).

Relevant literature already used by the bench:

- A. Connes and C. Consani, *Weil positivity and Trace formula, the
  archimedean place*, arXiv:2006.13771.
- A. Connes, C. Consani, H. Moscovici, *Zeta zeros and prolate wave
  operators*, arXiv:2310.18423.
