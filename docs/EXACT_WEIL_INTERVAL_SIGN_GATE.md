# Exact conditional sign checks for finite Weil matrices

Status: **EXACT integer algebra**, conditional on the correctness of the input
rational enclosures. It is neither a certified integration result nor a
complete-space Weil positivity or Riemann-hypothesis proof.

Implementation: src/weil_rational_gate.rs.

## Motivation

The prior finite sign and block-Schur modules evaluate their margins in
binary64 floating point. This implementation decides signs with checked
128-bit integer operations only, eliminating roundoff from the final
decision. The unresolved upstream problem is still to generate *proven*
enclosures for the actual Weil matrix entries.

## Input contract

Let A be the unknown exact, real-symmetric finite Weil pairing. Provide
an integer S > 0 and a closed integer interval [L_ij,U_ij] for each entry
with the independently justified enclosure

    L_ij / S <= A_ij <= U_ij / S.

All intervals must be symmetric. Supply a nonzero integer coefficient
direction v and, optionally, an integer split 0 < k < n. The function
does not accept f64 values for these input enclosures, and must not be
fed optimistically rounded quadrature estimates without rigorous error
budgets.

## Exact deductions

Set m_ij = max(|L_ij|,|U_ij|). The full Gershgorin bound is

    g = min_i (L_ii - sum_(j!=i) m_ij).

If g>0 then A is positive definite on the finite span.

For a nonzero integer v, define the upper bound for S*(v^T A v):

    r = sum_(i,j) v_i*v_j * (U_ij when v_i*v_j>=0, else L_ij).

If r<0 then v is an exact negative witness on this finite span,
conditional on the declared enclosure.

For the split at k, calculate the Gershgorin block lower numerators
alpha and delta and the upper Frobenius-coupling squared numerator
beta2 = sum_(i<k,j>=k) m_ij^2. The strict conditions

    alpha>0 and alpha*delta-beta2>0

suffice for finite block positivity by the Schur complement. All
numerators have exact common denominators; the deciding comparison
requires neither square root nor floating-point division.

The result is Undetermined if neither implication is available. Checked
i128 operations reject overflows, including the magnitude of i128::MIN
when required. Larger rational certificates will need arbitrary
precision arithmetic.

## Soundness regression

The Rust tests exercise exact Gershgorin positivity, a 3x3 matrix
where full Gershgorin fails but Schur passes, increased uncertainty,
a negative Rayleigh direction, malformed inputs and overflow guards.
An exhaustive 2x2 test uses centers in [-2,2] and interval radii 0 or
1, checking every positive verdict against exact Sylvester minors and
every negative witness against its exact quadratic form.

    cargo fmt --all -- --check
    cargo clippy --locked --all-targets --all-features -- -D warnings
    cargo test --locked weil_rational_gate
    cargo test --locked --all-targets

## Still open

1. Construct genuinely validated, source-locked entry enclosures for
   the actual compact Weil pairing, including integration, approximation
   and rounding errors.
2. Run this exact checker on those enclosures and retain provenance.
3. Identify the correct near-self-dual projector and infinite-dimensional
   Sonin complement; prove uniform coercivity and coupling estimates.
   The finite Legendre split used here is not the infinite complement.
4. Close the global Weil-positivity implication in the dependency audit.

See docs/RH_PROOF_DEPENDENCY_AUDIT.md and
docs/CONDITIONAL_WEIL_SCHUR_GATE.md. The operator literature and
assumptions are cited there, particularly Andrews, Connes and La Gatta,
arXiv:2611.13129 (9 November 2026).
