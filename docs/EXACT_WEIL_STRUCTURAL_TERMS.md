# Exact structural terms in the compact finite Weil pairing

**Evidence class: EXACT mathematical identity for the ideal smooth test
family**. Remaining archimedean and in-support prime-power integrals are
**NUMERICAL**, not certified. This is not a general Weil-positivity theorem
or a proof of RH.

Code: [src/weil_quadratic_matrix.rs](../src/weil_quadratic_matrix.rs).

## 1. Source-locked pole cancellation

Connes and Consani, *Weil positivity and Trace formula, the archimedean
place*, 25 June 2020, Lemma 3.3(ii), p. 18 of the arXiv PDF:
https://arxiv.org/pdf/2006.13771#page=18

Define the multiplicative logarithmic coordinate `x = log(rho)`.
For every ideal smooth compact generator

    g_j(rho) = bump(rho) * P_j(2*(rho-a)/(b-a)-1)

the present basis is

    h_j = Q g_j,     Q = -(rho d/drho)^2 + 1/4.

In `x`, the operator is `Q = -d_x^2+1/4`. Two integrations by parts
have no boundary contributions because `g_j(e^x)` and its derivatives
are compactly supported. Hence for every complex `s`:

    integral h_j(e^x) e^(s*x) dx
      = (1/4 - s^2) * integral g_j(e^x) e^(s*x) dx.

At `s = +/-1/2`, both critical moments vanish **exactly**. The mixed
pole term of the functional is a sum of products of critical moments:

    pole(i,j)
       = M_(h_i)(-1/2) * M_(h_j)(+1/2)
       + M_(h_i)(+1/2) * M_(h_j)(-1/2)
       = 0.

Before this change, the finite pairing accumulated a tiny but
conceptually unnecessary *numerically integrated* pole residual. Now
the source-level exact pole is set to 0 directly in each mixed entry.

The independent `boundary_residuals` diagnostic still integrates
each test function at the chosen diagnostic quadrature order.
Those values measure numerical consistency and **are not interpreted
as the mathematical Mellin moments of the ideal function**.

In particular, changing the `boundary_order` alone can now change
the reported numerical diagnostics but cannot change matrix entries.

### Scalar Weil functional

The existing `audit_finite_weil_functional` evaluates the scalar test
`h_1=Qg` from `CompactWeilTestFunction`. Its exact critical pole term
is now also zero. The scalar source value, like the Legendre matrix,
cannot change merely by refining the *diagnostic* Mellin quadrature
order; the `boundary()` audit itself remains available and numerical.

### Generic compact and projected-sine routes

The generic `CompactWeilEvaluand` trait now has the source-provenance
method `exact_q_image_mellin_zero()`, defaulting to **false** for
unqualified implementations. The known `CompactWeilBasisFunction` and
`CompactSineWeilFunction` implementors and the finite linear combinations
`ProjectedLegendreWeilFunction` and `ProjectedSineWeilFunction` opt
in because they all represent ideal `Qg` test functions.

In the generic mixed pairing, a guaranteed Q-image on **either** side
annihilates both critical pole summands, irrespective of the other factor.
The numerical moments are still calculated for diagnostic residuals.
If *neither* side has Q-image provenance, the generic route retains its
legacy numerical pole calculation; a small observed residual never upgrades
an arbitrary function into an exact Q-image. This preserves the generic
API's ability to inspect unqualified user-defined evaluands.

The direct sine and projected Legendre comparisons now share the same
source-level pole convention as the parent Legendre matrix; the
cross-route regression tests remain strict rather than being loosened.

## 2. Exact finite prime-power support

If `h_i,h_j` are both supported in `[a,b]` with `0<a<b`, then

    (h_i^* * h_j)(rho) = 0 outside [a/b, b/a].

For smooth compact bumps, the correlation also vanishes at its support
endpoints. Consequently, only integers `n >= 2` that are prime
powers and satisfy the *strict* inequality

    n < b/a

can contribute to the finite nonarchimedean explicit-formula sum.

The support endpoints in `CompactArchimedeanBump` are positive
rationals. The code computes the largest strict integer
`N = ceil(b/a)-1` **exactly** as

    N = floor(b/a)-1    if b/a is an integer,
        floor(b/a)      otherwise.

No floating-point log or endpoint comparison is used to decide
which prime powers might be active. Unlike the previous multiplication
`n * denominator`, the decision uses a ratio remainder and cannot
overflow from that multiplication.

In the particularly useful support window `b/a <= 2`, one has
`N <= 1`, so **all prime-power contributions are identically zero**.
The matrix audit exposes `max_strict_prime_argument()` and
`prime_power_sum_provably_zero()` as source-level support facts.

These statements concern only exact support and pole **terms**. The
archimedean integral still has to be evaluated and enclosed.

## 3. Numerical versus exact representations

The Rust function `CompactWeilBasisFunction::value_inside` uses a
binary64 exponential and a `LOG_UNDERFLOW_GUARD` for numerical
stability. It approximates the exact ideal `C_c^infinity` profile,
and is not literally the exact function at every point.

The preceding identities apply to the **ideal mathematical profile**
that the numerical implementation is intended to approximate.
Substituting an exact pole term does not make the other numerical
integrals rigorously enclosed. In particular, small diagnostics or
convergence of Gauss-Legendre quadrature do not justify setting
entrywise error budgets to zero.

## 4. Verification

    cargo fmt --all -- --check
    cargo clippy --locked --all-targets --all-features -- -D warnings
    cargo test --locked --all-targets
    cargo run --locked --release

The dedicated Rust tests verify that all diagonal and mixed pole terms
are exactly zero, that matrix coefficients are invariant when only
the *diagnostic* boundary quadrature order changes, and that rational
support windows with ratios 2, 5/2, 3 and 7/2 exclude or include the
correct **strict** integer cutoff. The ratio-2 test also checks that
the actual calculated nonarchimedean total is exactly zero.

The actual remaining quadrature and floating-point error budgets
are still open. Their independent certification is the next
research obligation. Even a certified finite sign cannot alone imply
RH; see [RH_PROOF_DEPENDENCY_AUDIT.md](RH_PROOF_DEPENDENCY_AUDIT.md).
