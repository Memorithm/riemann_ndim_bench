# Conditional Weil block-Schur research gate

Status: finite-dimensional Rust code plus an abstract **conditional** operator
lemma. This work is neither a full Weil-positivity result nor a proof of RH.

Implementation: [src/weil_schur_gate.rs](../src/weil_schur_gate.rs).

## Precise mathematical statement

For a real-symmetric (or self-adjoint) block operator

    M = [ A   B  ]
        [ B*  D  ]

on H_head (+) H_tail, assume:

1. A >= alpha I, where alpha > 0;
2. D >= delta I;
3. B is bounded, with ||B||_2 <= beta.

Then A is invertible and

    D - B* A^-1 B >= (delta - beta^2 / alpha) I.

Consequently alpha > 0 and delta - beta^2 / alpha > 0 imply
strict positivity of M on its natural quadratic-form domain, under the
standard domain/boundedness assumptions. This is a sufficient condition,
not a necessary one. The code evaluates its *numerical margin*; it does not
prove that these operator assumptions hold.

For a **finite** Weil pairing matrix, partition the explicitly declared
basis into leading k coordinates and the remaining n-k coordinates.
The adapter computes conservative Gershgorin lower bounds for the two
diagonal blocks, using entrywise error budgets, and an upper bound for
the cross-block operator norm using its Frobenius norm. This can succeed
when plain Gershgorin on the full matrix is inconclusive.

## Entrywise error contract and rounding

Let B be the symmetrized floating-point matrix obtained from the existing
Weil quadrature implementation and A_exact the exact matrix on the same
finite basis. The caller MUST justify E_ij such that

    |(A_exact)_ij - B_ij| <= E_ij.

Both diagonal and off-diagonal E values must cover all errors, including
quadrature, source approximation, symmetrization, and binary64 arithmetic.

The finite adapter uses

    alpha_lower = min_{i<k} [ B_ii - E_ii
       - sum_{j<k,j!=i} (|B_ij| + max(E_ij,E_ji)) ];

    delta_lower = min_{i>=k} [ B_ii - E_ii
       - sum_{j>=k,j!=i} (|B_ij| + max(E_ij,E_ji)) ];

    beta_upper = sqrt(
       sum_{i<k,j>=k} (|B_ij| + max(E_ij,E_ji))^2
    ).

If alpha_lower is positive, the code evaluates the remaining margin

    delta_lower - beta_upper^2 / alpha_lower.

No conclusion is drawn when this is nonpositive or alpha_lower is
nonpositive. A negative or zero diagnostic margin alone does NOT establish
that the actual Weil form has a negative direction.

**Important:** the gate uses ordinary IEEE-754 binary64 operations, not
outward-rounded interval arithmetic, so all reported margins are conditional
floating-point diagnostics. Near-zero signs must not be treated as rigorous
certificates. A future certified implementation must enclose all operations,
not only the underlying matrix entries.

## Connection to the latest research

Ronnie Andrews, Alain Connes, Antonio La Gatta, *The uncertainty operator,
the Weil ground state and the zeros of zeta: a numerical study and conditional
positivity of the Weil quadratic form*, arXiv:2611.13129, 9 November 2026.

Primary source: https://arxiv.org/abs/2611.13129

The abstract states that, under an additional hypothesis valid at the
archimedean place, Weil positivity can be reduced to a Schur complement
involving a finite-dimensional near-self-dual space and a positive
complement. This is **not** the same as taking an arbitrary leading block of
our Legendre quadrature matrix. The paper's decomposition, its Hilbert
space, the tail coercivity constant and the cross-block norm must each be
derived or sourced separately before the abstract lemma can be applied.

In particular, our finite matrix's "tail block" consists only of the
remaining finite basis vectors. It is not the infinite-dimensional positive
complement in that paper.

## Operational tests

    cargo fmt --all -- --check
    cargo clippy --locked --all-targets --all-features -- -D warnings
    cargo test --locked weil_schur_gate
    cargo test --locked --all-targets
    cargo run --locked --release

Synthetic tests cover (i) a block decomposition for which global
Gershgorin fails but Schur passes, (ii) large uncertainty budgets,
(iii) positive decoupled blocks, (iv) invalid data and (v) integration
with the existing compact Weil matrix builder. They do not constitute a
certified bound for a real Weil pairing.

## Next research milestones

1. Derive actual, independently verified matrix-entry error bounds,
   preferably with outward-rounded/ball-arithmetic integration.
2. Match the paper's near-self-dual projection and the correct Sonin-tail
   operator to concrete source-locked inputs.
3. Prove (rather than assume) coercivity of the tail and a valid
   cross-block coupling norm, uniformly in the needed cutoff parameter.
4. Transfer a valid block bound into the compact-support semilocal
   positivity obligation recorded in
   [RH_PROOF_DEPENDENCY_AUDIT.md](RH_PROOF_DEPENDENCY_AUDIT.md).

Until 1-4 hold, do not assert a global Weil positivity statement or RH.
