# Symbolic mu² agent gate

This layer connects the exact symbolic helper introduced by PR #12 to the
read-only research harness introduced by PR #10.

It deliberately keeps the base runners unchanged. Two wrapper entrypoints
install an idempotent extension of the existing `verify_math` tool:

- `tools/riemann_symbolic_research_agent.py`
- `tools/riemann_dual_symbolic_research_agent.py`

The original entrypoints remain available for experiments that should stop at
the recurrence/Gamma/asymptotic layer.

## Added verifier modes

### `symbolic_forcing_ratio`

This mode closes the provenance gap between the exact second-order recurrence
and a proposed hypergeometric forcing family.

Starting from

```text
A_n * Delta[n+1] - B_n * Delta[n] = F_n
```

let `h` be a homogeneous increment satisfying

```text
h[n+1] / h[n] = B_n / A_n.
```

For the variation-of-constants forcing

```text
T_n = F_n / (A_n * h[n+1]),
```

the unknown normalization of `h` cancels and one obtains the exact quotient

```text
T[n+1] / T[n] = F[n+1] * A[n] / (F[n] * B[n+1]).
```

The agent supplies:

- `A(j)`;
- `B(j)`;
- the exact inhomogeneous forcing `F(j,n)`;
- the rational lattice offset in `j = n + offset`;
- its proposed rational quotient for `T[n+1]/T[n]`.

The helper performs the affine substitution and rational simplification exactly.
A successful result contains

```text
candidate_status=PROVED_EQUAL
exact_status=PROVED_BY_VARIATION_OF_CONSTANTS_QUOTIENT
```

A wrong proposal contains

```text
candidate_status=MISMATCH
exact_status=REFUTED_FORCING_QUOTIENT
```

This is the required bridge from the verified perturbative recurrence to the
next hypergeometric stage.

### `symbolic_hypergeometric`

The agent supplies:

- rational numerator Pochhammer shifts;
- rational denominator shifts;
- a rational base;
- an explicit candidate quotient `t[k+1]/t[k]`;
- an optional number of displayed initial terms.

The extension dispatches to `tools/symbolic_mu2.py hypergeometric` through a
fixed argv list. The candidate quotient is mandatory at the agent boundary.
This prevents the final gate from being satisfied merely by asking the helper
to construct an arbitrary hypergeometric sequence.

A successful result contains

```text
candidate_status=PROVED_EQUAL
exact_status=PROVED_BY_POCHHAMMER_QUOTIENT
```

A wrong candidate contains

```text
candidate_status=MISMATCH
exact_status=REFUTED_CANDIDATE_RATIO
```

and is classified by `ProofLedger` as `refuted`.

### `symbolic_finite_part`

The agent supplies:

- a base generating expression in `z`;
- a polynomial in the theta index `k`, representing `theta = z d/dz`;
- an optional extra expression;
- a bounded Puiseux expansion order.

The helper performs exact SymPy algebra and substitutes

```text
x = sqrt(1-z),  z = 1-x^2
```

before extracting the coefficient of `x^0`.

A successful result contains

```text
exact_status=PROVED_BY_EXACT_THETA_ALGEBRA_AND_PUISEUX_SERIES
```

## ProofLedger behavior

The ledger understands all three post-perturbative symbolic modes directly.

A candidate mismatch is a hard `refuted` record. The convenience predicate

```python
ledger.has_exact_symbolic_mu2_chain()
```

is true only if all three stages have exact successful records:

1. recurrence-derived forcing quotient;
2. hypergeometric/Pochhammer quotient;
3. local finite-part extraction.

The dual symbolic wrapper adds all three modes to the existing final
required-mode and exact-mode sets. Consequently the final synthesis cannot pass
merely because recurrence transformation, perturbative extraction, Gamma
arithmetic and asymptotic fitting succeeded.

## Provenance chain

The final report is expected to maintain the following chain explicitly:

1. exact raw-to-staggered recurrence normalization;
2. exact formal perturbative extraction through `mu^2`;
3. exact `symbolic_forcing_ratio` derivation for the relevant parity/site
   offset;
4. exact `symbolic_hypergeometric` audit of the coefficient family identified
   from that recurrence-derived quotient;
5. derivation of the weighted generating expression from the verified sequence;
6. exact `symbolic_finite_part` audit;
7. Gamma and normalization assembly only after the preceding steps are
   justified.

The new forcing-ratio stage removes a previous weak point: the agent may no
longer jump directly from an inhomogeneous recurrence to a visually familiar
binomial family without proving the quotient implied by variation of constants.

A finite-part success on an arbitrary expression is still not enough to
establish the semilocal limiting constant. The report must justify the weighted
generating expression itself and preserve any remaining gap.

## Blindness boundary

No final limiting constant is encoded by this extension. The regression suite
for `symbolic_mu2.py` contains exact algebraic checkpoints for the
post-perturbative chain, but the agent-facing extension accepts generic safe
expressions and does not reveal a target constant.

## Scientific boundary

This is research infrastructure for exact-vs-asymptotic-vs-numerical auditing.
It does not identify finite semilocal spectra with Riemann-zeta zeros and makes
no claim of proving the Riemann hypothesis.
