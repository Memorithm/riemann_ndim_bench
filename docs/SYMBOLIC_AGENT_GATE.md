# Symbolic mu² agent gate

This layer connects the exact symbolic helper introduced by PR #12 to the
read-only research harness introduced by PR #10.

It deliberately keeps the base runners unchanged.  Two wrapper entrypoints
install an idempotent extension of the existing `verify_math` tool:

- `tools/riemann_symbolic_research_agent.py`
- `tools/riemann_dual_symbolic_research_agent.py`

The original entrypoints remain available for experiments that should stop at
the recurrence/Gamma/asymptotic layer.

## Added verifier modes

### `symbolic_hypergeometric`

The agent supplies:

- rational numerator Pochhammer shifts;
- rational denominator shifts;
- a rational base;
- an explicit candidate quotient `t[k+1]/t[k]`;
- an optional number of displayed initial terms.

The extension dispatches to `tools/symbolic_mu2.py hypergeometric` through a
fixed argv list.  The candidate quotient is mandatory at the agent boundary.
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

The ledger now understands the two symbolic modes directly.

A symbolic candidate mismatch is a hard `refuted` record.  The convenience
predicate

```python
ledger.has_exact_symbolic_mu2_chain()
```

is true only if both the hypergeometric quotient and finite-part stages have
exact successful records.

The dual symbolic wrapper adds both modes to the existing final required-mode
and exact-mode sets.  Consequently the final synthesis cannot pass merely
because recurrence transformation, perturbative extraction, Gamma arithmetic
and asymptotic fitting succeeded.

## Provenance boundary

The symbolic tools prove statements about the expressions supplied to them.
They do **not** prove that an agent obtained those expressions from the correct
semilocal recurrence.

The agent prompts therefore require a separate provenance chain:

1. exact raw-to-staggered recurrence normalization;
2. exact formal perturbative extraction through `mu^2`;
3. derivation of the parity/subsequence forcing quotient from that verified
   inhomogeneous equation;
4. exact `symbolic_hypergeometric` audit of the proposed quotient;
5. derivation of the weighted generating expression from the same sequence;
6. exact `symbolic_finite_part` audit;
7. Gamma and normalization assembly only after the preceding steps are
   justified.

A finite-part success on an arbitrary expression is therefore not enough to
establish the semilocal limiting constant.  The final report must state the
provenance argument explicitly and preserve any gap in it.

## Blindness boundary

No final limiting constant is encoded by this extension.  The regression suite
for `symbolic_mu2.py` contains exact known algebraic checkpoints for the
post-perturbative chain, but the agent-facing extension accepts generic safe
expressions and does not reveal a target constant.

## Scientific boundary

This is research infrastructure for exact-vs-asymptotic-vs-numerical auditing.
It does not identify finite semilocal spectra with Riemann-zeta zeros and makes
no claim of proving the Riemann hypothesis.
