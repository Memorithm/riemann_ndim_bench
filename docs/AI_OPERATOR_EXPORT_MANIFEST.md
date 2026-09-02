# RB5 — Structured-operator export manifest

Status: **research export contract — no AI utility or performance claim**

This document defines a conservative handoff surface from `riemann_ndim_bench`
to downstream Memorithm research benches. It does not change the mathematical
status of any RiemannBench result, does not claim an attention mechanism is
useful, and does not create a direct production-promotion path.

The governing rule is:

```text
RiemannBench mathematical hypothesis
        ↓
ADA semantic candidate + independent oracle/falsification
        ↓
ITD/TDI independent mechanistic experiments
        ↓
possible downstream promotion only after destination-owned gates
```

Riemann/zeta interpretation does not transfer through this interface.

## 1. Evidence classes

Every exported operator record MUST carry exactly one primary evidence class.
The class describes what RiemannBench currently supports; downstream use does
not strengthen it.

### `exact_identity`

A complete finite mathematical identity or derivation is available under stated
assumptions. Numerical agreement is supporting validation, not part of the
proof status.

### `source_locked_reproduction`

The construction reproduces a published formula or benchmark under recorded
conventions and tolerances. This does not make downstream extrapolations exact.

### `numerical_validation`

The construction has finite-dimensional numerical evidence with declared
method, precision/dimension scope and cross-checks. Repetition does not promote
it to a theorem.

### `formal_asymptotic`

A derivation contains an explicitly open uniformity, limit-interchange, domain
or related proof gap. The gap must travel with the export.

### `heuristic`

The object or relationship is exploratory. Fit quality, high dimension or
repeated agreement does not upgrade this status.

## 2. Required export record

A downstream-consumable operator hypothesis must identify:

```text
operator_id
operator_version
evidence_class
mathematical_definition
source_or_derivation
finite_domain_and_dimensions
parameter_domain
proved_properties
numerically_validated_properties
open_gaps
reference_fixtures
non_transferable_interpretations
preferred_downstream_route
```

The record must distinguish proved properties from numerically observed ones.
An empty `open_gaps` field is valid only when the exported statement itself has
no known mathematical gap; it does not assert usefulness for AI.

## 3. Initial conservative inventory

The entries below are **candidate export families**, not adopted AI mechanisms.
Before a downstream implementation depends on one, the exact formula and the
supporting RiemannBench document/fixture must be pinned in a versioned record.

### RB5-OP-TOEPLITZ — Toeplitz / triangular Toeplitz operators

Current basis in this repository:

- symmetric Toeplitz constructions and spectral solvers are part of the
  numerical foundation;
- triangular operator structure is also present in the explicit Hardy/Copson
  factorisation work.

Permitted downstream hypothesis:

- use a precisely defined finite Toeplitz or triangular Toeplitz operator as a
  structured sequence/operator candidate.

Not transferred:

- any Riemann/zeta interpretation;
- any claim of attention quality, long-context superiority or computational
  advantage.

### RB5-OP-PROLATE — finite prolate / concentration operators

Current basis in this repository:

- source-locked finite prolate constructions are part of the Phase 3/4 lineage;
- finite spectral computations and reference checks exist for the declared
  RiemannBench problems.

Permitted downstream hypothesis:

- test a frozen finite prolate/concentration operator or parameterisation as an
  independently evaluated structured mixer or projection.

Not transferred:

- identification of finite crossings with zeta zeros;
- evidence that spectral concentration improves an AI task;
- low-rank, memory or runtime claims without separate measurement.

### RB5-OP-SPECTRAL-FLOW — controlled finite operator perturbation

Current basis in this repository:

- Phase 4 studies a source-derived first-order perturbation of a finite operator;
- eigenvalue/eigenvector and sign behaviour are analysed within the documented
  finite-dimensional scope.

Permitted downstream hypothesis:

- use an explicitly defined finite one-parameter operator deformation as a
  candidate mechanism for controlled context-dependent operator change.

Not transferred:

- semilocal-prime meaning;
- causal or adaptive-compute usefulness;
- any claim that the RiemannBench perturbation is the correct AI parameterisation.

### RB5-OP-GROUND-STATE — weighted ground-state factorisation

Current basis in this repository:

- an exact ground-state factorisation is documented in the Phase 4 lineage;
- the associated free-Laplacian/soft-edge programme has separate numerical and
  formal components whose evidence classes must remain distinct.

Permitted downstream hypothesis:

- export the finite factorised operator under its exact stated assumptions as a
  structured positive/difference-operator candidate where those properties are
  actually proved.

Not transferred:

- unresolved asymptotic conclusions;
- attention stability, training stability or quality claims.

### RB5-OP-GREEN — finite inverse / Green-kernel construction

Current basis in this repository:

- the Phase 4 lineage records an explicit finite inverse-kernel construction
  associated with the ground-state factorisation.

Permitted downstream hypothesis:

- evaluate a pinned finite kernel as a long-range structured interaction
  candidate after ADA or another downstream bench builds its own oracle and
  adversarial tests.

Not transferred:

- a claim that the kernel is a universal long-range memory mechanism;
- runtime or complexity claims;
- Riemann relevance to the downstream task.

### RB5-OP-HARDY-COPSON — triangular factor family

Current basis in this repository:

- the active uniform-trace programme uses an explicit triangular Hardy/Copson
  factor in the factorisation of the finite inverse.

Permitted downstream hypothesis:

- export a finite triangular factor with its exact indexing, coefficients and
  boundary conventions for independent study.

Not transferred:

- the unresolved uniform trace-control objective;
- any downstream conditioning or optimization advantage.

### RB5-OP-SOFT-EDGE — low-mode / soft-edge finite structures

Current basis in this repository:

- RiemannBench contains finite numerical and formal work on free-Laplacian
  low-mode/soft-edge behaviour and fixed-rank corrections.

Export restriction:

- only a finite mathematical construction with an independent definition may be
  exported;
- asymptotic formulas must remain labelled `formal_asymptotic` or `heuristic`
  whenever the corresponding uniform proof gap remains open.

Not transferred:

- an asymptotic theorem stronger than the repository has proved;
- a claim that low-mode structure improves sequence modelling.

## 4. Downstream responsibilities

### ADA

ADA may convert one pinned operator record into a bounded semantic candidate.
It must start from zero for attention/sequence-specific correctness, numerical,
adversarial, task-quality and prior-art evidence. RiemannBench evidence is not
an ADA qualification result.

### TDI

TDI may use a pinned operator as a newly defined reference mechanism in a new
preregistered experiment. RiemannBench does not provide TDI dynamic-information
semantics or intervention/recovery evidence.

### ITD Simulator

ITD may compare structural behaviour of a pinned operator in its isolated AI
research namespace. Frozen ITD V29.18 is not modified by this export.

### SciRust

A primitive that proves genuinely general and reusable should be considered for
promotion to SciRust with independent tests and review. Riemann-specific
interpretation stays here.

### FLAT-ATTENTION / NNIS / ElasticXxx

There is no direct promotion path from this manifest. Production execution,
hardware performance and runtime adaptation require their own destination-owned
qualification.

## 5. Promotion checklist for one concrete record

Before an operator moves from this inventory to a versioned downstream record:

1. pin the exact defining document/formula and repository commit;
2. state the finite domain, dimensions and parameter bounds;
3. label the weakest applicable evidence class;
4. list proved and numerically validated properties separately;
5. include deterministic reference fixtures when available;
6. list all known mathematical gaps relevant to the exported statement;
7. state explicitly which Riemann/zeta interpretations do not transfer;
8. select one bounded downstream research question;
9. require the downstream repository to build an independent oracle/evidence
   ladder;
10. preserve a negative downstream result as a valid outcome.

This document is an interface specification only. It does not itself export or
qualify a new AI mechanism.
