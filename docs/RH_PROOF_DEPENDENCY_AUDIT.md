# RH proof-dependency audit

Status: research-scope and proof-dependency audit.

This document recenters `riemann_ndim_bench` on the Riemann hypothesis (RH). It is intentionally stricter than a numerical roadmap: every arrow is classified by what is actually established, source-backed, conjectural, or absent.

## Status vocabulary

- `EXACT` — finite algebraic identity proved directly in the bench or reproduced from a source without an asymptotic limit.
- `SOURCE THEOREM` — theorem from the cited mathematical literature; the bench does not claim authorship.
- `SOURCE CONJECTURE` — conjectural implication stated in the cited literature.
- `CONDITIONAL` — implication valid only if explicitly stated hypotheses are established.
- `FORMAL ASYMPTOTIC` — derived asymptotic without the uniform remainder needed for a theorem.
- `NUMERICAL EVIDENCE` — finite numerical observation or regression.
- `OPEN BRIDGE` — logical implication currently missing from the project.
- `REFUTED` — a tested candidate implication or specificity claim that failed its controls.

## 1. The proof target

The final target is not a finite prolate spectrum, a fitted constant, or a numerical match to zeta zeros. The target is an implication to RH through a mathematically recognized criterion.

The relevant source-level endpoint is Weil positivity for the explicit formula. In the operator-theoretic program used by the project, the global/semi-local positivity problem is the bridge to RH.

The project therefore uses the following top-level dependency graph:

```text
RH
^ 
|  SOURCE THEOREM / Weil criterion
|
global Weil positivity for the required test-function class
^
|  OPEN BRIDGE in this bench
|
semilocal operator inequality / trace positivity sufficient for that class
^
|  SOURCE CONJECTURE / source framework plus missing proof obligations
|
semilocal Sonin/prolate construction as finite S grows
^
|  SOURCE THEOREM for the construction; bench regressions below
|
finite-S Jacobi/prolate operators and their source-derived coefficients
^
|  EXACT + NUMERICAL EVIDENCE in this bench
|
finite compressions, q-series perturbations, spectral diagnostics
```

The important conclusion is that the current bottleneck specific to Riemann is **not** the generic numerical control of a tridiagonal resolvent. It is the missing implication between the semilocal prolate/Jacobi constructions that the bench can compute and a Weil-positivity statement strong enough to invoke the global RH criterion.

## 2. External source anchors

### 2.1 Weil positivity and RH

`SOURCE THEOREM`.

The operator-theoretic literature used by the project takes Weil positivity for the explicit formula as the criterion tied to RH. The bench does not re-prove that classical equivalence.

Primary background used here:

- Connes, Consani, Marcolli, *The Weil proof and the geometry of the adeles class space*.
- Connes, Consani, *The Scaling Hamiltonian*, arXiv:1910.14368.
- Connes, Consani, *Weil positivity and Trace formula, the archimedean place*, arXiv:2006.13771.

### 2.2 Archimedean positivity mechanism

`SOURCE THEOREM / SOURCE-BACKED CONSTRUCTION`, but **not sufficient by itself for RH**.

arXiv:2006.13771 analyzes the single archimedean place and expresses the difference between the Weil distribution and the Sonin trace through prolate spheroidal data, with Hermitian Toeplitz matrices controlling the difference.

The paper explicitly presents the general semilocal setting as the place where Weil positivity would imply RH. The archimedean calculation is therefore a source-locked building block, not the full global proof.

### 2.3 Semilocal sufficiency for bounded support

`SOURCE CONJECTURE`.

In *The Scaling Hamiltonian*, the semilocal operator-theoretic program proposes that a finite set of places containing the archimedean place and sufficiently many primes should suffice to establish the Weil inequality for test functions with a corresponding bounded multiplicative support.

This is the conceptual bridge that would allow a growing finite set `S` to address the global criterion, but it is conjectural in the cited source. The bench must not silently promote it to a theorem.

### 2.4 Semilocal prolate operators

`SOURCE THEOREM / SOURCE-BACKED CONSTRUCTION`.

arXiv:2310.18423 introduces a semilocal analogue of the prolate wave operator, proves stability properties of the semilocal Sonin space as the finite set of places grows, and connects the construction to the spectral-realization program for zeta zeros.

This source statement does **not** justify identifying the finite generalized crossing parameters computed by this bench with zeta zeros.

### 2.5 Semilocal Jacobi q-series

`SOURCE THEOREM` for the source formulas; `EXACT` for the finite consequences derived in the bench.

arXiv:2403.01247 associates Jacobi data to finite sets of places and treats `S={infinity,p}` with `q=1/p`. The moments, orthogonal polynomials, and Jacobi coefficients admit q-series expansions. The current Phase 4 first-order matrix derivative is derived from these source coefficients.

## 3. What Phase 3 actually establishes

Phase 3 reproduces the archimedean `Q epsilon` / prolate Toeplitz benchmark.

### Established in the bench

- `NUMERICAL EVIDENCE / SOURCE REGRESSION`: the published boundary contributions `t(0)..t(4)` are reproduced in `tests/prolate_boundary_probe.rs`.
- `NUMERICAL EVIDENCE / SOURCE REGRESSION`: `epsilon'(1+) ~= 22.9965` is reproduced.
- `EXACT IMPLEMENTATION IDENTITY`: the `Q epsilon` implementation enforces the source identity `Q epsilon(1)=0` termwise through the source formula used by the regression.
- `NUMERICAL EVIDENCE / SOURCE REGRESSION`: the 11-mode source construction fed into the Toeplitz discretization reproduces approximately `1.05177` and `0.687925` in `tests/qepsilon_spectrum.rs`.

### Not yet established as an independent theorem in the bench

- an independent proof of the source's uniform `1e-11` remainder for the 11-term truncation;
- a theorem passing from the finite Toeplitz compression to the full semilocal/global Weil criterion;
- any direct implication from the reproduced two finite eigenvalues to RH.

Thus Phase 3 is a successful source reproduction, not a proof bridge.

## 4. What Phase 4 actually establishes

### Finite algebra already established

`EXACT` at fixed finite block size, modulo standard floating-point evaluation of the closed formulas:

- exact archimedean Jacobi coefficients;
- source-derived first q derivative at `q=0`;
- exact finite tridiagonal `K(0)` and `K'(0)` formulas;
- ordinary symmetric perturbation formulas for finite simple eigenvalues;
- finite sign lemma for the first derivative of the generalized crossings;
- exact ground-state factorisation identities;
- exact zero-shift left Schur-cavity closed form;
- exact first shift derivative of that cavity.

### Numerically validated large-m structure

`NUMERICAL EVIDENCE` and, where explicitly documented, `FORMAL ASYMPTOTIC`:

- high-dimensional pairwise eigenvalue response;
- total-variation checkpoints through `m=16384`;
- soft-edge fixed-rank scaling;
- candidate/logarithmic trace-growth laws;
- finite-size corrections;
- local frozen-cavity asymptotics.

### What Phase 4 does not establish

- finite generalized-prolate crossings = zeta zeros;
- monotonicity of these finite crossings => Weil positivity;
- sign of the first q derivative => positivity for an enlarged semilocal set;
- convergence of the sequence of finite compressions to a global RH criterion;
- sufficiency of the `S={infinity,p}` first-order perturbation for arbitrary compactly supported Weil test functions.

These missing arrows dominate the logical distance to RH.

## 5. Why the generic resolvent program leaves RiemannBench

The work on generic Jacobi resolvents, Schur cavities, frozen Toeplitz models, Green functions, and uniform slowly-varying tridiagonal estimates is mathematically useful, but its core statement can be formulated without mentioning zeta or RH.

It therefore belongs to the TDI-10.xx operator/resolvent research line once extracted generically.

RiemannBench may retain:

- its specific semilocal coefficients;
- adapters that instantiate a generic theorem;
- checks that the Riemann-specific coefficients satisfy the theorem hypotheses;
- the Riemann-specific consequence of that theorem.

It should not continue extending a general tridiagonal-operator theory merely because that theory was first encountered here.

## 6. The main Riemann-specific open bridge

The next central question is:

> What exact semilocal operator inequality, positivity statement, or convergence theorem would be sufficient to imply Weil positivity on a declared test-function support window, and how does the current semilocal prolate construction enter that statement?

This must be answered before more effort is spent on secondary asymptotics.

A useful target should have the shape:

```text
Given a finite set S and a declared support window I(S),
prove that operator statement P(S) implies
Weil(h) >= 0 for every admissible h supported in I(S).
```

Then one must separately establish either:

```text
P(S) for the required growing family of S,
```

or a theorem passing to an exhausting limit of support windows.

Until such a proposition is source-backed or proved, the following implication is an `OPEN BRIDGE`:

```text
semilocal prolate spectral behavior
        =>
Weil positivity on the matching support window.
```

## 7. Immediate research priorities

### R1 — source-lock the semilocal Weil bridge

Highest priority.

Extract from the primary sources the precise function spaces, support conditions, finite set of places, cutoff projections, Sonin space, and operator positivity statement that participate in the proposed semilocal proof route.

Deliverable: a mathematical specification with no numerical surrogate.

### R2 — identify theorem vs conjecture at every arrow

Highest priority.

In particular, do not treat the bounded-support semilocal sufficiency statement as proved if the source labels it conjectural.

### R3 — map existing code to the source bridge

For each implemented object (`Q epsilon`, Toeplitz form, semilocal Jacobi matrix, prolate operator, generalized crossing), state exactly which source object it represents and whether it occurs in a known implication toward Weil positivity.

Any object with no demonstrated role in the bridge is diagnostic, not a proof-critical target.

### R4 — define the first falsifiable Riemann-specific proposition

The next new theorem attempted in this repository should be Riemann-specific. It should connect an implemented semilocal quantity to a recognized Weil-positivity obligation.

If that proposition reduces to a generic Jacobi/resolvent theorem, extract the generic theorem to TDI-10.xx and keep only the application here.

## 8. Extraction rule for future work

For every new algorithm or lemma discovered while working on RH, ask:

```text
Can the statement and its tests be written without mentioning
Riemann, zeta, Weil, Sonin, the semilocal place set S,
or the specific prolate construction?
```

If yes, it is a candidate for transfer to the appropriate research bench.

Current default routing:

- generic Jacobi/resolvent/soft-edge operator theory -> TDI-10.xx;
- mature generic scientific numerical kernels -> SciRust;
- generic algorithm search/optimization machinery -> Forge;
- CUDA/NVIDIA kernels -> NNIS;
- orchestration/provenance -> SciRust Hub.

RiemannBench keeps the Riemann-specific adapter, hypotheses, provenance, and consequence.

## 9. Stop conditions

The project must stop or redirect a line of work when one of the following occurs:

1. the next lemma has no documented implication toward a recognized RH criterion;
2. the same lemma can be stated generically and has a natural owner elsewhere;
3. a numerical signature survives more precision but no logical bridge is found;
4. a source conjecture is being treated operationally as if it were a theorem;
5. a finite compression is being interpreted as a zeta-zero realization without a proof of that identification.

## 10. Current conclusion

The bench has strong source reproduction and finite-dimensional operator mathematics, but the logical gap to RH remains large.

The highest-value Riemann-specific task is now to formalize the semilocal-to-Weil bridge and determine exactly which operator statement would be sufficient for the global criterion.

The generic uniform-resolvent problem is no longer the central roadmap item of `riemann_ndim_bench`; it is an extracted supporting problem for TDI-10.xx.

Nothing in this audit asserts a proof of RH.
