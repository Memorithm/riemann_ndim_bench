# Semilocal Weil bridge — source-locked specification

Status: Riemann-specific research specification.

Purpose: state exactly which source theorems, conjectures and operator objects connect the current semilocal/prolate bench to Weil positivity, and therefore to the Riemann hypothesis (RH). This file deliberately contains no new theorem.

## Status labels

- `SOURCE THEOREM` — proved in the cited literature.
- `SOURCE FACT` — explicit finite/logical consequence stated in the cited source.
- `SOURCE CONJECTURE` — conjecture in the cited literature; not promoted here.
- `SOURCE STRATEGY` — a proposed route or expectation, not a theorem.
- `BENCH REGRESSION` — source quantity numerically reproduced in this repository.
- `OPEN IMPLEMENTATION` — source object not yet represented faithfully in the bench.
- `OPEN BRIDGE` — mathematical implication still missing.
- `REFUTED ROUTE` — proposed route explicitly shown to fail in the cited source.

## 1. Weil criterion: the actual RH endpoint

Primary source for the formulation used here:

- Alain Connes, Caterina Consani, *The Scaling Hamiltonian*, arXiv:1910.14368, §4.
- Alain Connes, Caterina Consani, *Weil positivity and Trace formula, the archimedean place*, arXiv:2006.13771, Introduction and Appendix C.

`SOURCE THEOREM`.

The cited papers recall Weil's criterion: RH is equivalent to the required sign of the Riemann-Weil explicit formula on convolution squares subject to the two vanishing conditions. In the symmetric multiplicative notation of *The Scaling Hamiltonian*, one works with

```text
h = h1 * h2,
h2(x) = overline(h1(x^-1)),
```

and the boundary conditions

```text
integral_0^infinity h1(x) x^(+1/2) d*x = 0,
integral_0^infinity h1(x) x^(-1/2) d*x = 0.
```

The source also recalls that compactly supported test functions suffice.

This is the proof endpoint that matters to RiemannBench. Finite prolate eigenvalues, generalized crossing parameters and asymptotic constants are not substitutes for this criterion.

## 2. Compact support implies finitely many places

Primary source: *The Scaling Hamiltonian*, Fact 3.2.

`SOURCE FACT`.

For a compactly supported test function, the Riemann-Weil explicit formula involves only finitely many places: the contribution of a prime depends on values at non-zero powers of that prime, so all sufficiently large primes are absent for a fixed compact support.

This fact is the logical reason a semilocal finite-place calculation can be relevant to a global RH criterion one support window at a time.

It does **not** by itself prove that the available semilocal Hilbert-space operator yields the required Weil sign.

## 3. Semilocal trace formula

Primary source: *The Scaling Hamiltonian*, Theorem 2.5 and equations (2.28)–(2.30).

`SOURCE THEOREM`.

Let `S` be a finite set of places containing the archimedean place. The source constructs the semilocal Hilbert space `L^2(X_S)`, infrared and ultraviolet cutoffs, and a cutoff operator `R_Lambda`. For compactly supported `f`, Theorem 2.5 gives, as `Lambda -> infinity`, a trace formula of the form

```text
Tr(theta_a(f) R_Lambda)
  = 2 f(1) log Lambda
    + sum_{v in S} local_v(f)
    + o(1).
```

In the symmetric `h` normalization, equation (2.28) expresses the local sum in terms of the quantized logarithmic differential of the product of the local unitaries.

This is a genuine operator-theoretic representation of the required finite-place distributions. It is a source theorem and must be kept distinct from the finite Jacobi matrices currently used by Phase 4.

### Current bench status

`OPEN IMPLEMENTATION`.

RiemannBench currently does **not** expose the source semilocal Hilbert space `L^2(X_S)`, its cutoff operators, or Theorem 2.5 as a verified numerical/operator object. The current semilocal Jacobi/prolate matrices encode a different, source-related representation and must not silently be identified with `R_Lambda`.

## 4. A tempting direct sign route is known to fail

Primary source: *The Scaling Hamiltonian*, equations (3.1), (3.2), Lemma 3.4 and Fact 3.6.

`REFUTED ROUTE`.

The paper examines the natural attempt to prove the finite-place Weil inequality by forcing a constant sign on the logarithmic differential operator appearing in equation (3.2). Fact 3.6 states explicitly that inequality (3.1) does **not** hold in general; the single archimedean place already gives a counterexample through the non-monotonicity of the Riemann-Siegel angular function.

Therefore the following roadmap item is forbidden in this repository unless the hypotheses are materially changed and proved:

```text
prove a global constant sign for u^-1 du
=> conclude Weil positivity.
```

A numerical sign pattern in the finite generalized-prolate crossings cannot repair a source-level counterexample to the naïve operator inequality.

## 5. The precise semilocal support conjecture

Primary source: *The Scaling Hamiltonian*, Conjecture 4.1.

`SOURCE CONJECTURE`.

For a parameter `q`, the source defines

```text
S(q) = {infinity} union { primes p with p < q }.
```

Conjecture 4.1 states that the semilocal operator-theoretic framework for `S(q)` suffices to prove the Weil inequality for all test functions supported in

```text
(q^(-1/2), q^(1/2)).
```

This is the cleanest source statement connecting a growing finite set of places to growing compact-support windows.

It is a **conjecture**, not an established implication. RiemannBench must never write

```text
semilocal positivity for S(q) => RH
```

without specifying whether the step being used is proved independently or is invoking Conjecture 4.1 conditionally.

The same section emphasizes the Poisson formula and the normalization of the basic additive characters as central to any solution of this conjecture.

## 6. Archimedean model: what is actually proved

Primary source: *Weil positivity and Trace formula, the archimedean place*, arXiv:2006.13771.

### 6.1 Positive Sonin trace

`SOURCE THEOREM`.

Let `S` denote the projection onto the archimedean Sonin space, the orthogonal complement of the cutoff ranges in position and Fourier space. For convolution squares, the trace

```text
Tr(theta(g) S theta(g)^*)
```

is non-negative by construction.

The source proves an identity of the form

```text
Tr(theta(f) S)
  = W_infinity(f) + integral f(rho) epsilon(rho) d*rho,
```

where `epsilon` is expressed through prolate spheroidal data.

This identity is the source of the `epsilon` / `Q epsilon` construction reproduced in Phase 3.

### 6.2 Support and boundary operator Q

`SOURCE THEOREM`.

The paper introduces

```text
Q = -(rho d/drho)^2 + 1/4,
```

which imposes the required vanishing at the two boundary Mellin/Fourier points while preserving compact support. Applying `Q` to the trace-remainder transforms the support-constrained positivity problem into a compact-operator problem.

This is conceptually much closer to the RH bridge than the generalized crossing parameters of Phase 4.

### 6.3 Compact Toeplitz obstruction

`SOURCE THEOREM` plus `BENCH REGRESSION`.

For the interval associated with support near a factor of two, the source reduces the remaining archimedean problem to the spectrum of a compact operator represented numerically by a Hermitian Toeplitz matrix. The first two reported eigenvalues are approximately

```text
1.05177
0.687925.
```

RiemannBench reproduces these values in `tests/qepsilon_spectrum.rs` and reproduces the boundary quantity

```text
epsilon'(1+) ~= 22.9965
```

and its leading contributions in `tests/prolate_boundary_probe.rs`.

Interpretation: the finite Toeplitz regression confirms the source's archimedean obstruction structure. It is **not** a zeta-zero computation and does not prove RH.

### 6.4 Archimedean lower bound

`SOURCE THEOREM`.

The source proves, for the declared compact support and vanishing conditions, a lower bound comparing the archimedean Weil functional with the positive Sonin trace, with at most an explicit rank-one/evaluation correction. With the additional vanishing at zero, this yields the required archimedean Weil positivity on that support window.

This is the model that a successful semilocal extension would need to generalize in substance, not merely reproduce spectrally.

## 7. Semilocal prolate/Sonin structure

Primary source:

- Alain Connes, Caterina Consani, Henri Moscovici, *Zeta zeros and prolate wave operators*, arXiv:2310.18423.

`SOURCE THEOREM` for the constructions and stability statements; `SOURCE STRATEGY` for their proposed role toward RH.

The paper introduces a semilocal analogue of the prolate wave operator and proves stability properties of the semilocal Sonin space when the finite place set grows. It also relates this structure to the spectral-realization program for zeta zeros.

What RiemannBench may safely use:

- the semilocal prolate operator is source-backed;
- the Sonin-space structure is source-backed;
- increasing finite place sets is mathematically built into the semilocal framework.

What RiemannBench may **not** infer without an additional proof:

- finite generalized crossings are zeta zeros;
- their first derivative under `q=1/p` controls the Weil functional;
- convergence of these finite spectra establishes the Weil criterion.

## 8. Semilocal Jacobi q-series

Primary source:

- Alain Connes, Caterina Consani, Henri Moscovici, *On q-series and the moment problem associated to local factors*, arXiv:2403.01247.

`SOURCE THEOREM` for the q-series coefficients; `EXACT` for the finite algebra derived from them in this bench.

For `S={infinity,p}`, with `q=1/p`, the source provides a q-series for the Jacobi coefficients. Phase 4 uses the first source coefficient to construct `K'(0)` and the finite generalized-crossing derivatives.

This is a legitimate semilocal deformation diagnostic.

The missing Riemann-specific arrow is:

```text
Jacobi/prolate q-deformation
      ?
      v
variation of the semilocal Weil/Sonin trace inequality.
```

No source theorem currently encoded in the bench supplies this arrow.

## 9. Map of source objects to current repository objects

| Source object / obligation | Status in RiemannBench | Current representation |
|---|---|---|
| Weil convolution-square criterion | documented only | `RH_PROOF_DEPENDENCY_AUDIT.md` |
| compact-support finite-place reduction | documented only | this specification |
| semilocal trace formula on `L^2(X_S)` | missing | `OPEN IMPLEMENTATION` |
| basic additive-character / Poisson normalization | missing | `OPEN IMPLEMENTATION` |
| archimedean prolate basis | reproduced in tests | `tests/support/prolate_basis.rs` |
| archimedean `epsilon'(1+)` | reproduced | `tests/prolate_boundary_probe.rs`, `tests/qepsilon_spectrum.rs` |
| archimedean `Q epsilon` kernel | reproduced, test-local | `tests/qepsilon_spectrum.rs` |
| source Toeplitz eigenvalues | reproduced | `tests/qepsilon_spectrum.rs` |
| reusable source-locked `Q epsilon` API | missing | should be upstreamed from test-only code |
| semilocal Jacobi coefficients | implemented | `src/semilocal.rs` |
| finite semilocal prolate compression | implemented | `src/semilocal.rs` |
| first q derivative at `q=0` | implemented / exact finite algebra | `src/semilocal.rs` + tests |
| generic resolvent/cavity machinery | implemented but generic | extraction target TDI-10.xx |
| semilocal Sonin projection/operator | missing | `OPEN IMPLEMENTATION` |
| semilocal counterpart of archimedean trace-remainder comparison | missing | `OPEN IMPLEMENTATION` |
| implementation of Conjecture 4.1 | impossible to mark complete without proof | `OPEN BRIDGE` |

## 10. Immediate implementation roadmap

### R1 — upstream the source-locked archimedean prolate / Q-epsilon kernel

Priority: high.

The current verified implementation lives only in integration-test support code. Move it into a focused Riemann-specific library module without changing formulas:

```text
prolate basis
-> epsilon'(1+)
-> C_n(rho)
-> Q epsilon(rho)
-> normalized Toeplitz kernel.
```

Keep independent regressions against the published boundary contributions and the two published Toeplitz eigenvalues.

Reason: the archimedean trace-remainder is a genuine node of the Weil bridge and should not remain hidden as test scaffolding.

### R2 — make the Weil support/boundary transform executable

Priority: high after R1.

Implement a source-locked representation of the multiplicative test-function conventions and the differential operator

```text
Q = -(rho d/drho)^2 + 1/4
```

with manufactured exact tests verifying the two boundary Mellin/Fourier zeros and support preservation.

Do not infer positivity from this transform alone.

### R3 — specify the semilocal trace object before implementing it

Priority: high.

Write the exact data model required by Theorem 2.5:

- finite set of places `S`;
- basic character normalization;
- semilocal test function;
- cutoff parameter `Lambda`;
- trace and local-distribution sides.

Only after this source-level API is fixed should implementation begin.

### R4 — test the first Riemann-specific bridge proposition

Candidate form:

```text
For a declared support window I and finite place set S,
show that a source-identified semilocal trace-remainder inequality
implies the Weil inequality on I.
```

The first deliverable may be a conditional proposition or a counterexample. It must not assume Conjecture 4.1 under a different name.

## 11. Explicit non-goals

The following do not count as direct progress on the RH bridge without a further theorem:

- pushing generalized crossings to larger matrix sizes;
- fitting additional constants to those crossings;
- proving generic resolvent estimates that do not mention the Weil/Sonin source objects;
- observing agreement with known zeta zeros;
- proving sign patterns only for finite compressed eigenvalue derivatives.

Such work may still be valuable, but generic components should be transferred to the appropriate bench and Riemann-specific diagnostics should remain subordinate to a documented proof dependency.

## 12. Current conclusion

The source-locked Riemann route is now:

```text
Weil criterion for compact convolution squares
        |
        v
compact support -> finitely many places
        |
        v
semilocal trace formula / normalized local factors
        |
        v
positive Sonin/cutoff trace + controlled trace remainder
        |
        v
semilocal support-window inequality
        |
        v
exhaust support windows
        |
        v
RH
```

The first two arrows are source theorems/facts. The archimedean instance of the positivity mechanism is source-proved and numerically reproduced in part by this bench. The general semilocal support-window sufficiency is Conjecture 4.1 and remains open.

Therefore the next code in RiemannBench should strengthen the faithful representation of the Weil/Sonin source objects, not continue generic resolvent theory.

Nothing in this specification claims a proof of RH.
