# Phase 3 checklist — prolate and Q-epsilon

Primary source: Connes–Consani, arXiv:2006.13771.

This checklist tracks **source reproduction**, not progress toward a proof of RH.

Status legend:

- `[x]` reproduced by an explicit bench regression;
- `[~]` implemented and indirectly constrained by downstream source regressions, but still missing a focused independent acceptance test;
- `[ ]` open.

## Source normalization and prolate data

- [~] Fix the exact normalization of the prolate eigenvalues `lambda(n)` and vectors `xi_n`, `zeta_n` used around equation (84).
  - Current evidence: `tests/support/prolate_basis.rs` implements the source-normalized cosine integral operator and its modes.
  - The downstream boundary and Toeplitz spectrum regressions strongly constrain the normalization.
  - Remaining work: move the implementation out of test-only support code and add focused normalization regressions rather than declaring the normalization proved from downstream matches alone.
- [~] Reproduce the analytic continuation `xi_n^an` and derivative used in equation (99).
  - Current implementation: `Basis::value(mode, y)` and `Basis::derivative(mode, y)` evaluate the source continuation through the prolate integral representation, including `y>1` values used by `Q epsilon`.
  - Remaining work: dedicated source regression for the continuation itself.

## Boundary quantities

- [x] Reproduce the first published boundary contributions `t(0)..t(4)`.
  - Regression: `tests/prolate_boundary_probe.rs`.
  - Targets: `11.9719`, `8.77574`, `2.20528`, `0.0433983`, `0.000125459`.
- [x] Reproduce `epsilon'(1+) ~= 22.9965` with a stated numerical tolerance.
  - Regressions: `tests/prolate_boundary_probe.rs` and `tests/qepsilon_spectrum.rs`.
- [~] Verify `Q epsilon(1)=0`.
  - The current source-form implementation returns the source value zero at `rho=1` termwise in `tests/qepsilon_spectrum.rs`.
  - Remaining work: add a focused public-kernel regression after `Q epsilon` is moved into reusable source-locked code.

## Q-epsilon kernel

- [x] Implement `Q epsilon` from the source equation, not from a fitted surrogate.
  - Current implementation: `c_n(...)` and `q_epsilon(...)` in `tests/qepsilon_spectrum.rs`.
  - Limitation: the implementation still lives in integration-test code rather than a reusable library module.
- [ ] Independently check the first-11-term truncation against the source's published uniform `1e-11` remainder statement.
  - The current test uses 11 modes because the source supplies the bound.
  - The bench has **not** independently reproduced/proved that uniform remainder.

## Toeplitz source benchmark

- [x] Feed the source-derived kernel into the existing `q^Z` Toeplitz discretization.
  - Regression: `tests/qepsilon_spectrum.rs`.
- [x] At `omega=1e-3` and `a=log(2)`, reproduce the published largest eigenvalue `~= 1.05177`.
  - Regression tolerance currently `8e-5`.
- [x] Reproduce the published second eigenvalue `~= 0.687925`.
  - Regression tolerance currently `8e-5`.

## Still-open Phase 3 work

- [ ] Move the verified prolate / `Q epsilon` implementation from test scaffolding into a focused Riemann-specific library module while keeping an independent regression path.
- [ ] Add direct tests for `Q epsilon(1)=0` and the source continuation/normalization conventions.
- [ ] Independently reproduce or otherwise formally account for the source's uniform 11-term truncation error bound.
- [ ] Study convergence as the logarithmic discretization parameter tends to zero (`q -> 1+`) only after the source kernel API and its acceptance tests are stable.
- [ ] Connect the reproduced archimedean trace-remainder object to the source-locked semilocal Weil bridge described in `SEMILOCAL_WEIL_BRIDGE_SPEC.md`.

## Historical experimental branch

- [x] Keep the experimental pi-radial coordinate logically separate from the source-locked Weil/prolate chain.
  - Phase 4 controls showed that the tested scalar deformation was not specific to `pi`; it must not be promoted as a distinguished RH mechanism.

No successful item in this checklist constitutes a proof of RH. In particular, reproducing the Toeplitz eigenvalues validates the archimedean source calculation; it does not identify those eigenvalues, or any generalized-prolate crossings, with zeta zeros.
