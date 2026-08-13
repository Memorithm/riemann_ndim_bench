# Phase 3 checklist — prolate and Q-epsilon

Primary source: Connes–Consani, arXiv:2006.13771.

- [ ] Fix the exact normalization of the prolate eigenvalues lambda(n) and vectors xi_n, zeta_n used in equation (84).
- [ ] Reproduce the analytic continuation xi_n^an used in equation (99).
- [ ] Reproduce the first published boundary contributions t(0)..t(4).
- [ ] Reproduce epsilon'(1+) approximately 22.9965 with a stated numerical tolerance.
- [ ] Verify Q-epsilon(1) = 0.
- [ ] Implement Q-epsilon from equation (99), not from a fitted surrogate.
- [ ] Check the 11-term truncation against the published uniform 1e-11 approximation statement.
- [ ] Feed the verified kernel into the existing q^Z Toeplitz discretization.
- [ ] At omega=1e-3 and a=log(2), reproduce the published largest eigenvalue approximately 1.05177.
- [ ] Reproduce the published second eigenvalue approximately 0.687925.
- [ ] Only after all upstream regressions pass, study convergence as q tends to 1.
- [ ] Keep the experimental pi-radial coordinate logically separate and test whether it contributes any invariant not already encoded by the completed zeta factors.

No successful item in this checklist constitutes a proof of RH; each is an implementation/regression milestone.
