#!/usr/bin/env python3
"""Dual-agent entrypoint with exact post-perturbative mu^2 gates."""

from __future__ import annotations

import riemann_research_agent as base
from symbolic_mu2_agent_extension import install


install(base)

import riemann_dual_research_agent as dual


SYMBOLIC_REQUIRED_MODES = {
    "symbolic_forcing_ratio",
    "symbolic_hypergeometric",
    "symbolic_finite_part",
    "symbolic_assembly",
}

dual.FINAL_REQUIRED_MODES.update(SYMBOLIC_REQUIRED_MODES)
dual.FINAL_REQUIRED_EXACT_MODES.update(SYMBOLIC_REQUIRED_MODES)

dual.RESEARCHER_SYSTEM += r"""

Post-perturbative symbolic requirements:
- once the exact mu^2 forcing is established, derive the normalized
  variation-of-constants quotient for each parity/shift directly from the
  verified recurrence and audit it with verify_math symbolic_forcing_ratio;
- only after that provenance check, identify any Pochhammer/hypergeometric
  family and audit its proposed coefficient quotient with verify_math
  symbolic_hypergeometric;
- derive the generating function and theta-polynomial weighting from the same
  recurrence-verified forcing sequence;
- audit the local x=sqrt(1-z) finite part with verify_math
  symbolic_finite_part;
- derive the Gamma normalization factors independently with exact Gamma tools;
- assemble those exact Gamma factors and finite parts only at the end, using
  verify_math symbolic_assembly to audit both the algebraic combination and the
  proposed closed form;
- keep separate algebraic verification of supplied expressions from the proof
  that those expressions were actually derived from the recurrence.
"""

dual.CRITIC_SYSTEM += r"""

Independently audit the post-perturbative chain. Reject a hypergeometric family
recognized only numerically or by pattern matching. First require an exact
symbolic_forcing_ratio check derived from A, B, F and the parity/site offset;
then require an exact symbolic_hypergeometric quotient check. Inspect whether
the generating expression really follows from that sequence. A finite-part
tool success proves the supplied expression only, not its provenance. Finally,
check that the exact Gamma normalizations and finite parts are combined with the
correct signs and normalization, and require an exact symbolic_assembly success
for any published closed-form residual constant.
"""

dual.FINAL_SYSTEM += r"""

The final synthesis is additionally gated on exact successes for
symbolic_forcing_ratio, symbolic_hypergeometric, symbolic_finite_part and
symbolic_assembly. These gates do not license a provenance shortcut: state
explicitly how the normalized forcing quotient follows from the verified mu^2
recurrence, how the coefficient family follows from that quotient, how the
weighted generating expression follows from the verified sequence, how the
Gamma factors were established, and why the final algebraic combination has the
claimed normalization and signs.
"""


if __name__ == "__main__":
    dual.main()
