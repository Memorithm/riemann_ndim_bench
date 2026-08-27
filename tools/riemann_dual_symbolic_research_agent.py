#!/usr/bin/env python3
"""Dual-agent entrypoint with exact post-perturbative mu^2 gates."""

from __future__ import annotations

import riemann_research_agent as base
from symbolic_mu2_agent_extension import install


install(base)

import riemann_dual_research_agent as dual


SYMBOLIC_REQUIRED_MODES = {
    "symbolic_hypergeometric",
    "symbolic_finite_part",
}

dual.FINAL_REQUIRED_MODES.update(SYMBOLIC_REQUIRED_MODES)
dual.FINAL_REQUIRED_EXACT_MODES.update(SYMBOLIC_REQUIRED_MODES)

dual.RESEARCHER_SYSTEM += r"""

Post-perturbative symbolic requirements:
- once the exact mu^2 forcing is established, derive any parity/subsequence
  quotient from that forcing rather than guessing a familiar binomial family;
- audit the proposed coefficient quotient with verify_math
  symbolic_hypergeometric;
- derive the generating function and theta-polynomial weighting from the same
  forcing sequence;
- audit the local x=sqrt(1-z) finite part with verify_math
  symbolic_finite_part;
- keep separate algebraic verification of a supplied expression from the proof
  that the expression was actually derived from the recurrence.
"""

dual.CRITIC_SYSTEM += r"""

Independently audit the post-perturbative chain. Reject a hypergeometric family
recognized only numerically or by pattern matching. Require an exact
symbolic_hypergeometric quotient check, then inspect whether the generating
expression really follows from that sequence. A finite-part tool success proves
the supplied expression only, not its provenance.
"""

dual.FINAL_SYSTEM += r"""

The final synthesis is additionally gated on exact successes for
symbolic_hypergeometric and symbolic_finite_part. These gates do not license a
provenance shortcut: state explicitly how the coefficient sequence and weighted
generating expression were derived from the verified mu^2 forcing.
"""


if __name__ == "__main__":
    dual.main()
