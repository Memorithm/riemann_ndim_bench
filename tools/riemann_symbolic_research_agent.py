#!/usr/bin/env python3
"""Single-agent entrypoint with exact post-perturbative mu^2 verification."""

from __future__ import annotations

import riemann_research_agent as base
from symbolic_mu2_agent_extension import install


install(base)


if __name__ == "__main__":
    base.main()
