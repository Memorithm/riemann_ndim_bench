#!/usr/bin/env python3
"""Install the exact symbolic mu^2 verifier into the research-agent tool surface.

The base runner remains unchanged. This module extends its existing
``verify_math`` tool with fail-closed modes that dispatch only to
``symbolic_mu2.py`` using fixed argv construction.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType


SYMBOLIC_MODES = {
    "symbolic_forcing_ratio",
    "symbolic_hypergeometric",
    "symbolic_finite_part",
}


def _bounded_int(args: dict, name: str, *, default: int, lower: int, upper: int) -> int:
    value = args.get(name, default)
    if not isinstance(value, int) or not lower <= value <= upper:
        raise RuntimeError(f"{name} must be an integer between {lower} and {upper}")
    return value


def install(base: ModuleType) -> None:
    """Idempotently extend ``riemann_research_agent`` in-place."""
    if getattr(base, "_SYMBOLIC_MU2_EXTENSION_INSTALLED", False):
        return

    original_verify_math = base.TOOL_IMPL["verify_math"]
    script = Path(base.__file__).resolve().with_name("symbolic_mu2.py")

    def extended_verify_math(root: Path, args: dict) -> str:
        mode = args.get("mode", "")
        if mode not in SYMBOLIC_MODES:
            return original_verify_math(root, args)

        argv = [sys.executable, str(script)]

        if mode == "symbolic_forcing_ratio":
            a_expr = base.text_arg(args, "A", max_len=2048)
            b_expr = base.text_arg(args, "B", max_len=2048)
            forcing = base.text_arg(args, "forcing", max_len=2048)
            offset = base.text_arg(args, "offset", max_len=128)
            candidate = base.text_arg(args, "candidate_ratio", max_len=2048)
            argv += [
                "forcing-ratio",
                "--A",
                a_expr,
                "--B",
                b_expr,
                "--forcing",
                forcing,
                "--offset",
                offset,
                "--candidate-ratio",
                candidate,
            ]

        elif mode == "symbolic_hypergeometric":
            # A candidate ratio is mandatory here. Without it the symbolic
            # helper can derive a Pochhammer quotient, but the agent gate must
            # audit the agent's *proposed* quotient rather than merely prove
            # that some supplied shifts define a hypergeometric sequence.
            numerator = base.text_arg(args, "numerator_shifts", max_len=1024)
            denominator = base.text_arg(args, "denominator_shifts", max_len=1024)
            base_value = base.text_arg(
                args,
                "hypergeometric_base",
                required=False,
                default="1",
                max_len=128,
            )
            candidate = base.text_arg(args, "candidate_ratio", max_len=2048)
            terms = _bounded_int(args, "terms", default=6, lower=1, upper=16)
            argv += [
                "hypergeometric",
                "--numerator-shifts",
                numerator,
                "--denominator-shifts",
                denominator,
                "--base",
                base_value,
                "--candidate-ratio",
                candidate,
                "--terms",
                str(terms),
            ]

        elif mode == "symbolic_finite_part":
            base_expr = base.text_arg(args, "base_expr", max_len=4096)
            theta_polynomial = base.text_arg(
                args,
                "theta_polynomial",
                required=False,
                default="1",
                max_len=2048,
            )
            extra_expr = base.text_arg(
                args,
                "extra_expr",
                required=False,
                default="0",
                max_len=4096,
            )
            order = _bounded_int(args, "order", default=10, lower=3, upper=20)
            argv += [
                "finite-part",
                "--base-expr",
                base_expr,
                "--theta-polynomial",
                theta_polynomial,
                "--extra-expr",
                extra_expr,
                "--order",
                str(order),
            ]

        return base.run_process(argv, root, timeout=60)

    base.TOOL_IMPL["verify_math"] = extended_verify_math

    mode_enum = base.VERIFY_PROPERTIES["mode"]["enum"]
    for mode in sorted(SYMBOLIC_MODES):
        if mode not in mode_enum:
            mode_enum.append(mode)

    base.VERIFY_PROPERTIES.update(
        {
            "forcing": {"type": "string"},
            "offset": {"type": "string"},
            "numerator_shifts": {"type": "string"},
            "denominator_shifts": {"type": "string"},
            "hypergeometric_base": {"type": "string"},
            "candidate_ratio": {"type": "string"},
            "terms": {"type": "integer"},
            "base_expr": {"type": "string"},
            "theta_polynomial": {"type": "string"},
            "extra_expr": {"type": "string"},
            "order": {"type": "integer"},
        }
    )

    for tool in base.TOOLS:
        function = tool.get("function", {})
        if function.get("name") == "verify_math":
            function["description"] += (
                " symbolic_forcing_ratio derives the variation-of-constants "
                "forcing quotient directly from exact recurrence coefficients; "
                "symbolic_hypergeometric exactly checks a proposed "
                "Pochhammer/hypergeometric coefficient quotient; "
                "symbolic_finite_part exactly applies a theta polynomial and "
                "extracts the x=sqrt(1-z) finite part by Puiseux expansion."
            )
            break

    base.SYSTEM_PROMPT += r"""

12. After the exact second-order recurrence has been extracted, do not jump
    directly to a familiar binomial or hypergeometric family. For each parity or
    shifted lattice under study, derive the candidate normalized forcing
    quotient and audit it first with verify_math symbolic_forcing_ratio using
    the verified A, B, forcing F and exact site offset. A MISMATCH or REFUTED
    verdict is a hard contradiction.
13. Only after the forcing quotient is recurrence-verified may you identify a
    Pochhammer/hypergeometric coefficient family. Audit that proposed family
    with verify_math symbolic_hypergeometric.
14. If a generating-function expression is derived, use verify_math
    symbolic_finite_part to apply the theta polynomial and extract the local
    finite part in x=sqrt(1-z). This proves only the supplied symbolic
    expression; separately justify why that expression follows from the source
    recurrence and the verified forcing sequence.
"""

    base._SYMBOLIC_MU2_EXTENSION_INSTALLED = True
