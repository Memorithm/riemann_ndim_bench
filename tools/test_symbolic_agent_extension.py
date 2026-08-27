#!/usr/bin/env python3
"""Regression tests for symbolic_mu2_agent_extension."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import riemann_research_agent as base
from proof_ledger import EvidenceStatus, ProofLedger
from symbolic_mu2_agent_extension import install


ROOT = Path(__file__).resolve().parents[1]


class SymbolicAgentExtensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        install(base)

    def test_modes_are_registered(self) -> None:
        modes = base.VERIFY_PROPERTIES["mode"]["enum"]
        self.assertIn("symbolic_hypergeometric", modes)
        self.assertIn("symbolic_finite_part", modes)

    def test_hypergeometric_dispatch_proves_exact_candidate(self) -> None:
        name, output = base.execute_tool(
            ROOT,
            {
                "function": {
                    "name": "verify_math",
                    "arguments": {
                        "mode": "symbolic_hypergeometric",
                        "numerator_shifts": "1/4,3/4",
                        "denominator_shifts": "1/2,1",
                        "hypergeometric_base": "1",
                        "candidate_ratio": "((k+1/4)*(k+3/4))/((k+1/2)*(k+1))",
                        "terms": 5,
                    },
                }
            },
        )
        self.assertEqual(name, "verify_math")
        self.assertIn("exit_status=0", output)
        self.assertIn("candidate_status=PROVED_EQUAL", output)
        self.assertIn("exact_status=PROVED_BY_POCHHAMMER_QUOTIENT", output)

        ledger = ProofLedger()
        record = ledger.add_verifier_output("symbolic_hypergeometric", output)
        self.assertEqual(record.status, EvidenceStatus.PROVED_EXACT)

    def test_hypergeometric_dispatch_refutes_wrong_candidate(self) -> None:
        _, output = base.execute_tool(
            ROOT,
            {
                "function": {
                    "name": "verify_math",
                    "arguments": {
                        "mode": "symbolic_hypergeometric",
                        "numerator_shifts": "1/2",
                        "denominator_shifts": "1",
                        "candidate_ratio": "(k+3/2)/(k+1)",
                    },
                }
            },
        )
        self.assertIn("candidate_status=MISMATCH", output)
        self.assertIn("exact_status=REFUTED_CANDIDATE_RATIO", output)

        ledger = ProofLedger()
        record = ledger.add_verifier_output("symbolic_hypergeometric", output)
        self.assertEqual(record.status, EvidenceStatus.REFUTED)

    def test_finite_part_dispatch_is_exact(self) -> None:
        _, output = base.execute_tool(
            ROOT,
            {
                "function": {
                    "name": "verify_math",
                    "arguments": {
                        "mode": "symbolic_finite_part",
                        "base_expr": "1/sqrt(1-z)",
                        "theta_polynomial": "1",
                        "extra_expr": "0",
                        "order": 8,
                    },
                }
            },
        )
        self.assertIn("finite_part=0", output)
        self.assertIn(
            "exact_status=PROVED_BY_EXACT_THETA_ALGEBRA_AND_PUISEUX_SERIES",
            output,
        )

    def test_hypergeometric_requires_candidate_ratio(self) -> None:
        _, output = base.execute_tool(
            ROOT,
            {
                "function": {
                    "name": "verify_math",
                    "arguments": {
                        "mode": "symbolic_hypergeometric",
                        "numerator_shifts": "1/2",
                        "denominator_shifts": "1",
                    },
                }
            },
        )
        self.assertIn("TOOL ERROR", output)
        self.assertIn("candidate_ratio", output)


if __name__ == "__main__":
    unittest.main()
