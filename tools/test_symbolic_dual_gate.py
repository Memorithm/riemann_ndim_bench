#!/usr/bin/env python3
"""Regression tests for the symbolic dual-agent final gate."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import riemann_dual_symbolic_research_agent  # installs the extension and gate
import riemann_dual_research_agent as dual
from proof_ledger import ProofLedger


class SymbolicDualGateTests(unittest.TestCase):
    def base_final_ledger(self) -> ProofLedger:
        ledger = ProofLedger()
        ledger.add_verifier_output(
            "recurrence_transform",
            "exact_status=PROVED_BY_EXACT_AFFINE_SUBSTITUTION_AND_SIGN_NORMALIZATION",
        )
        ledger.add_verifier_output(
            "perturbative_recurrence",
            "\n".join(
                [
                    "mu1_status=PROVED_EXACT_SOLUTION",
                    "mu2_forcing_rhs=32/3*j*n",
                    "exact_status=PROVED_BY_FORMAL_COEFFICIENT_EXTRACTION",
                ]
            ),
        )
        ledger.add_verifier_output(
            "asymptotic_power",
            "best_power=1/2\nwarning=ASYMPTOTIC_FIT_IS_NUMERICAL_EVIDENCE_NOT_PROOF",
        )
        ledger.add_verifier_output(
            "gamma_quotient",
            "exact_status=PROVED_BY_GAMMA_RECURRENCE_AND_SPECIAL_IDENTITIES",
        )
        return ledger

    def test_symbolic_modes_are_required(self) -> None:
        self.assertIn("symbolic_hypergeometric", dual.FINAL_REQUIRED_MODES)
        self.assertIn("symbolic_finite_part", dual.FINAL_REQUIRED_MODES)
        self.assertIn("symbolic_hypergeometric", dual.FINAL_REQUIRED_EXACT_MODES)
        self.assertIn("symbolic_finite_part", dual.FINAL_REQUIRED_EXACT_MODES)

    def test_gate_rejects_missing_symbolic_chain(self) -> None:
        failures = dual.final_gate_failures(self.base_final_ledger())
        rendered = "\n".join(failures)
        self.assertIn("symbolic_hypergeometric", rendered)
        self.assertIn("symbolic_finite_part", rendered)

    def test_gate_accepts_complete_exact_symbolic_chain(self) -> None:
        ledger = self.base_final_ledger()
        ledger.add_verifier_output(
            "symbolic_hypergeometric",
            "candidate_status=PROVED_EQUAL\nexact_status=PROVED_BY_POCHHAMMER_QUOTIENT",
        )
        ledger.add_verifier_output(
            "symbolic_finite_part",
            "finite_part=0\nexact_status=PROVED_BY_EXACT_THETA_ALGEBRA_AND_PUISEUX_SERIES",
        )
        self.assertEqual(dual.final_gate_failures(ledger), [])

    def test_gate_rejects_refuted_symbolic_candidate(self) -> None:
        ledger = self.base_final_ledger()
        ledger.add_verifier_output(
            "symbolic_hypergeometric",
            "candidate_status=MISMATCH\nexact_status=REFUTED_CANDIDATE_RATIO",
        )
        ledger.add_verifier_output(
            "symbolic_finite_part",
            "finite_part=0\nexact_status=PROVED_BY_EXACT_THETA_ALGEBRA_AND_PUISEUX_SERIES",
        )
        failures = dual.final_gate_failures(ledger)
        self.assertTrue(
            any("without an exact successful result" in failure for failure in failures)
        )


if __name__ == "__main__":
    unittest.main()
