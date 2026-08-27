#!/usr/bin/env python3

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from proof_ledger import EvidenceStatus, ProofLedger


class ProofLedgerTests(unittest.TestCase):
    def test_exact_perturbative_success_satisfies_gate(self) -> None:
        ledger = ProofLedger()
        record = ledger.add_verifier_output(
            "perturbative_recurrence",
            "\n".join(
                [
                    "mode=perturbative_recurrence",
                    "mu1_status=PROVED_EXACT_SOLUTION",
                    "mu2_forcing_rhs=32/3*j*n",
                    "exact_status=PROVED_BY_FORMAL_COEFFICIENT_EXTRACTION",
                ]
            ),
        )
        self.assertEqual(record.status, EvidenceStatus.PROVED_EXACT)
        self.assertTrue(ledger.has_successful_perturbative_extraction())
        self.assertEqual(
            ledger.gate_failures(require_perturbative_success=True),
            [],
        )

    def test_failed_first_order_candidate_is_refuted(self) -> None:
        ledger = ProofLedger()
        record = ledger.add_verifier_output(
            "perturbative_recurrence",
            "\n".join(
                [
                    "mu1_status=CANDIDATE_U_FAILS",
                    "exact_status=PROVED_BY_FORMAL_COEFFICIENT_EXTRACTION",
                ]
            ),
        )
        self.assertEqual(record.status, EvidenceStatus.REFUTED)
        self.assertFalse(ledger.has_successful_perturbative_extraction())

    def test_unresolved_gamma_is_not_exact(self) -> None:
        ledger = ProofLedger()
        record = ledger.add_verifier_output(
            "gamma_quotient",
            "exact_status=UNRESOLVED_GAMMA_BASES\nunresolved=Gamma(1/3)^-1",
        )
        self.assertEqual(record.status, EvidenceStatus.UNRESOLVED)
        self.assertTrue(ledger.unresolved_gamma_seen())
        self.assertFalse(ledger.has_exact_success("gamma_quotient"))

    def test_asymptotic_best_power_is_preserved(self) -> None:
        ledger = ProofLedger()
        record = ledger.add_verifier_output(
            "asymptotic_power",
            "best_power=1/2\nwarning=ASYMPTOTIC_FIT_IS_NUMERICAL_EVIDENCE_NOT_PROOF",
        )
        self.assertEqual(record.status, EvidenceStatus.ASYMPTOTIC_EVIDENCE)
        self.assertEqual(ledger.best_asymptotic_power(), "1/2")

    def test_required_modes_fail_closed(self) -> None:
        ledger = ProofLedger()
        failures = ledger.gate_failures(
            required_modes={"gamma_quotient", "asymptotic_power"},
            require_perturbative_success=True,
            require_index_transform=True,
        )
        self.assertEqual(len(failures), 3)


if __name__ == "__main__":
    unittest.main()
