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

    def test_exact_mode_gate_rejects_unresolved_gamma(self) -> None:
        ledger = ProofLedger()
        ledger.add_verifier_output(
            "gamma_quotient",
            "exact_status=UNRESOLVED_GAMMA_BASES\nunresolved=Gamma(1/3)^-1",
        )
        failures = ledger.gate_failures(
            required_modes={"gamma_quotient"},
            require_exact_modes={"gamma_quotient"},
        )
        self.assertEqual(len(failures), 1)
        self.assertIn("without an exact successful result", failures[0])

    def test_exact_mode_gate_accepts_proved_gamma(self) -> None:
        ledger = ProofLedger()
        ledger.add_verifier_output(
            "gamma_quotient",
            "exact_status=PROVED_BY_GAMMA_RECURRENCE_AND_SPECIAL_IDENTITIES\nunresolved=none",
        )
        self.assertEqual(
            ledger.gate_failures(
                required_modes={"gamma_quotient"},
                require_exact_modes={"gamma_quotient"},
            ),
            [],
        )

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

    def test_symbolic_ratio_mismatch_is_refuted(self) -> None:
        ledger = ProofLedger()
        record = ledger.add_verifier_output(
            "symbolic_hypergeometric",
            "\n".join(
                [
                    "candidate_status=MISMATCH",
                    "exact_status=REFUTED_CANDIDATE_RATIO",
                ]
            ),
        )
        self.assertEqual(record.status, EvidenceStatus.REFUTED)
        self.assertFalse(ledger.has_exact_success("symbolic_hypergeometric"))

    def test_forcing_provenance_mismatch_is_refuted(self) -> None:
        ledger = ProofLedger()
        record = ledger.add_verifier_output(
            "symbolic_forcing_ratio",
            "\n".join(
                [
                    "candidate_status=MISMATCH",
                    "exact_status=REFUTED_FORCING_QUOTIENT",
                ]
            ),
        )
        self.assertEqual(record.status, EvidenceStatus.REFUTED)
        self.assertFalse(ledger.has_exact_success("symbolic_forcing_ratio"))

    def test_component_assembly_mismatch_is_refuted(self) -> None:
        ledger = ProofLedger()
        record = ledger.add_verifier_output(
            "symbolic_assembly",
            "\n".join(
                [
                    "assembled_value=sqrt(pi)/6",
                    "candidate_status=MISMATCH",
                    "exact_status=REFUTED_COMPONENT_ASSEMBLY",
                ]
            ),
        )
        self.assertEqual(record.status, EvidenceStatus.REFUTED)
        self.assertFalse(ledger.has_exact_success("symbolic_assembly"))

    def test_complete_symbolic_chain_satisfies_exact_gate(self) -> None:
        ledger = ProofLedger()
        ledger.add_verifier_output(
            "symbolic_forcing_ratio",
            "\n".join(
                [
                    "candidate_status=PROVED_EQUAL",
                    "derived_ratio=(n+1)/(n+2)",
                    "exact_status=PROVED_BY_VARIATION_OF_CONSTANTS_QUOTIENT",
                ]
            ),
        )
        ledger.add_verifier_output(
            "symbolic_hypergeometric",
            "\n".join(
                [
                    "candidate_status=PROVED_EQUAL",
                    "exact_status=PROVED_BY_POCHHAMMER_QUOTIENT",
                ]
            ),
        )
        ledger.add_verifier_output(
            "symbolic_finite_part",
            "\n".join(
                [
                    "finite_part=-sqrt(2)/9",
                    "exact_status=PROVED_BY_EXACT_THETA_ALGEBRA_AND_PUISEUX_SERIES",
                ]
            ),
        )
        ledger.add_verifier_output(
            "symbolic_assembly",
            "\n".join(
                [
                    "assembled_value=sqrt(pi)/6",
                    "candidate_status=PROVED_EQUAL",
                    "exact_status=PROVED_BY_EXACT_COMPONENT_ASSEMBLY",
                ]
            ),
        )
        self.assertTrue(ledger.has_exact_symbolic_mu2_chain())
        self.assertEqual(
            ledger.gate_failures(
                required_modes={
                    "symbolic_forcing_ratio",
                    "symbolic_hypergeometric",
                    "symbolic_finite_part",
                    "symbolic_assembly",
                },
                require_exact_modes={
                    "symbolic_forcing_ratio",
                    "symbolic_hypergeometric",
                    "symbolic_finite_part",
                    "symbolic_assembly",
                },
                require_symbolic_mu2_chain=True,
            ),
            [],
        )

    def test_chain_without_forcing_provenance_fails_closed(self) -> None:
        ledger = ProofLedger()
        ledger.add_verifier_output(
            "symbolic_hypergeometric",
            "candidate_status=PROVED_EQUAL\nexact_status=PROVED_BY_POCHHAMMER_QUOTIENT",
        )
        ledger.add_verifier_output(
            "symbolic_finite_part",
            "finite_part=0\nexact_status=PROVED_BY_EXACT_THETA_ALGEBRA_AND_PUISEUX_SERIES",
        )
        ledger.add_verifier_output(
            "symbolic_assembly",
            "assembled_value=0\ncandidate_status=PROVED_EQUAL\nexact_status=PROVED_BY_EXACT_COMPONENT_ASSEMBLY",
        )
        self.assertFalse(ledger.has_exact_symbolic_mu2_chain())
        failures = ledger.gate_failures(require_symbolic_mu2_chain=True)
        self.assertEqual(len(failures), 1)
        self.assertIn("symbolic mu2 chain is incomplete", failures[0])
        self.assertIn("symbolic_forcing_ratio", failures[0])

    def test_chain_without_final_assembly_fails_closed(self) -> None:
        ledger = ProofLedger()
        ledger.add_verifier_output(
            "symbolic_forcing_ratio",
            "candidate_status=PROVED_EQUAL\nexact_status=PROVED_BY_VARIATION_OF_CONSTANTS_QUOTIENT",
        )
        ledger.add_verifier_output(
            "symbolic_hypergeometric",
            "candidate_status=PROVED_EQUAL\nexact_status=PROVED_BY_POCHHAMMER_QUOTIENT",
        )
        ledger.add_verifier_output(
            "symbolic_finite_part",
            "finite_part=0\nexact_status=PROVED_BY_EXACT_THETA_ALGEBRA_AND_PUISEUX_SERIES",
        )
        self.assertFalse(ledger.has_exact_symbolic_mu2_chain())
        failures = ledger.gate_failures(require_symbolic_mu2_chain=True)
        self.assertEqual(len(failures), 1)
        self.assertIn("symbolic_assembly", failures[0])


if __name__ == "__main__":
    unittest.main()
