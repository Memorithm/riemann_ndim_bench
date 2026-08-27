#!/usr/bin/env python3
"""Regression tests for tools/verify_math.py."""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("verify_math.py")


def run_verifier(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        text=True,
        capture_output=True,
        check=check,
        timeout=20,
    )


class VerifyMathTests(unittest.TestCase):
    def test_rational_exact(self) -> None:
        result = run_verifier("rational", "--expr", "9/8 + 1/8")
        self.assertIn("exact_value=5/4", result.stdout)
        self.assertIn("PROVED_BY_RATIONAL_ARITHMETIC", result.stdout)

    def test_gamma_ratio_reduces_exactly(self) -> None:
        result = run_verifier(
            "gamma-quotient",
            "--numerator",
            "9/4,1/2",
            "--denominator",
            "5/4,3/2",
        )
        self.assertIn("rational_factor=5/2", result.stdout)
        self.assertIn("two_power=0", result.stdout)
        self.assertIn("pi_power=0", result.stdout)
        self.assertIn(
            "exact_status=PROVED_BY_GAMMA_RECURRENCE_AND_SPECIAL_IDENTITIES",
            result.stdout,
        )

    def test_quarter_gamma_pair_is_closed_exactly(self) -> None:
        result = run_verifier(
            "gamma-quotient",
            "--numerator",
            "5/4,7/4",
            "--denominator",
            "1/2,1",
        )
        self.assertIn("rational_factor=3/16", result.stdout)
        self.assertIn("two_power=1/2", result.stdout)
        self.assertIn("pi_power=1/2", result.stdout)
        self.assertIn("unresolved=none", result.stdout)
        self.assertIn(
            "exact_status=PROVED_BY_GAMMA_RECURRENCE_AND_SPECIAL_IDENTITIES",
            result.stdout,
        )

    def test_perturbative_recurrence_extracts_second_order_forcing(self) -> None:
        result = run_verifier(
            "perturbative-recurrence",
            "--A=4*j**2 + 3*j + 5/16",
            "--B=4*j**2 - 3*j + 5/16",
            "--D=8*j",
            "--u=-4*n/3",
        )
        self.assertIn("mu1_residual=0", result.stdout)
        self.assertIn("mu1_status=PROVED_EXACT_SOLUTION", result.stdout)
        self.assertIn("mu2_forcing_rhs=32/3*j*n", result.stdout)
        self.assertIn(
            "exact_status=PROVED_BY_FORMAL_COEFFICIENT_EXTRACTION",
            result.stdout,
        )

    def test_bad_first_order_candidate_is_rejected(self) -> None:
        result = run_verifier(
            "perturbative-recurrence",
            "--A=4*j**2 + 3*j + 5/16",
            "--B=4*j**2 - 3*j + 5/16",
            "--D=8*j",
            "--u=4*n/3",
        )
        self.assertIn("mu1_status=CANDIDATE_U_FAILS", result.stdout)
        self.assertNotIn("mu1_residual=0\n", result.stdout)

    def test_affine_index_transform_proves_staggered_coefficients(self) -> None:
        result = run_verifier(
            "recurrence-transform",
            "--A=(d+1/2)*(d+3/2)",
            "--B=d*(d-1)",
            "--D=4*d+1",
            "--raw-var=d",
            "--site-var=j",
            "--raw-in-site=2*j-1/4",
            "--source-orientation=current_minus_neighbors_equals_mu",
        )
        self.assertIn("A_site=4*j^2 + 3*j + 5/16", result.stdout)
        self.assertIn("B_site=4*j^2 - 3*j + 5/16", result.stdout)
        self.assertIn("D_site=8*j", result.stdout)
        self.assertIn(
            "exact_status=PROVED_BY_AFFINE_SUBSTITUTION_AND_POLYNOMIAL_ARITHMETIC",
            result.stdout,
        )

    def test_asymptotic_power_prefers_half_power_for_blind_sequence(self) -> None:
        result = run_verifier(
            "asymptotic-power",
            "--n",
            "256,512,1024,2048,4096,8192,16384",
            "--y",
            "0.3060529601386244,0.3028866422527186,0.3006727188464453,"
            "0.2991193299086591,0.2980268264679914,0.2972572173673145,"
            "0.2967144592235840",
            "--powers",
            "1/2,1,3/2,2",
        )
        self.assertIn("best_power=1/2", result.stdout)
        self.assertIn("ASYMPTOTIC_FIT_IS_NUMERICAL_EVIDENCE_NOT_PROOF", result.stdout)

    def test_expression_injection_is_rejected(self) -> None:
        result = run_verifier(
            "rational",
            "--expr",
            "__import__('os').system('echo nope')",
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
