#!/usr/bin/env python3
"""Regression tests for tools/symbolic_mu2.py."""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("symbolic_mu2.py")


def run_symbolic(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        text=True,
        capture_output=True,
        check=check,
        timeout=30,
    )


class SymbolicMu2Tests(unittest.TestCase):
    def test_plus_variation_forcing_ratio_is_derived_exactly(self) -> None:
        candidate = (
            "((4*n+1)*(4*n+3)*(8*n+9))"
            "/(8*n*(2*n+1)*(8*n+1))"
        )
        result = run_symbolic(
            "forcing-ratio",
            "--A",
            "4*j**2+3*j+5/16",
            "--B",
            "4*j**2-3*j+5/16",
            "--forcing",
            "(32/3)*j*n",
            "--offset",
            "1/8",
            "--candidate-ratio",
            candidate,
        )
        self.assertIn("candidate_difference=0", result.stdout)
        self.assertIn("candidate_status=PROVED_EQUAL", result.stdout)
        self.assertIn(
            "exact_status=PROVED_BY_VARIATION_OF_CONSTANTS_QUOTIENT",
            result.stdout,
        )

    def test_minus_variation_forcing_ratio_is_derived_exactly(self) -> None:
        candidate = (
            "((4*n+3)*(4*n+5)*(8*n+13))"
            "/(8*n*(2*n+3)*(8*n+5))"
        )
        result = run_symbolic(
            "forcing-ratio",
            "--A",
            "4*j**2+3*j+5/16",
            "--B",
            "4*j**2-3*j+5/16",
            "--forcing",
            "(32/3)*j*n",
            "--offset",
            "5/8",
            "--candidate-ratio",
            candidate,
        )
        self.assertIn("candidate_difference=0", result.stdout)
        self.assertIn("candidate_status=PROVED_EQUAL", result.stdout)
        self.assertIn(
            "exact_status=PROVED_BY_VARIATION_OF_CONSTANTS_QUOTIENT",
            result.stdout,
        )

    def test_wrong_variation_forcing_ratio_is_refuted(self) -> None:
        result = run_symbolic(
            "forcing-ratio",
            "--A",
            "4*j**2+3*j+5/16",
            "--B",
            "4*j**2-3*j+5/16",
            "--forcing",
            "(32/3)*j*n",
            "--offset",
            "1/8",
            "--candidate-ratio",
            "(n+1)/(n+2)",
        )
        self.assertIn("candidate_status=MISMATCH", result.stdout)
        self.assertIn("exact_status=REFUTED_FORCING_QUOTIENT", result.stdout)
        self.assertNotIn(
            "exact_status=PROVED_BY_VARIATION_OF_CONSTANTS_QUOTIENT",
            result.stdout,
        )

    def test_hypergeometric_ratio_and_terms_are_exact(self) -> None:
        result = run_symbolic(
            "hypergeometric",
            "--numerator-shifts",
            "1/2",
            "--denominator-shifts",
            "1",
            "--base",
            "1",
            "--candidate-ratio",
            "(k+1/2)/(k+1)",
            "--terms",
            "5",
        )
        self.assertIn("first_terms=1,1/2,3/8,5/16,35/128", result.stdout)
        self.assertIn("candidate_difference=0", result.stdout)
        self.assertIn("candidate_status=PROVED_EQUAL", result.stdout)
        self.assertIn("exact_status=PROVED_BY_POCHHAMMER_QUOTIENT", result.stdout)

    def test_hypergeometric_wrong_ratio_is_refuted_fail_closed(self) -> None:
        result = run_symbolic(
            "hypergeometric",
            "--numerator-shifts",
            "1/2",
            "--denominator-shifts",
            "1",
            "--candidate-ratio",
            "(k+3/2)/(k+1)",
        )
        self.assertIn("candidate_status=MISMATCH", result.stdout)
        self.assertNotIn("candidate_difference=0\n", result.stdout)
        self.assertIn("exact_status=REFUTED_CANDIDATE_RATIO", result.stdout)
        self.assertNotIn("exact_status=PROVED_BY_POCHHAMMER_QUOTIENT", result.stdout)

    def test_hypergeometric_series_is_recognized(self) -> None:
        result = run_symbolic(
            "hypergeometric",
            "--numerator-shifts",
            "1/3,2/3",
            "--denominator-shifts",
            "1,1",
            "--base",
            "1/4",
        )
        self.assertIn("hypergeometric_series=_2F_1", result.stdout)
        self.assertIn("exact_status=PROVED_BY_POCHHAMMER_QUOTIENT", result.stdout)

    def test_mu2_forcing_subsequence_a_k_is_exact_pochhammer_sequence(self) -> None:
        result = run_symbolic(
            "hypergeometric",
            "--numerator-shifts",
            "1/4,3/4",
            "--denominator-shifts",
            "1/2,1",
            "--base",
            "1",
            "--candidate-ratio",
            "((k+1/4)*(k+3/4))/((k+1/2)*(k+1))",
            "--terms",
            "6",
        )
        self.assertIn(
            "first_terms=1,3/8,35/128,231/1024,6435/32768,46189/262144",
            result.stdout,
        )
        self.assertIn("candidate_difference=0", result.stdout)
        self.assertIn("hypergeometric_series=_2F_1([1/4,3/4];[1/2];1*z)", result.stdout)
        self.assertIn("exact_status=PROVED_BY_POCHHAMMER_QUOTIENT", result.stdout)

    def test_simple_puiseux_finite_part(self) -> None:
        result = run_symbolic(
            "finite-part",
            "--base-expr",
            "1/sqrt(1-z)",
            "--theta-polynomial",
            "1",
            "--extra-expr",
            "0",
            "--order",
            "8",
        )
        self.assertIn("leading_power=-1", result.stdout)
        self.assertIn("finite_part=0", result.stdout)
        self.assertIn(
            "exact_status=PROVED_BY_EXACT_THETA_ALGEBRA_AND_PUISEUX_SERIES",
            result.stdout,
        )

    def test_theta_operator_and_constant_extra_term(self) -> None:
        result = run_symbolic(
            "finite-part",
            "--base-expr",
            "1/(1-z)",
            "--theta-polynomial",
            "k",
            "--extra-expr",
            "3/5",
            "--order",
            "8",
        )
        self.assertIn("leading_power=-4", result.stdout)
        self.assertIn("finite_part=3/5", result.stdout)

    def test_nested_radicals_are_supported(self) -> None:
        result = run_symbolic(
            "finite-part",
            "--base-expr",
            "sqrt(1+sqrt(z))",
            "--theta-polynomial",
            "1",
            "--extra-expr",
            "0",
            "--order",
            "8",
        )
        self.assertIn("finite_part=sqrt(2)", result.stdout)
        self.assertIn("leading_power=0", result.stdout)

    def test_mu2_plus_generating_finite_part_is_exact(self) -> None:
        base = "(1/sqrt(1-sqrt(z))+1/sqrt(1+sqrt(z)))/2"
        result = run_symbolic(
            "finite-part",
            "--base-expr",
            base,
            "--theta-polynomial",
            "(16/9)*(8*k**2+k)",
            "--extra-expr",
            "0",
            "--order",
            "10",
        )
        self.assertIn("leading_power=-5", result.stdout)
        self.assertIn("finite_part=-sqrt(2)/9", result.stdout)
        self.assertIn(
            "exact_status=PROVED_BY_EXACT_THETA_ALGEBRA_AND_PUISEUX_SERIES",
            result.stdout,
        )

    def test_mu2_minus_generating_finite_part_is_exact(self) -> None:
        base = "(1/sqrt(1-sqrt(z))+1/sqrt(1+sqrt(z)))/2"
        extra = "(8/45)*(sqrt(1+sqrt(z))-sqrt(1-sqrt(z)))/sqrt(z)"
        result = run_symbolic(
            "finite-part",
            "--base-expr",
            base,
            "--theta-polynomial",
            "(16/45)*(16*k**2+6*k-1/2)",
            "--extra-expr",
            extra,
            "--order",
            "10",
        )
        self.assertIn("leading_power=-5", result.stdout)
        self.assertIn("finite_part=2*sqrt(2)/45", result.stdout)
        self.assertIn(
            "exact_status=PROVED_BY_EXACT_THETA_ALGEBRA_AND_PUISEUX_SERIES",
            result.stdout,
        )

    def test_code_injection_is_rejected(self) -> None:
        result = run_symbolic(
            "finite-part",
            "--base-expr",
            "__import__('os').system('echo nope')",
            "--theta-polynomial",
            "1",
            "--extra-expr",
            "0",
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
