#!/usr/bin/env python3
"""Machine-readable evidence ledger for research-agent verification gates."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class EvidenceStatus(str, Enum):
    PROVED_EXACT = "proved_exact"
    ASYMPTOTIC_EVIDENCE = "asymptotic_evidence"
    NUMERICAL_EVIDENCE = "numerical_evidence"
    UNRESOLVED = "unresolved"
    REFUTED = "refuted"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class EvidenceRecord:
    mode: str
    status: EvidenceStatus
    fields: dict[str, str]
    raw_output: str


@dataclass
class ProofLedger:
    records: list[EvidenceRecord] = field(default_factory=list)

    def add_verifier_output(self, mode: str, output: str) -> EvidenceRecord:
        fields: dict[str, str] = {}
        for raw_line in output.splitlines():
            line = raw_line.strip()
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            fields[key.strip()] = value.strip()

        status = classify_verifier_output(fields, output)
        record = EvidenceRecord(
            mode=mode,
            status=status,
            fields=fields,
            raw_output=output,
        )
        self.records.append(record)
        return record

    def modes_used(self) -> set[str]:
        return {record.mode for record in self.records}

    def has_exact_success(self, mode: str) -> bool:
        return any(
            record.mode == mode and record.status == EvidenceStatus.PROVED_EXACT
            for record in self.records
        )

    def has_successful_perturbative_extraction(self) -> bool:
        return any(
            record.mode == "perturbative_recurrence"
            and record.fields.get("mu1_status") == "PROVED_EXACT_SOLUTION"
            and record.fields.get("exact_status")
            == "PROVED_BY_FORMAL_COEFFICIENT_EXTRACTION"
            for record in self.records
        )

    def has_exact_symbolic_mu2_chain(self) -> bool:
        """Return true only when both post-perturbative exact stages succeeded."""
        return self.has_exact_success(
            "symbolic_hypergeometric"
        ) and self.has_exact_success("symbolic_finite_part")

    def unresolved_gamma_seen(self) -> bool:
        return any(
            record.mode == "gamma_quotient"
            and record.status == EvidenceStatus.UNRESOLVED
            for record in self.records
        )

    def best_asymptotic_power(self) -> str | None:
        for record in reversed(self.records):
            if record.mode == "asymptotic_power":
                return record.fields.get("best_power")
        return None

    def gate_failures(
        self,
        *,
        required_modes: set[str] | None = None,
        require_exact_modes: set[str] | None = None,
        require_perturbative_success: bool = False,
        require_index_transform: bool = False,
        require_symbolic_mu2_chain: bool = False,
    ) -> list[str]:
        failures: list[str] = []
        required_modes = required_modes or set()
        require_exact_modes = require_exact_modes or set()

        missing = sorted(required_modes - self.modes_used())
        if missing:
            failures.append("missing verifier modes: " + ", ".join(missing))

        non_exact = sorted(
            mode for mode in require_exact_modes if not self.has_exact_success(mode)
        )
        if non_exact:
            failures.append(
                "verifier modes without an exact successful result: "
                + ", ".join(non_exact)
            )

        if require_perturbative_success and not self.has_successful_perturbative_extraction():
            failures.append(
                "no perturbative_recurrence call proved the first-order candidate "
                "and extracted the second-order equation"
            )

        if require_index_transform and not self.has_exact_success("recurrence_transform"):
            failures.append(
                "no recurrence_transform call exactly verified the index/sign normalization"
            )

        if require_symbolic_mu2_chain and not self.has_exact_symbolic_mu2_chain():
            failures.append(
                "the exact post-perturbative symbolic mu2 chain is incomplete: "
                "both symbolic_hypergeometric and symbolic_finite_part must succeed"
            )

        return failures

    def public_summary(self) -> str:
        lines = ["DETERMINISTIC EVIDENCE LEDGER"]
        if not self.records:
            lines.append("- no verifier evidence recorded")
            return "\n".join(lines)

        for index, record in enumerate(self.records, start=1):
            detail = ""
            if record.mode == "asymptotic_power" and record.fields.get("best_power"):
                detail = f" best_power={record.fields['best_power']}"
            elif record.mode == "perturbative_recurrence":
                detail = (
                    f" mu1_status={record.fields.get('mu1_status', 'unknown')}"
                    f" forcing={record.fields.get('mu2_forcing_rhs', 'unknown')}"
                )
            elif record.mode == "gamma_quotient":
                detail = f" exact_status={record.fields.get('exact_status', 'unknown')}"
            elif record.mode == "symbolic_hypergeometric":
                detail = (
                    f" candidate_status={record.fields.get('candidate_status', 'unknown')}"
                    f" exact_status={record.fields.get('exact_status', 'unknown')}"
                )
            elif record.mode == "symbolic_finite_part":
                detail = (
                    f" finite_part={record.fields.get('finite_part', 'unknown')}"
                    f" exact_status={record.fields.get('exact_status', 'unknown')}"
                )
            lines.append(
                f"- {index}: mode={record.mode} status={record.status.value}{detail}"
            )
        return "\n".join(lines)


def classify_verifier_output(fields: dict[str, str], output: str) -> EvidenceStatus:
    if fields.get("mu1_status") == "CANDIDATE_U_FAILS":
        return EvidenceStatus.REFUTED
    if fields.get("numeric_status") == "MISMATCH":
        return EvidenceStatus.REFUTED
    if fields.get("candidate_status") == "MISMATCH":
        return EvidenceStatus.REFUTED

    exact_status = fields.get("exact_status", "")
    if exact_status.startswith("REFUTED_"):
        return EvidenceStatus.REFUTED
    if exact_status.startswith("PROVED_"):
        return EvidenceStatus.PROVED_EXACT
    if exact_status == "UNRESOLVED_GAMMA_BASES":
        return EvidenceStatus.UNRESOLVED

    if fields.get("numeric_status") == "MATCH_WITHIN_TOLERANCE":
        return EvidenceStatus.NUMERICAL_EVIDENCE
    if "ASYMPTOTIC_FIT_IS_NUMERICAL_EVIDENCE_NOT_PROOF" in output:
        return EvidenceStatus.ASYMPTOTIC_EVIDENCE

    return EvidenceStatus.UNKNOWN
