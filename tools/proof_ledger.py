#!/usr/bin/env python3
"""Machine-readable evidence ledger for research-agent verification gates."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from enum import Enum


class EvidenceStatus(str, Enum):
    PROVED_EXACT = "proved_exact"
    ASYMPTOTIC_EVIDENCE = "asymptotic_evidence"
    NUMERICAL_EVIDENCE = "numerical_evidence"
    UNRESOLVED = "unresolved"
    REFUTED = "refuted"
    UNKNOWN = "unknown"


class EvidenceAuthority(str, Enum):
    SCRIPT_VERIFIED_EXACT = "script_verified_exact"
    SCRIPT_NUMERICAL = "script_numerical"
    SCRIPT_ASYMPTOTIC = "script_asymptotic"
    SCRIPT_REFUTED = "script_refuted"
    SCRIPT_UNRESOLVED = "script_unresolved"
    INVALID_OR_UNBOUND = "invalid_or_unbound"


EXPECTED_OUTPUT_MODES = {
    "symbolic_forcing_ratio": "forcing_ratio",
    "symbolic_hypergeometric": "hypergeometric",
    "symbolic_finite_part": "finite_part",
    "symbolic_assembly": "assembly",
}

REQUIRED_FIELDS = {
    "rational": {"mode", "exact_status"},
    "gamma_quotient": {"mode", "exact_status"},
    "numeric_identity": {"mode", "numeric_status", "warning"},
    "asymptotic_power": {"mode", "best_power", "warning"},
    "perturbative_recurrence": {"mode", "mu1_status", "exact_status"},
    "recurrence_transform": {"mode", "exact_status"},
    "symbolic_forcing_ratio": {"mode", "candidate_status", "exact_status"},
    "symbolic_hypergeometric": {"mode", "candidate_status", "exact_status"},
    "symbolic_finite_part": {"mode", "finite_part", "exact_status"},
    "symbolic_assembly": {"mode", "candidate_status", "exact_status"},
}


def _canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class EvidenceRecord:
    mode: str
    status: EvidenceStatus
    authority: EvidenceAuthority
    proposition: str
    proposition_sha256: str
    input_sha256: str
    verifier_id: str
    verifier_sha256: str
    source_sha: str
    execution_id: str
    exit_code: int | None
    schema_valid: bool
    identity_complete: bool
    validation_errors: tuple[str, ...]
    prooflab_receipt_id: str | None
    gate_target: bool
    fields: dict[str, str]
    raw_output: str


@dataclass
class ProofLedger:
    records: list[EvidenceRecord] = field(default_factory=list)
    source_sha: str | None = None
    _allow_test_defaults: bool = False
    _test_counter: int = 0
    gate_targets: dict[str, tuple[str, str]] = field(default_factory=dict)

    @classmethod
    def testing(cls) -> "ProofLedger":
        """Create an explicitly synthetic ledger for unit tests only."""
        return cls(source_sha="0" * 40, _allow_test_defaults=True)

    def add_verifier_output(
        self,
        mode: str,
        output: str,
        *,
        arguments: dict | None = None,
        proposition: str | None = None,
        execution_id: str | None = None,
        verifier_id: str | None = None,
        verifier_sha256: str | None = None,
    ) -> EvidenceRecord:
        if self._allow_test_defaults:
            self._test_counter += 1
            arguments = arguments or {"mode": mode, "gate_target": True}
            proposition = proposition or f"unit-test proposition for {mode}"
            execution_id = execution_id or f"unit-test-execution-{self._test_counter}"
            verifier_id = verifier_id or "unit-test-verifier"
            verifier_sha256 = verifier_sha256 or "1" * 64
            if not any(line.startswith("exit_status=") for line in output.splitlines()):
                output = "exit_status=0\n" + output
            if not any(line.startswith("mode=") for line in output.splitlines()):
                output = (
                    f"mode={EXPECTED_OUTPUT_MODES.get(mode, mode)}\n" + output
                )

        arguments = arguments if isinstance(arguments, dict) else {}
        fields: dict[str, str] = {}
        validation_errors: list[str] = []
        for raw_line in output.splitlines():
            line = raw_line.strip()
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            if key in fields:
                validation_errors.append(f"duplicate output field: {key}")
            fields[key] = value.strip()

        exit_code: int | None = None
        try:
            exit_code = int(fields.get("exit_status", ""))
        except ValueError:
            validation_errors.append("exit_status is missing or not an integer")

        required = REQUIRED_FIELDS.get(mode)
        if required is None:
            validation_errors.append(f"unsupported evidence mode: {mode}")
        else:
            missing = sorted(required - fields.keys())
            if missing:
                validation_errors.append("missing output fields: " + ", ".join(missing))
            expected_mode = EXPECTED_OUTPUT_MODES.get(mode, mode)
            if fields.get("mode") != expected_mode:
                validation_errors.append(
                    f"output mode mismatch: expected {expected_mode!r}, "
                    f"got {fields.get('mode')!r}"
                )
        if exit_code != 0:
            validation_errors.append(f"verifier exit code is not zero: {exit_code!r}")

        identity_complete = all(
            (
                isinstance(self.source_sha, str) and len(self.source_sha) == 40,
                isinstance(proposition, str) and bool(proposition.strip()),
                isinstance(execution_id, str) and bool(execution_id.strip()),
                isinstance(verifier_id, str) and bool(verifier_id.strip()),
                isinstance(verifier_sha256, str) and len(verifier_sha256) == 64,
                bool(arguments),
            )
        )
        if not identity_complete:
            validation_errors.append("evidence identity is incomplete")

        proposition_text = proposition.strip() if isinstance(proposition, str) else ""
        input_arguments = {
            key: value
            for key, value in arguments.items()
            if key not in {"proposition", "gate_target"}
        }
        proposition_sha256 = _canonical_sha256(proposition_text)
        input_sha256 = _canonical_sha256(input_arguments)
        gate_target = arguments.get("gate_target") is True
        if gate_target and identity_complete:
            identity = (proposition_sha256, input_sha256)
            existing = self.gate_targets.get(mode)
            if existing is None:
                self.gate_targets[mode] = identity
            elif existing != identity:
                validation_errors.append(
                    "gate target identity is already locked for this mode"
                )

        schema_valid = not validation_errors
        status = (
            classify_verifier_output(fields, output)
            if schema_valid and identity_complete
            else EvidenceStatus.UNKNOWN
        )
        authority = evidence_authority(status, schema_valid and identity_complete)

        record = EvidenceRecord(
            mode=mode,
            status=status,
            authority=authority,
            proposition=proposition_text,
            proposition_sha256=proposition_sha256,
            input_sha256=input_sha256,
            verifier_id=verifier_id or "",
            verifier_sha256=verifier_sha256 or "",
            source_sha=self.source_sha or "",
            execution_id=execution_id or "",
            exit_code=exit_code,
            schema_valid=schema_valid,
            identity_complete=identity_complete,
            validation_errors=tuple(validation_errors),
            prooflab_receipt_id=None,
            gate_target=gate_target,
            fields=fields,
            raw_output=output,
        )
        self.records.append(record)
        return record

    def modes_used(self) -> set[str]:
        return {
            mode
            for mode in self.gate_targets
            if (record := self.gate_record(mode))
            and record.schema_valid
            and record.identity_complete
        }

    def latest_record(self, mode: str) -> EvidenceRecord | None:
        return next(
            (record for record in reversed(self.records) if record.mode == mode),
            None,
        )

    def gate_record(self, mode: str) -> EvidenceRecord | None:
        identity = self.gate_targets.get(mode)
        if identity is None:
            return None
        proposition_sha256, input_sha256 = identity
        return next(
            (
                record
                for record in reversed(self.records)
                if record.mode == mode
                and record.gate_target
                and record.proposition_sha256 == proposition_sha256
                and record.input_sha256 == input_sha256
            ),
            None,
        )

    def has_exact_success(
        self,
        mode: str,
        *,
        proposition_sha256: str | None = None,
        input_sha256: str | None = None,
    ) -> bool:
        if proposition_sha256 is None and input_sha256 is None:
            record = self.gate_record(mode)
        else:
            record = next(
                (
                    candidate
                    for candidate in reversed(self.records)
                    if candidate.mode == mode
                    and (
                        proposition_sha256 is None
                        or candidate.proposition_sha256 == proposition_sha256
                    )
                    and (
                        input_sha256 is None
                        or candidate.input_sha256 == input_sha256
                    )
                ),
                None,
            )
        return bool(
            record
            and record.status == EvidenceStatus.PROVED_EXACT
            and record.authority == EvidenceAuthority.SCRIPT_VERIFIED_EXACT
            and (
                proposition_sha256 is None
                or record.proposition_sha256 == proposition_sha256
            )
            and (input_sha256 is None or record.input_sha256 == input_sha256)
        )

    def has_prooflab_acceptance(self, mode: str) -> bool:
        record = self.gate_record(mode)
        return bool(record and record.prooflab_receipt_id)

    def has_successful_perturbative_extraction(self) -> bool:
        record = self.gate_record("perturbative_recurrence")
        return bool(
            record
            and record.fields.get("mu1_status") == "PROVED_EXACT_SOLUTION"
            and record.fields.get("exact_status")
            == "PROVED_BY_FORMAL_COEFFICIENT_EXTRACTION"
            and record.status == EvidenceStatus.PROVED_EXACT
        )

    def has_exact_symbolic_mu2_chain(self) -> bool:
        """Return true only when all post-perturbative exact stages succeeded."""
        return (
            self.has_exact_success("symbolic_forcing_ratio")
            and self.has_exact_success("symbolic_hypergeometric")
            and self.has_exact_success("symbolic_finite_part")
            and self.has_exact_success("symbolic_assembly")
        )

    def unresolved_gamma_seen(self) -> bool:
        return any(
            record.mode == "gamma_quotient"
            and record.status == EvidenceStatus.UNRESOLVED
            for record in self.records
        )

    def best_asymptotic_power(self) -> str | None:
        record = self.gate_record("asymptotic_power")
        if not record or not record.schema_valid or not record.identity_complete:
            return None
        return record.fields.get("best_power")

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
                "symbolic_forcing_ratio, symbolic_hypergeometric, "
                "symbolic_finite_part and symbolic_assembly must all succeed"
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
            elif record.mode == "symbolic_forcing_ratio":
                detail = (
                    f" candidate_status={record.fields.get('candidate_status', 'unknown')}"
                    f" ratio={record.fields.get('derived_ratio', 'unknown')}"
                    f" exact_status={record.fields.get('exact_status', 'unknown')}"
                )
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
            elif record.mode == "symbolic_assembly":
                detail = (
                    f" assembled_value={record.fields.get('assembled_value', 'unknown')}"
                    f" candidate_status={record.fields.get('candidate_status', 'unknown')}"
                    f" exact_status={record.fields.get('exact_status', 'unknown')}"
                )
            lines.append(
                f"- {index}: mode={record.mode} status={record.status.value} "
                f"authority={record.authority.value} "
                f"proposition_sha256={record.proposition_sha256} "
                f"input_sha256={record.input_sha256} "
                f"verifier={record.verifier_id}@sha256:{record.verifier_sha256} "
                f"source_sha={record.source_sha} execution_id={record.execution_id} "
                f"exit_code={record.exit_code} schema_valid={record.schema_valid} "
                f"gate_target={record.gate_target} prooflab=not_submitted{detail}"
            )
            if record.validation_errors:
                lines.append("  validation_errors=" + "; ".join(record.validation_errors))
        return "\n".join(lines)


def evidence_authority(
    status: EvidenceStatus,
    is_valid: bool,
) -> EvidenceAuthority:
    if not is_valid:
        return EvidenceAuthority.INVALID_OR_UNBOUND
    return {
        EvidenceStatus.PROVED_EXACT: EvidenceAuthority.SCRIPT_VERIFIED_EXACT,
        EvidenceStatus.NUMERICAL_EVIDENCE: EvidenceAuthority.SCRIPT_NUMERICAL,
        EvidenceStatus.ASYMPTOTIC_EVIDENCE: EvidenceAuthority.SCRIPT_ASYMPTOTIC,
        EvidenceStatus.REFUTED: EvidenceAuthority.SCRIPT_REFUTED,
        EvidenceStatus.UNRESOLVED: EvidenceAuthority.SCRIPT_UNRESOLVED,
        EvidenceStatus.UNKNOWN: EvidenceAuthority.INVALID_OR_UNBOUND,
    }[status]


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
