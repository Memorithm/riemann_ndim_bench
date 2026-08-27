#!/usr/bin/env python3
"""Adversarial Qwen/Nemotron research harness with deterministic final gates."""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path

import riemann_research_agent as base
from proof_ledger import ProofLedger


RESEARCHER_SYSTEM = r"""
You are the primary mathematical researcher in an adversarial blind experiment.
Use only the supplied read-only tools. Derive claims from the workspace rather
than guessing a hidden target constant.

For a mu^2 recurrence problem:
- locate the raw recurrence and every index convention;
- use verify_math recurrence_transform to prove any affine raw-index to
  staggered/site-index conversion and sign normalization;
- derive the first-order coefficient from the recurrence/boundary condition;
- use perturbative_recurrence and treat CANDIDATE_U_FAILS as a hard refutation;
- use asymptotic_power before choosing a finite-size exponent;
- use gamma_quotient for Gamma recurrence claims;
- classify exact, asymptotic, numerical and unresolved statements separately.
Never claim a proof of the Riemann hypothesis.
""".strip()


CRITIC_SYSTEM = r"""
You are an adversarial mathematical referee. Independently inspect the workspace
and try to refute the researcher's public proposal.

Audit index conversions, signs, normalization, Gamma arithmetic and finite-size
powers with verify_math rather than trusting prose. For a mu^2 claim, reject any
argument that substitutes the homogeneous zero mode for the inhomogeneous
second-order equation. A failed first-order residual is a hard contradiction.
Never claim a proof of the Riemann hypothesis.
""".strip()


FINAL_SYSTEM = r"""
You are producing the final public synthesis of an adversarial mathematical
experiment. Preserve unresolved gaps. The deterministic evidence ledger is
binding: numerical/asymptotic evidence cannot be promoted to symbolic proof.
For the blind mu^2 task, the final phase must independently execute the required
verifier modes and pass the exact index-transform, Gamma and perturbative gates
before its report can be accepted.
""".strip()


FINAL_REQUIRED_MODES = {
    "recurrence_transform",
    "perturbative_recurrence",
    "asymptotic_power",
    "gamma_quotient",
}

FINAL_REQUIRED_EXACT_MODES = {
    "recurrence_transform",
    "gamma_quotient",
}


def append_jsonl(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")


def final_gate_failures(ledger: ProofLedger) -> list[str]:
    return ledger.gate_failures(
        required_modes=FINAL_REQUIRED_MODES,
        require_exact_modes=FINAL_REQUIRED_EXACT_MODES,
        require_perturbative_success=True,
        require_index_transform=True,
    )


def run_phase(
    *,
    root: Path,
    model: str,
    system_prompt: str,
    assignment: str,
    phase: str,
    transcript: Path,
    max_tool_turns: int,
    enforce_final_gate: bool = False,
) -> tuple[str, ProofLedger]:
    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": assignment},
    ]
    ledger = ProofLedger()
    final_audit_sent = False
    turn_limit = max_tool_turns + 5 if enforce_final_gate else max_tool_turns

    print("\n" + "=" * 78)
    print(f"PHASE: {phase}")
    print(f"MODEL: {model}")
    print("=" * 78)
    append_jsonl(
        transcript,
        {
            "event": "phase_start",
            "phase": phase,
            "model": model,
            "assignment": assignment,
            "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        },
    )

    for turn in range(1, turn_limit + 1):
        print(f"\n--- {phase}: model turn {turn}/{turn_limit} ---")
        try:
            result = base.call_ollama(model, messages)
        except Exception as exc:
            append_jsonl(
                transcript,
                {
                    "event": "ollama_error",
                    "phase": phase,
                    "turn": turn,
                    "model": model,
                    "error": str(exc),
                },
            )
            return f"OLLAMA_ERROR: {exc}", ledger

        message = result.get("message", {})
        content = message.get("content", "") or ""
        calls = message.get("tool_calls", []) or []
        if content:
            print(content)
        messages.append(base.public_assistant_message(message))
        append_jsonl(
            transcript,
            {
                "event": "assistant_public",
                "phase": phase,
                "turn": turn,
                "model": model,
                "content": content,
                "tool_calls": calls,
                "prompt_eval_count": result.get("prompt_eval_count"),
                "eval_count": result.get("eval_count"),
            },
        )

        if not calls:
            if enforce_final_gate:
                failures = final_gate_failures(ledger)
                if failures:
                    gate = (
                        "VERIFICATION GATE REJECTED THIS DRAFT. Resolve all of "
                        "the following with actual verify_math calls before "
                        "trying to finish:\n- "
                        + "\n- ".join(failures)
                        + "\nThe deterministic ledger is binding."
                    )
                    print("\n--- VERIFICATION GATE ---")
                    print(gate)
                    append_jsonl(
                        transcript,
                        {
                            "event": "verification_gate_reject",
                            "phase": phase,
                            "turn": turn,
                            "model": model,
                            "failures": failures,
                            "draft": content,
                        },
                    )
                    messages.append({"role": "user", "content": gate})
                    continue

                if not final_audit_sent:
                    final_audit_sent = True
                    audit = (
                        "FINAL EVIDENCE ATTESTATION. The gate now passes. Before "
                        "the report can be accepted, revise it once against this "
                        "machine ledger. Do not state a different finite-size "
                        "best power than the ledger for the same data, do not "
                        "describe numerical/asymptotic evidence as proof, and do "
                        "not promote a Gamma identity that lacks an exact success.\n\n"
                        + ledger.public_summary()
                    )
                    print("\n--- FINAL EVIDENCE ATTESTATION ---")
                    print(audit)
                    append_jsonl(
                        transcript,
                        {
                            "event": "final_evidence_attestation",
                            "phase": phase,
                            "turn": turn,
                            "model": model,
                            "draft": content,
                            "ledger": ledger.public_summary(),
                        },
                    )
                    messages.append({"role": "user", "content": audit})
                    continue

            print(f"\n--- PUBLIC REPORT: {phase} ---")
            print(content)
            print(ledger.public_summary())
            return content, ledger

        for call in calls:
            name, output = base.execute_tool(root, call)
            print(f"\n--- tool: {name} ---")
            print(output)
            messages.append({"role": "tool", "tool_name": name, "content": output})
            arguments = call.get("function", {}).get("arguments", {}) or {}
            append_jsonl(
                transcript,
                {
                    "event": "tool",
                    "phase": phase,
                    "turn": turn,
                    "model": model,
                    "tool": name,
                    "arguments": arguments,
                    "output": output,
                },
            )
            if name == "verify_math" and isinstance(arguments, dict):
                mode = arguments.get("mode", "")
                if mode:
                    ledger.add_verifier_output(mode, output)

    if enforce_final_gate:
        failures = final_gate_failures(ledger)
        if failures:
            report = "VERIFICATION GATE FAILED: final synthesis withheld; " + "; ".join(failures)
            print("\n--- VERIFICATION GATE FAILURE ---")
            print(report)
            append_jsonl(
                transcript,
                {
                    "event": "verification_gate_failure",
                    "phase": phase,
                    "model": model,
                    "failures": failures,
                    "content": report,
                },
            )
            return report, ledger

    messages.append(
        {
            "role": "user",
            "content": (
                "Tool budget is exhausted. Do not call another tool. Give the "
                "public scientific report now, preserving the deterministic "
                "evidence classifications and all unresolved gaps.\n\n"
                + ledger.public_summary()
            ),
        }
    )
    try:
        result = base.call_ollama(model, messages, tools=False)
        content = result.get("message", {}).get("content", "") or ""
    except Exception as exc:
        content = f"OLLAMA_ERROR_DURING_FORCED_REPORT: {exc}"

    append_jsonl(
        transcript,
        {
            "event": "forced_public_report",
            "phase": phase,
            "model": model,
            "content": content,
        },
    )
    print(f"\n--- FORCED PUBLIC REPORT: {phase} ---")
    print(content)
    print(ledger.public_summary())
    return content, ledger


def read_challenge(root: Path, challenge_path: str) -> str:
    path = base.make_safe_path(root, challenge_path)
    if not path.is_file():
        raise RuntimeError(f"challenge file does not exist: {challenge_path}")
    return path.read_text(encoding="utf-8")


def collaborative_run(
    *,
    root: Path,
    researcher: str,
    critic: str,
    rounds: int,
    max_tool_turns: int,
    transcript: Path,
    challenge_path: str,
) -> str:
    challenge = read_challenge(root, challenge_path)
    researcher_reports: list[str] = []
    critic_reports: list[str] = []
    prior_critique = ""

    append_jsonl(
        transcript,
        {
            "event": "experiment_start",
            "root": str(root),
            "researcher": researcher,
            "critic": critic,
            "rounds": rounds,
            "challenge_path": challenge_path,
            "challenge": challenge,
            "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        },
    )

    for round_no in range(1, rounds + 1):
        if round_no == 1:
            assignment = (
                "BLIND RESEARCH TASK\n\n"
                + challenge
                + "\n\nThis is round 1. Investigate independently and produce a public derivation."
            )
        else:
            assignment = (
                "BLIND RESEARCH TASK\n\n"
                + challenge
                + "\n\nPREVIOUS PUBLIC CRITIQUE\n\n"
                + prior_critique
                + f"\n\nThis is research round {round_no}. Correct actual errors and deepen the derivation."
            )

        q_report, _ = run_phase(
            root=root,
            model=researcher,
            system_prompt=RESEARCHER_SYSTEM,
            assignment=assignment,
            phase=f"research-round-{round_no}",
            transcript=transcript,
            max_tool_turns=max_tool_turns,
        )
        researcher_reports.append(q_report)

        critic_assignment = (
            "BLIND RESEARCH TASK\n\n"
            + challenge
            + f"\n\nPUBLIC PROPOSAL ROUND {round_no}\n\n"
            + q_report
            + "\n\nIndependently inspect the workspace and try to refute this proposal."
        )
        n_report, _ = run_phase(
            root=root,
            model=critic,
            system_prompt=CRITIC_SYSTEM,
            assignment=critic_assignment,
            phase=f"critique-round-{round_no}",
            transcript=transcript,
            max_tool_turns=max_tool_turns,
        )
        critic_reports.append(n_report)
        prior_critique = n_report

    dossier_parts = ["BLIND RESEARCH TASK", challenge]
    for index, (research, critique) in enumerate(zip(researcher_reports, critic_reports), start=1):
        dossier_parts.extend(
            [
                f"RESEARCH REPORT ROUND {index}",
                research,
                f"CRITIQUE ROUND {index}",
                critique,
            ]
        )
    dossier = "\n\n".join(dossier_parts)
    final_assignment = (
        "ADVERSARIAL PUBLIC DOSSIER\n\n"
        + dossier
        + "\n\nProduce the final synthesis. Reinspect the workspace and independently run "
        "all deterministic checks required by the final gate. Preserve every unresolved gap."
    )
    final_report, final_ledger = run_phase(
        root=root,
        model=researcher,
        system_prompt=FINAL_SYSTEM,
        assignment=final_assignment,
        phase="final-synthesis",
        transcript=transcript,
        max_tool_turns=max_tool_turns,
        enforce_final_gate=True,
    )
    append_jsonl(
        transcript,
        {
            "event": "experiment_final",
            "final_report": final_report,
            "ledger": final_ledger.public_summary(),
        },
    )
    print("\n" + "#" * 78)
    print("FINAL COLLABORATIVE REPORT")
    print("#" * 78)
    print(final_report)
    print(f"TRANSCRIPT={transcript}")
    return final_report


def main() -> None:
    parser = argparse.ArgumentParser(description="Adversarial dual-agent Riemann research harness")
    parser.add_argument("--root", default=".")
    parser.add_argument("--researcher", default="qwen3.8:latest")
    parser.add_argument("--critic", default="nemotron-3.5-lightning:30b")
    parser.add_argument("--rounds", type=int, default=4)
    parser.add_argument("--max-tool-turns", type=int, default=10)
    parser.add_argument("--challenge", default="BLIND_MU2_CHALLENGE.md")
    parser.add_argument("--transcript")
    args = parser.parse_args()

    if not 1 <= args.rounds <= 12:
        raise RuntimeError("--rounds must be between 1 and 12")
    if not 1 <= args.max_tool_turns <= 30:
        raise RuntimeError("--max-tool-turns must be between 1 and 30")
    root = Path(args.root).resolve()
    if not root.is_dir():
        raise RuntimeError(f"workspace root does not exist: {root}")

    if args.transcript:
        transcript = Path(args.transcript).resolve()
    else:
        stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
        transcript = root / "agent_runs" / f"dual-research-{stamp}.jsonl"

    collaborative_run(
        root=root,
        researcher=args.researcher,
        critic=args.critic,
        rounds=args.rounds,
        max_tool_turns=args.max_tool_turns,
        transcript=transcript,
        challenge_path=args.challenge,
    )


if __name__ == "__main__":
    main()
