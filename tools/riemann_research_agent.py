#!/usr/bin/env python3
"""Read-only multi-turn mathematical research agent.

The runner exposes a deliberately small tool surface, records public model
messages/tool outputs as JSONL, can resume a previous transcript, and maintains
a deterministic ProofLedger for verify_math calls.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

from proof_ledger import ProofLedger


OLLAMA_URL = os.environ.get("OLLAMA_URL", "http://127.0.0.1:11434/api/chat")
OLLAMA_TIMEOUT = int(os.environ.get("RIEMANN_AGENT_TIMEOUT", "1800"))
OLLAMA_RETRIES = int(os.environ.get("RIEMANN_AGENT_RETRIES", "2"))
MAX_TOOL_OUTPUT = 40_000
MAX_READ_LINES = 200


def fail(message: str) -> None:
    raise RuntimeError(message)


def limited(text: str, limit: int = MAX_TOOL_OUTPUT) -> str:
    if len(text) <= limit:
        return text
    return text[:limit] + f"\n...[truncated {len(text) - limit} chars]"


def make_safe_path(root: Path, relative: str) -> Path:
    if not isinstance(relative, str) or not relative.strip():
        fail("path must be a non-empty relative string")
    candidate = Path(relative)
    if candidate.is_absolute():
        fail("absolute paths are refused")
    resolved = (root / candidate).resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError as exc:
        raise RuntimeError("path escapes research workspace") from exc
    return resolved


def run_process(argv: list[str], cwd: Path, timeout: int = 120) -> str:
    result = subprocess.run(
        argv,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        check=False,
    )
    return f"exit_status={result.returncode}\n{result.stdout}".rstrip()


def tool_read(root: Path, args: dict) -> str:
    path = make_safe_path(root, args.get("path", ""))
    if not path.is_file():
        fail("read requires a regular file")
    start = args.get("start", 1)
    count = args.get("count", 120)
    if not isinstance(start, int) or start < 1:
        fail("read start must be a positive integer")
    if not isinstance(count, int) or count < 1 or count > MAX_READ_LINES:
        fail(f"read count must be between 1 and {MAX_READ_LINES}")
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    end = min(len(lines), start - 1 + count)
    return "\n".join(f"{index + 1:6d}: {lines[index]}" for index in range(start - 1, end))


def tool_search(root: Path, args: dict) -> str:
    query = args.get("query", "")
    relative = args.get("path", ".")
    max_results = args.get("max_results", 80)
    if not isinstance(query, str) or not query:
        fail("search query must be a non-empty string")
    if len(query) > 512:
        fail("search query too long")
    if not isinstance(max_results, int) or not 1 <= max_results <= 200:
        fail("max_results must be between 1 and 200")
    base = make_safe_path(root, relative)
    paths = [base] if base.is_file() else sorted(base.rglob("*"))
    matches: list[str] = []
    for path in paths:
        if len(matches) >= max_results:
            break
        if not path.is_file():
            continue
        if any(part in {".git", "target", "agent_runs"} for part in path.parts):
            continue
        try:
            lines = path.read_text(encoding="utf-8").splitlines()
        except (UnicodeDecodeError, OSError):
            continue
        for line_no, line in enumerate(lines, start=1):
            if query in line:
                rel = path.relative_to(root)
                matches.append(f"{rel}:{line_no}:{line}")
                if len(matches) >= max_results:
                    break
    return "\n".join(matches) if matches else "NO_MATCHES"


def tool_status(root: Path, _args: dict) -> str:
    branch = run_process(["git", "branch", "--show-current"], root, timeout=20)
    status = run_process(["git", "status", "--short"], root, timeout=20)
    stat = run_process(["git", "diff", "--stat"], root, timeout=20)
    return f"=== BRANCH ===\n{branch}\n=== STATUS ===\n{status}\n=== DIFF STAT ===\n{stat}"


def tool_test(root: Path, args: dict) -> str:
    target = args.get("target", "")
    commands = {
        "cargo_check": ["cargo", "check", "--all-targets"],
        "cargo_test": ["cargo", "test", "--all-targets"],
        "cargo_fmt_check": ["cargo", "fmt", "--all", "--", "--check"],
        "verifier_tests": [sys.executable, "-m", "unittest", "-v", "tools/test_verify_math.py", "tools/test_proof_ledger.py"],
    }
    argv = commands.get(target)
    if argv is None:
        fail("unsupported test target")
    return run_process(argv, root, timeout=900)


def text_arg(args: dict, name: str, *, required: bool = True, default=None, max_len: int = 4096):
    value = args.get(name, default)
    if value is None:
        if required:
            fail(f"verify_math argument is required: {name}")
        return None
    if not isinstance(value, str):
        fail(f"verify_math argument must be a string: {name}")
    if required and not value.strip():
        fail(f"verify_math argument is empty: {name}")
    if len(value) > max_len:
        fail(f"verify_math argument too long: {name}")
    return value


def tool_verify_math(root: Path, args: dict) -> str:
    mode = args.get("mode", "")
    allowed = {
        "rational",
        "gamma_quotient",
        "numeric_identity",
        "asymptotic_power",
        "perturbative_recurrence",
        "recurrence_transform",
    }
    if mode not in allowed:
        fail("unsupported verify_math mode")

    script = Path(__file__).resolve().with_name("verify_math.py")
    argv = [sys.executable, str(script)]

    if mode == "rational":
        argv += ["rational", "--expr", text_arg(args, "expr", max_len=1024)]
    elif mode == "gamma_quotient":
        argv += [
            "gamma-quotient",
            "--numerator",
            text_arg(args, "numerator", max_len=2048),
            "--denominator",
            text_arg(args, "denominator", max_len=2048),
        ]
    elif mode == "numeric_identity":
        argv += [
            "numeric-identity",
            "--lhs",
            text_arg(args, "lhs", max_len=2048),
            "--rhs",
            text_arg(args, "rhs", max_len=2048),
            "--tolerance",
            text_arg(args, "tolerance", required=False, default="1e-60", max_len=64),
        ]
    elif mode == "asymptotic_power":
        argv += [
            "asymptotic-power",
            "--n",
            text_arg(args, "n", max_len=8192),
            "--y",
            text_arg(args, "y", max_len=8192),
            "--powers",
            text_arg(args, "powers", required=False, default="1/2,1,3/2,2", max_len=1024),
        ]
    elif mode == "perturbative_recurrence":
        argv += [
            "perturbative-recurrence",
            f"--A={text_arg(args, 'A', max_len=2048)}",
            f"--B={text_arg(args, 'B', max_len=2048)}",
            f"--D={text_arg(args, 'D', max_len=2048)}",
            f"--u={text_arg(args, 'u', max_len=2048)}",
        ]
    elif mode == "recurrence_transform":
        argv += [
            "recurrence-transform",
            f"--A={text_arg(args, 'A', max_len=2048)}",
            f"--B={text_arg(args, 'B', max_len=2048)}",
            f"--D={text_arg(args, 'D', max_len=2048)}",
            f"--raw-var={text_arg(args, 'raw_var', max_len=64)}",
            f"--site-var={text_arg(args, 'site_var', max_len=64)}",
            f"--raw-in-site={text_arg(args, 'raw_in_site', max_len=2048)}",
            f"--source-orientation={text_arg(args, 'source_orientation', max_len=64)}",
        ]

    return run_process(argv, root, timeout=30)


TOOL_IMPL = {
    "read": tool_read,
    "search": tool_search,
    "status": tool_status,
    "test": tool_test,
    "verify_math": tool_verify_math,
}


VERIFY_PROPERTIES = {
    "mode": {
        "type": "string",
        "enum": [
            "rational",
            "gamma_quotient",
            "numeric_identity",
            "asymptotic_power",
            "perturbative_recurrence",
            "recurrence_transform",
        ],
    },
    "expr": {"type": "string"},
    "numerator": {"type": "string"},
    "denominator": {"type": "string"},
    "lhs": {"type": "string"},
    "rhs": {"type": "string"},
    "tolerance": {"type": "string"},
    "n": {"type": "string"},
    "y": {"type": "string"},
    "powers": {"type": "string"},
    "A": {"type": "string"},
    "B": {"type": "string"},
    "D": {"type": "string"},
    "u": {"type": "string"},
    "raw_var": {"type": "string"},
    "site_var": {"type": "string"},
    "raw_in_site": {"type": "string"},
    "source_orientation": {
        "type": "string",
        "enum": ["current_minus_neighbors_equals_mu", "verifier_canonical"],
    },
}


TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "read",
            "description": "Read at most 200 lines from a relative text file in the workspace.",
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {"type": "string"},
                    "start": {"type": "integer"},
                    "count": {"type": "integer"},
                },
                "required": ["path"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "search",
            "description": "Literal text search inside the research workspace. Paths must be relative.",
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {"type": "string"},
                    "path": {"type": "string"},
                    "max_results": {"type": "integer"},
                },
                "required": ["query"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "status",
            "description": "Inspect Git branch/status/diff statistics without modifying files.",
            "parameters": {"type": "object", "properties": {}},
        },
    },
    {
        "type": "function",
        "function": {
            "name": "test",
            "description": "Run a fixed whitelisted validation target.",
            "parameters": {
                "type": "object",
                "properties": {
                    "target": {
                        "type": "string",
                        "enum": ["cargo_check", "cargo_test", "cargo_fmt_check", "verifier_tests"],
                    }
                },
                "required": ["target"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "verify_math",
            "description": (
                "Deterministically audit rational/Gamma/numerical/asymptotic algebra, "
                "exact affine recurrence index transforms, and formal mu perturbation. "
                "For recurrence_transform prove the raw-to-site variable substitution "
                "before perturbative_recurrence. Numerical/asymptotic matches are not proof."
            ),
            "parameters": {
                "type": "object",
                "properties": VERIFY_PROPERTIES,
                "required": ["mode"],
            },
        },
    },
]


SYSTEM_PROMPT = r"""
You are an independent mathematical research auditor operating on a read-only
Rust numerical research workspace concerning the Riemann hypothesis.

Do not claim a proof of the Riemann hypothesis. Inspect evidence and classify
every important step as exact algebra, asymptotic reasoning, numerical evidence,
or unresolved.

Rules:
1. Never modify files, execute arbitrary shell commands, or use the Internet.
2. Paths supplied to tools must be relative to the workspace.
3. Keep read requests to at most 200 lines and use search to locate sections.
4. Use verify_math rational for fragile rational/index arithmetic.
5. Use gamma_quotient for Gamma recurrence claims. UNRESOLVED is not proof.
6. Use asymptotic_power before selecting a finite-size exponent from data.
7. Before perturbative_recurrence, use recurrence_transform when the source
   recurrence changes from a raw/Favard index to a shifted/staggered site index.
8. Normalize the source sign convention exactly before formal perturbation.
9. A CANDIDATE_U_FAILS verdict is a hard contradiction; revise the candidate.
10. Never upgrade numerical/asymptotic verifier output into symbolic proof.
11. If deterministic verification contradicts hand algebra, the hand algebra is
    refuted until the discrepancy is located and explained.
""".strip()


def public_assistant_message(message: dict) -> dict:
    result = {"role": "assistant", "content": message.get("content", "") or ""}
    calls = message.get("tool_calls", []) or []
    if calls:
        result["tool_calls"] = calls
    return result


def call_ollama(model: str, messages: list[dict], *, tools: bool = True) -> dict:
    payload = {
        "model": model,
        "stream": False,
        "messages": messages,
        "options": {"temperature": 0},
    }
    if tools:
        payload["tools"] = TOOLS

    last_error: Exception | None = None
    for attempt in range(OLLAMA_RETRIES + 1):
        request = urllib.request.Request(
            OLLAMA_URL,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
        )
        try:
            with urllib.request.urlopen(request, timeout=OLLAMA_TIMEOUT) as response:
                return json.load(response)
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")
            last_error = RuntimeError(
                f"Ollama HTTP {exc.code} {exc.reason}; url={OLLAMA_URL}; body={body[:4000]}"
            )
            # 4xx usually means a malformed request/model protocol and should not
            # be retried blindly, except 408/429 which can be transient.
            if exc.code not in {408, 429}:
                break
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = exc

        if attempt < OLLAMA_RETRIES:
            time.sleep(min(2**attempt, 8))

    raise RuntimeError(f"Ollama request failed after {attempt + 1} attempt(s): {last_error}")


def append_jsonl(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")


def load_resume_transcript(path: Path, *, expected_root: Path, expected_model: str):
    if not path.is_file():
        fail(f"resume transcript does not exist: {path}")
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            raw = raw.strip()
            if not raw:
                continue
            try:
                record = json.loads(raw)
            except json.JSONDecodeError as exc:
                fail(f"invalid JSONL at line {line_no}: {exc}")
            if isinstance(record, dict):
                records.append(record)

    starts = [record for record in records if record.get("event") == "start"]
    if len(starts) != 1:
        fail(f"resume transcript must contain exactly one start record, found {len(starts)}")
    start = starts[0]
    if Path(start.get("root", "")).resolve() != expected_root.resolve():
        fail("resume transcript workspace mismatch")
    if start.get("model") != expected_model:
        fail("resume transcript model mismatch")
    task = start.get("task", "")
    if not isinstance(task, str) or not task.strip():
        fail("resume transcript has no valid original task")

    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": task},
    ]
    ledger = ProofLedger()
    completed_turns = 0

    for record in records:
        event = record.get("event")
        if event == "assistant":
            message = {"role": "assistant", "content": record.get("content", "") or ""}
            calls = record.get("tool_calls", []) or []
            if calls:
                message["tool_calls"] = calls
            messages.append(message)
            turn = record.get("turn", 0)
            if isinstance(turn, int):
                completed_turns = max(completed_turns, turn)
        elif event == "tool":
            name = record.get("tool", "")
            output = record.get("output", "")
            messages.append({"role": "tool", "tool_name": name, "content": output})
            if name == "verify_math":
                mode = (record.get("arguments") or {}).get("mode", "")
                if mode:
                    ledger.add_verifier_output(mode, output)

    if completed_turns < 1:
        fail("resume transcript contains no completed assistant turn")
    return task, messages, completed_turns, ledger


def prepare_resume_transcript(source: Path, destination: Path) -> None:
    source = source.resolve()
    destination = destination.resolve()
    if source == destination:
        return
    if destination.exists():
        fail(f"refusing to overwrite existing continuation transcript: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")


def execute_tool(root: Path, call: dict) -> tuple[str, str]:
    function = call.get("function", {})
    name = function.get("name", "")
    args = function.get("arguments", {})
    if name not in TOOL_IMPL:
        return name, f"REFUSED UNKNOWN TOOL: {name}"
    if not isinstance(args, dict):
        return name, "REFUSED: tool arguments must be an object"
    try:
        output = TOOL_IMPL[name](root, args)
    except Exception as exc:  # Tool failures are evidence shown back to the model.
        output = f"TOOL ERROR: {exc}"
    return name, limited(output)


def run_agent(
    *,
    root: Path,
    model: str,
    task: str | None,
    max_turns: int,
    transcript: Path,
    resume_from: Path | None = None,
) -> int:
    if resume_from is not None:
        task, messages, completed_turns, ledger = load_resume_transcript(
            resume_from,
            expected_root=root,
            expected_model=model,
        )
        prepare_resume_transcript(resume_from, transcript)
        append_jsonl(
            transcript,
            {
                "event": "resume",
                "model": model,
                "root": str(root),
                "completed_turns": completed_turns,
                "resumed_from": str(resume_from.resolve()),
                "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
            },
        )
    else:
        completed_turns = 0
        ledger = ProofLedger()
        messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": task},
        ]
        append_jsonl(
            transcript,
            {
                "event": "start",
                "model": model,
                "root": str(root),
                "task": task,
                "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
            },
        )

    for additional_turn in range(1, max_turns + 1):
        turn = completed_turns + additional_turn
        print(f"\n=== AGENT TURN {turn} ===")
        try:
            result = call_ollama(model, messages)
        except Exception as exc:
            append_jsonl(
                transcript,
                {
                    "event": "ollama_error",
                    "turn": turn,
                    "model": model,
                    "error": str(exc),
                    "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
                },
            )
            print(f"OLLAMA_ERROR: {exc}")
            print(f"TRANSCRIPT={transcript}")
            return 5

        message = result.get("message", {})
        content = message.get("content", "") or ""
        calls = message.get("tool_calls", []) or []
        if content:
            print("--- assistant ---")
            print(content)
        messages.append(public_assistant_message(message))
        append_jsonl(
            transcript,
            {
                "event": "assistant",
                "turn": turn,
                "content": content,
                "tool_calls": calls,
                "prompt_eval_count": result.get("prompt_eval_count"),
                "eval_count": result.get("eval_count"),
            },
        )

        if not calls:
            print("\n=== DETERMINISTIC LEDGER ===")
            print(ledger.public_summary())
            print("\n=== AGENT FINAL ===")
            print(content)
            print(f"TRANSCRIPT={transcript}")
            return 0

        for call in calls:
            name, output = execute_tool(root, call)
            print(f"\n--- tool: {name} ---")
            print(output)
            messages.append({"role": "tool", "tool_name": name, "content": output})
            arguments = call.get("function", {}).get("arguments", {}) or {}
            append_jsonl(
                transcript,
                {
                    "event": "tool",
                    "turn": turn,
                    "tool": name,
                    "arguments": arguments,
                    "output": output,
                },
            )
            if name == "verify_math" and isinstance(arguments, dict):
                mode = arguments.get("mode", "")
                if mode:
                    ledger.add_verifier_output(mode, output)

    print("\nMAX_TURNS_REACHED_WITHOUT_FINAL_ANSWER")
    print(ledger.public_summary())
    print(f"TRANSCRIPT={transcript}")
    return 4


def main() -> None:
    parser = argparse.ArgumentParser(description="Read-only multi-turn Riemann research agent")
    parser.add_argument("--root", default=".")
    parser.add_argument("--model", default="qwen3.8:latest")
    parser.add_argument("--task")
    parser.add_argument("--max-turns", type=int, default=12)
    parser.add_argument("--transcript")
    parser.add_argument("--resume-from")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    if not root.is_dir():
        fail(f"workspace root does not exist: {root}")
    if not 1 <= args.max_turns <= 30:
        fail("--max-turns must be between 1 and 30")
    if args.resume_from and args.task:
        fail("--task cannot be combined with --resume-from")
    if not args.resume_from and not args.task:
        fail("--task is required unless --resume-from is used")

    resume_from = Path(args.resume_from).resolve() if args.resume_from else None
    if resume_from is not None and not resume_from.is_file():
        fail(f"resume transcript does not exist: {resume_from}")

    if args.transcript:
        transcript = Path(args.transcript).resolve()
    else:
        stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
        transcript = root / "agent_runs" / f"research-{stamp}.jsonl"

    rc = run_agent(
        root=root,
        model=args.model,
        task=args.task,
        max_turns=args.max_turns,
        transcript=transcript,
        resume_from=resume_from,
    )
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
