# Verified research-agent harness

This repository contains a read-only local research harness for adversarial
mathematical experiments. It is designed to make fragile algebra auditable and
to keep numerical evidence separate from proof-level claims.

The harness does **not** claim a proof of the Riemann hypothesis.

## Components

- `tools/verify_math.py` — deterministic verifier with restricted expression grammars.
- `tools/proof_ledger.py` — machine-readable evidence classification and final gates.
- `tools/riemann_research_agent.py` — resilient single-agent Ollama runner.
- `tools/riemann_dual_research_agent.py` — researcher/critic harness with a fail-closed final gate.

## Deterministic verifier modes

### `rational`

Exact integer/rational arithmetic using an AST whitelist. Model-supplied Python
is never evaluated.

### `gamma-quotient`

Reduces positive rational Gamma arguments by recurrence. Quarter/half bases are
also simplified with the exact identities

```text
Gamma(1/2) = sqrt(pi)
Gamma(1/4) Gamma(3/4) = pi sqrt(2)
```

The output explicitly distinguishes an exact reduction from unresolved Gamma
bases.

### `numeric-identity`

Independent 80-digit `mpmath` evaluation. A numerical match is labelled as
numerical evidence and is never promoted to symbolic proof.

### `asymptotic-power`

Fits candidate models `y_n = C + a n^{-p}` and reports the best candidate power
by SSE, together with successive-difference ratios. The result is explicitly
labelled numerical/asymptotic evidence.

### `recurrence-transform`

Verifies an affine change from a raw recurrence index to one common site index
and normalizes the recurrence sign convention. This prevents accidental mixing
of coefficients expressed in different variables.

Example:

```bash
python tools/verify_math.py recurrence-transform \
  --A='(d+1/2)*(d+3/2)' \
  --B='d*(d-1)' \
  --D='4*d+1' \
  --raw-var=d \
  --site-var=j \
  --raw-in-site='2*j-1/4' \
  --source-orientation=current_minus_neighbors_equals_mu
```

### `perturbative-recurrence`

Formally expands

```text
A(j)*(g[n+1]-g[n])
-B(j)*(g[n]-g[n-1])
+mu*D(j)*g[n] = 0
```

with

```text
g[n] = 1 + mu*u[n] + mu^2*v[n] + O(mu^3).
```

The candidate `u[n]` must make the first-order residual exactly zero. Otherwise
the verifier returns `CANDIDATE_U_FAILS`.

## Proof ledger

Every `verify_math` result can be classified as one of:

- `proved_exact`
- `asymptotic_evidence`
- `numerical_evidence`
- `unresolved`
- `refuted`
- `unknown`

The dual runner uses this ledger mechanically. For the blind mu^2 workflow the
final synthesis is withheld unless the final phase has executed:

- `recurrence_transform`
- `perturbative_recurrence`
- `asymptotic_power`
- `gamma_quotient`

and the recurrence transform plus first-order perturbative candidate have exact
success verdicts.

## Read-only tool boundary

Agents may only:

- read bounded line ranges;
- run literal searches inside the workspace;
- inspect Git status;
- run fixed validation targets;
- call deterministic `verify_math` modes.

They cannot modify files, run arbitrary shell commands, escape the workspace,
or access the Internet through the harness.

## Ollama resilience

The single-agent runner defaults to a 1800-second request timeout and records
HTTP response bodies for protocol errors. Timeout/retry behavior can be changed
with:

```text
RIEMANN_AGENT_TIMEOUT
RIEMANN_AGENT_RETRIES
OLLAMA_URL
```

A request failure is appended to the JSONL transcript before the runner exits,
so a long experiment is recoverable rather than silently lost.

## Transcript resume

A previous single-agent JSONL run can be replayed without re-executing its tools:

```bash
python tools/riemann_research_agent.py \
  --root /path/to/workspace \
  --model qwen3.8:latest \
  --resume-from agent_runs/previous.jsonl \
  --transcript agent_runs/continued.jsonl \
  --max-turns 8
```

Only public assistant messages, structured tool calls and recorded tool outputs
are restored. Hidden model reasoning is neither stored nor reconstructed.

## Dual-agent run

Given a challenge file in the workspace:

```bash
python tools/riemann_dual_research_agent.py \
  --root /path/to/workspace \
  --researcher qwen3.8:latest \
  --critic nemotron-3.5-lightning:30b \
  --rounds 4 \
  --max-tool-turns 10 \
  --challenge BLIND_MU2_CHALLENGE.md
```

The final synthesis fails closed if the deterministic evidence requirements are
not met.

## CI

GitHub Actions installs the pinned research dependency, compiles all harness
modules and runs verifier/ledger regression tests in addition to the existing
Rust checks.
