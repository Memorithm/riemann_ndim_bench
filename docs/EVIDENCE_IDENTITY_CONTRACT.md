# Evidence identity contract

`ProofLedger` records are admissible to deterministic gates only when one exact
verifier invocation is bound to all of the following identities:

- the canonical proposition text and its SHA-256 digest;
- the canonical JSON verifier arguments and their SHA-256 digest;
- the verifier filename and SHA-256 digest;
- the checked-out Git source commit;
- a unique execution identifier;
- an integer process exit code equal to zero; and
- the required output fields for the selected verifier mode, including an exact
  mode match and no duplicate keys.

Missing, malformed, non-zero, unbound, or mode-mismatched output is classified
`unknown` with authority `invalid_or_unbound` and cannot satisfy a gate. Gates
evaluate only a proposition/input pair explicitly submitted with
`gate_target=true`. The first such identity for a mode is locked. Later output
for an unrelated proposition cannot replace it, even if that output succeeds.

## Evidence authority

An exact success emitted by `verify_math.py` or `symbolic_mu2.py` is labelled
`script_verified_exact`. This means the pinned script verified the pinned input
at the pinned source revision. It is not a theorem review and is not a formal
ProofLab acceptance.

The ledger therefore exposes ProofLab acceptance separately. Local verifier
records have `prooflab_receipt_id = None`, and the public summary renders
`prooflab=not_submitted`. No local output string can manufacture a ProofLab
receipt. A future ProofLab integration must validate a signed receipt against
the same proposition, input, verifier, source and execution identities before
setting formal acceptance.

## Caller requirement

Every `verify_math` tool call must include a self-contained `proposition`
string. The runner hashes that statement and the exact arguments independently;
it also computes verifier and source digests rather than trusting stdout. New
transcripts persist the resulting execution and verifier identities so resumed
runs reconstruct the same fail-closed ledger. A fresh run receives a random
nonce that participates in every execution ID, and the runner refuses to append
a fresh run to an existing transcript. Legacy transcripts without their own
source SHA are rejected rather than rebound to the current checkout.
