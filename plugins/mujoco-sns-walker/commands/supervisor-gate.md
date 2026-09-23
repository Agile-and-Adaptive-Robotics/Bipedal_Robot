---
description: Spawn the SUPERVISOR_AGENT.md gatekeeper at milestone M1/M2/M3 - mandatory before declaring results, editing figures/tex/circuit, or claiming literature fidelity
---

# /supervisor-gate

The repo's AGENTS.md mandates: walker sessions MUST spawn a general-purpose agent
gated by `Code\MuJoCo_SNS\spinal\SUPERVISOR_AGENT.md` at three milestones, and act on
BLOCK verdicts before proceeding. This command makes that one step.

## Milestones

- **M1** - before declaring any tuning/result "done" or "improved"
- **M2** - before editing any figure, the dissertation .tex, or the circuit
- **M3** - before calling any wiring/connectome change "literature-faithful"

## Procedure

1. Read `<repo>\Code\MuJoCo_SNS\spinal\SUPERVISOR_AGENT.md` FRESH from the repo and use
   its FULL TEXT verbatim as the spawn prompt (never an embedded copy - it must not
   drift). Prepend one line naming the milestone and the session context:
   `Milestone: M1. Context: <one paragraph - what was done, what claim is being gated,
   file paths involved>.`
2. Spawn it as a **general-purpose** agent (it must read repo files itself; it is
   instructed NOT to trust the session's claims).
3. Wait for the verdict. `VERDICT: SHIP` -> proceed. `VERDICT: BLOCK` -> fix the
   listed items first, then re-gate.

## Notes

- The supervisor checks figure standards, dissertation-tex protection, connectome
  ownership (every wiring change traces to Ben's spec or a cited literature rule),
  literature fidelity markings, Ben's standing instructions (ankle posture, trunk
  upright, per-leg metrics, ISB conventions), and honest reporting (score
  decomposition, reproducibility gates).
- Skipping the gate to save time is exactly the failure mode the role was created
  for (2026-09-20/21: days lost to a mis-wired circuit + unauthorized edits).
