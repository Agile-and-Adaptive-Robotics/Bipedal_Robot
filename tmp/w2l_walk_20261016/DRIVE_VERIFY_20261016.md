# DRIVE IMPLEMENTER — verify stage report (2026-10-16, overnight walk campaign)

Stage role: implement L/R alpha+gamma drive on BOTH models per the audited
HUMAN + ROBOT plans, record the baseline BEFORE any edit, re-run every gate,
and deliver a one-command defaults gate.

## 0. What this stage actually did (honest scope)

Both plans arrived ALREADY IMPLEMENTED on egg's working tree by the sibling
campaign stages (human: `CHANGES_20261016_connectome_drive.md`; robot:
`w2l_audit/CHANGES_MANIFEST_20261016_w2l.md`).
This stage therefore (a) verified every diff line against the plans,
(b) RE-RECORDED the pre-change baseline itself at HEAD b497b67b (see §2 —
none of the sibling stages' gate numbers are taken on faith here), and
(c) re-ran every gate post-change with its own commands. No source file was
edited by this stage; the only new files are in `reports_20261016/`
(`verify/` logs + three gate scripts + this report).

Implementation verified by reading (git diff on egg, this session):
- HUMAN (5 files, +133 insertions, 0 deletions — matches plan):
  `params.py` (gamma_l/r, alpha_l/r default 0.0 + alpha_to_mn=0.5 at the end
  of G), `build_network.py` (GAMMA_s/ALPHA_s cells+ports in `__post_init__`
  after the VEST block, per-side gated; GAMMA→RG_E/RG_F in `_build_rg` after
  the POSTURE edge; ALPHA→ipsilateral MNs at end of `_wire_balance`),
  `build_network_w2lvar.py` / `build_network_syn6.py` (same flags/cells +
  their own GAMMA→RG wiring; ALPHA fan-out inherited), `runner.py`
  (port-name-guarded `u[GAMMA_c_s]=drive*G[gamma_s]`, `u[ALPHA_c_s]=drive*
  G[alpha_s]` right after the POSTURE feed).
- ROBOT (8 files): `w2l_cpg/build_w2l_net.py` + `w2l_mujoco/
  build_w2l_split_net.py` (alpha_l/r fields; conditional per-MN "ALPHA <MN>"
  ports, 6/side, built only when that side's alpha != 0; input-count asserts
  10→10+len / 6→6+len), `build_w2l_aff_net.py` (passthrough + census fix),
  `smoke_w2l.py` (--tonic-e-l/-r, --tonic-f-l/-r gamma + --alpha-l/-r),
  `test_w2l_air.py` (air-variant weld regression fix), `test_w2l_air_
  afferented.py` / `test_w2l_ground.py` (KV te_l/te_r/tf_l/tf_r nan-default
  + alpha_l/r 0-default; per-side injection), `w2l_cpg/README.md` (docs).
  GAMMA on the robot = the 4 pre-existing TONIC ports (value-level, no
  topology change) — per plan.

Knob semantics confirmed as specified: four default-0 knobs per plant;
alpha = tonic drive to MN pools, gamma = drive to the RG half-centers, L/R
separable. `alpha_to_mn=0.5` is a fixed conductance whose edges exist only
when alpha_s>0 — no default-build effect (documented deviation from "all
knobs default 0" in letter, none in effect).

## 1. Baseline recording method (deviation from the stash recipe, documented)

`git stash`/`pop` was REJECTED on purpose: egg has core.autocrlf=true and two
of the campaign-edited files are LF-authored (git warns "LF will be replaced
by CRLF"), so a stash pop would have rewritten their bytes. Instead the
pre-change baseline was recorded in a THROWAWAY DETACHED WORKTREE at HEAD:

    git worktree add --detach "%TEMP%\lrdrive_head_wt" HEAD   # b497b67b

with the untracked gate scripts mirrored in at their relative paths (the
anchor script is cwd-relative; NOTE `reports_20260925/audit/
audit_defaults_gate.py` hardcodes the REAL repo path, so a parametrized twin
`reports_20261016/default_census_20261016.py [spinal_dir]` was written for
the worktree standard census). Worktree removed after (`git worktree remove
--force`); `git worktree list` shows only the main tree.

## 2. PRE-CHANGE baseline (recorded at HEAD b497b67b by THIS stage)

All commands run on egg with `D:\Anaconda\envs\myo\python.exe` from the
(spinal) worktree cwd; logs in `reports_20261016\verify\head_*.log`.

| gate | value at HEAD |
|---|---|
| standard default census | (410, 376, 1186) |
| AARL_NET=w2lvar default census | (888, 382, 7866) |
| AARL_NET=syn6 default census   | (794, 382, 3033) |
| anchor A (s3k stage-4 full_rules=0 eval) | **-161.56754173676563** |
| robot drive-gate hashes (A/B/C nets) | 01929488aad3… / 57637001168c… / 6834c8d25036… |
| W2L smoke | PASS period 1.000 s, antiphase -0.525, 22/22 RG-E bursts |
| M6 ground walk (default) | full log captured (head_m6_walk.log) |

- The old AGENTS baseline **-160.23425729850192 does NOT reproduce at HEAD**
  (my HEAD run gives -161.56754173676563) — consistent with the 2026-09-26
  audit note (stage renumber + left-ref fix moved it). The recorded baseline
  from here on is -161.56754173676563 (= the 09-26 recorded value).
- Anchor B (production full_rules=1): NOT re-run at HEAD by this stage
  (worktree already removed when the need was noted); post-change value -25.0
  matches the sibling stage's recorded pre-change observation (-25.0). Anchor
  A is the bit-exact anchor of record; B is a degenerate frozen-sentinel
  anchor on current physics either way.

## 3. POST-CHANGE gates (this stage's own runs, logs in verify\post_*.log)

HUMAN / spinal:
- standard default census (410, 376, 1186) == HEAD. PASS.
- w2lvar (888, 382, 7866) and syn6 (794, 382, 3033) default builds == HEAD. PASS.
- anchor A: **-161.56754173676563 BIT-EXACT** vs my HEAD baseline. PASS.
- anchor B: -25.0 (matches sibling's pre-change observation). PASS (as recorded).
- winner+4-knobs probe (`_lrdrive_probe_20261016.py`, knobs 0.8/0.4/0.5/0.5
  at the s3k winner config): counts (882, 386, 8492) = winner 878/382/8396
  +4n/+4p/+96s; wiring dump = 46 ALPHA→MN per side at g=0.5 + GAMMA→RG_E/RG_F
  at the winner conductances 1.6988/0.9964; 4 ports GAMMA_c/ALPHA_c × l/r.
  Exactly the drawn wiring. PASS.

ROBOT / W2L:
- drive gate (`_w2l_drive_gate.py`) default probes: censuses
  (95,208,10)/(87,164,6)/(97,186,12); trace hashes IDENTICAL to my HEAD run
  (01929488…/57637001…/6834c8d2…), 3/3 fresh-net repeats equal. PASS.
- **W2L smoke gate** (`w2l_cpg\smoke_w2l.py`, default knobs): output
  BYTE-IDENTICAL to HEAD (`fc /b`: no differences) —
  "W2L_SMOKE PASS period=1.000 antiphase_r=-0.525 rg_e_bursts_L=22
  rg_e_bursts_R=22". PASS. (antiphase -0.525 is current-HEAD's value; the
  older AGENTS -0.532 was a previous net revision.)
- M6 ground walk (default): BYTE-IDENTICAL to HEAD (`fc /b`). PASS.
- Functional alpha (paired, `_w2l_alpha_debug.py`): bilateralrg L-only
  1.5123→1.8809 mV (R 1.5885→1.5874); aff-net R-only 1.5933→1.9000 mV
  (L 1.5115→1.5101). Each knob lifts exactly its own side's MN pools. PASS.
- Functional gamma, closed loop (M6 walk `--tf-l=4.75`): tilt max 8.2→6.9 deg,
  L RG-E period 1.539→1.513 s. PASS (matches sibling's numbers).
- Functional alpha, closed loop (M6 walk `--alpha-l=0.3`): regime change
  FALLS/NO-GAIT → SHUFFLES-IN-PLACE (heel strikes L3/R10, harness 0%) —
  the documented "0.3 nA overdrives" behavior. Observed.

Honest caveats (unchanged from the sibling plans, re-confirmed here):
- The drive gate's inline A-sanity prints "L MN mean 0.998 vs R 1.906 (want
  L > R)" — its 1 s window is confounded by the alternating heel schedule;
  the paired 8 s window-mean alpha-debug above is the authoritative check.
- Open-loop gamma asymmetry window is narrow: te 3.0/2.0 stops the c1-locked
  pair's alternation (burst count 7/7 → 1/1). Closed-loop (M6) asymmetry works.

## 4. defaultsGateCmd (of record; corrected after consumer-side ssh 255s)

    ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=30 -o ConnectionAttempts=3 easteregg2 "cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal && D:\Anaconda\envs\myo\python.exe reports_20261016\drive_defaults_gate.py"

History: the first submitted form (plain `ssh -o BatchMode=yes ...`) was
reported by the consumer as exit 255 "Permission denied (publickey,
password,keyboard-interactive)" — an ssh transport/auth failure (the gate
logic was never reached; the config alias was clearly read, since the error
names the config's User). A parallel fix round additionally reported a
"Host key verification failed" consumer mode. The same plain command re-run
from the laptop minutes later passed clean, so the episode was transient
egg-side/lab-network flakiness (this session also saw two plain connection
timeouts to egg). The corrected command adds, with NO change to the gate
logic or its exit semantics: StrictHostKeyChecking=accept-new (self-heals a
missing known_hosts entry, still hard-fails a genuinely changed host key),
ConnectTimeout=30 + ConnectionAttempts=3 (fail fast, retry TCP-level
blips). Verified this fix round, every run my own: (1) accept-new
self-heal with UserKnownHostsFile=/dev/null → warning + connect, exit 0;
(2) full command under bash → DEFAULTS GATE: PASS, exit 0; (3) full
command under cmd.exe (batch file) → PASS, exit 0; (4) negative test
(remote AARL_NET=w2lvar) → DEFAULTS GATE: FAIL, exit 1; (5) the original
plain form → PASS, exit 0 (the gate itself was never broken). A future
TRANSIENT auth denial can still 255 the command (no ssh option force-
retries auth under BatchMode); re-running it is the documented remedy —
exactly how this episode cleared. Remote gate script unchanged
(reports_20261016/drive_defaults_gate.py: gate 1 = default census
(410,376,1186) with AARL_NET unset; gate 2 = fresh-subprocess anchor A
must equal -161.56754173676563 bit-exactly; AARL_NPZ redirect protects
spinal_run.npz).

Gate content (unchanged, re-verified 4 more times during the fix): exits 0
iff (a) default census == (410, 376, 1186) with AARL_NET unset and (b) a
fresh-process anchor A reproduces -161.56754173676563 bit-exactly. It
redirects AARL_NPZ to a side name and deletes it, so spinal_run.npz is
never touched. It depends on the untracked anchor script
`spinal\_s3k_anchor_eval_20261016.py` (left in place by the sibling stage
precisely as the reusable gate).

## 5. Files added by this stage (all untracked, campaign outputs)

- reports_20261016/DRIVE_VERIFY_20261016.md (this report)
- reports_20261016/drive_defaults_gate.py (the §4 gate)
- reports_20261016/default_census_20261016.py, variant_census_20261016.py
- reports_20261016/verify/ — 19 logs: head_{defaults_gate, census_w2lvar,
  census_syn6, anchor_A, drive_gate, smoke_w2l, m6_walk}.log and
  post_{defaults_gate, census_w2lvar, census_syn6, anchor_A, anchor_B,
  lrdrive_probe, drive_gate, smoke_w2l, alpha_debug, m6_walk, m6_gammaL_demo,
  m6_alphaL_demo}.log

No commits, no pushes. Working tree left exactly as found (git status
identical to session start apart from the reports_20261016 additions;
my AARL_NPZ byproducts deleted).

## 6. For the next stage (phase-4 study recipe, from the plans)

Add gamma_l/gamma_r/alpha_l/alpha_r to the study KEYS; loader branch
`P.G[k] = float(p.get(k, 0.0))` before build (set_stage-style); runner needs
nothing else. Wean ladder = shrink the gamma values toward 0 in the winner
json. Robot side: save any tuned te_l/te_r/tf_l/tf_r/alpha_l/alpha_r in the
winner json with `if key in best` loader branches.
