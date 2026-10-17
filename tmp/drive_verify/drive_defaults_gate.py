r"""DEFAULTS GATE for the 2026-10-16 L/R alpha+gamma drive implementation
(overnight walk campaign, drive-implementer verify stage).

Exits 0 iff BOTH hold on the current working tree, with all four drive
knobs at their params.G defaults (gamma_l=gamma_r=alpha_l=alpha_r=0):
  1. Default build census == (410, 376, 1186) neurons/inputs/synapses
     (standard net, AARL_NET unset).
  2. The recorded baseline anchor eval reproduces BIT-EXACTLY:
     ANCHOR[A] kine_score == -161.56754173676563
     (recipe: _s3k_anchor_eval_20261016.py A = AGENTS 2026-09-26 anchor —
     set_stage(4, {v10 multipliers + curriculum_stage3.json params,
     full_rules:0}) + runner --eval at the winner drive. Baseline
     RE-RECORDED AT HEAD b497b67b by this verify stage in a detached
     git worktree (reports_20261016/verify/head_anchor_A.log), matching
     the 09-26 recorded value; the older -160.23425729850192 is stale on
     current HEAD per the 09-26 audit note.)

The anchor runs in a FRESH subprocess (params.G is module-global; one
flavor per process keeps the merge unambiguous) with AARL_NPZ redirected
so spinal_run.npz is never touched. The redirect npz is deleted on exit.

Run (from anywhere on easteregg2):
  D:\Anaconda\envs\myo\python.exe D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\drive_defaults_gate.py
"""
import os
import subprocess
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SPINAL = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, SPINAL)
os.chdir(SPINAL)

EXPECTED_COUNTS = (410, 376, 1186)
EXPECTED_ANCHOR = "-161.56754173676563"
NPZ_SIDENAME = "_defaults_gate_run.npz"

fail = []

# --- gate 1: default-build census ----------------------------------------
aarl = os.environ.get("AARL_NET", "")
if aarl != "":
    fail.append(f"AARL_NET must be unset for the defaults gate (got {aarl!r})")
else:
    import muscle_map as MM
    acts = []
    for base in MM._GROUPS_BY_NAME:
        acts.append(base + "_r")
        acts.append(base + "_l")
    import build_network as BN
    net = BN.build(acts)
    n = net.net
    counts = (n.get_num_neurons(), n.get_num_inputs_actual(),
              n.get_num_connections())
    print(f"gate1 default counts neurons/inputs/synapses = {counts}")
    if counts != EXPECTED_COUNTS:
        fail.append(f"counts {counts} != {EXPECTED_COUNTS}")
    else:
        print("gate1 COUNTS PASS")

# --- gate 2: baseline anchor eval, bit-exact, fresh subprocess ------------
env = dict(os.environ)
env["PYTHONIOENCODING"] = "utf-8"
env["AARL_NPZ"] = NPZ_SIDENAME
env.pop("AARL_NET", None)
env.pop("AARL_KY", None)
env.pop("AARL_PELVIS_TY", None)
proc = subprocess.run(
    [sys.executable, os.path.join(SPINAL, "_s3k_anchor_eval_20261016.py"), "A"],
    cwd=SPINAL, env=env, capture_output=True, text=True, encoding="utf-8",
    errors="replace")
out = (proc.stdout or "") + (proc.stderr or "")
anchor = None
for line in out.splitlines():
    if line.startswith("ANCHOR[A] kine_score:"):
        anchor = line.split(":", 1)[1].strip()
print(f"gate2 ANCHOR[A] kine_score = {anchor!r} (expected {EXPECTED_ANCHOR!r})")
if proc.returncode != 0:
    fail.append(f"anchor subprocess exit {proc.returncode}")
if anchor is None:
    fail.append("anchor eval produced no ANCHOR[A] kine_score line")
elif anchor != EXPECTED_ANCHOR:
    fail.append(f"anchor {anchor!r} != {EXPECTED_ANCHOR!r} (NOT bit-exact)")
else:
    print("gate2 ANCHOR BIT-EXACT PASS")

# --- cleanup the redirect npz (best effort) -------------------------------
try:
    os.remove(os.path.join(SPINAL, NPZ_SIDENAME))
except OSError:
    pass

if fail:
    print("DEFAULTS GATE: FAIL")
    for f in fail:
        print("  - " + f)
    sys.exit(1)
print("DEFAULTS GATE: PASS (counts (410, 376, 1186) + anchor "
      "-161.56754173676563 bit-exact with knobs at defaults)")
sys.exit(0)
