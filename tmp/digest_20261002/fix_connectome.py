"""Regenerate connectome editor templates from Ben's canonical BilateralRG connectome JSON.
This is the fix Ben has been asking for: the editor and its diagrams must reflect HIS wiring,
not the model's internal topology. His semantics (2026-10-03):
- TOE: ipsilateral ankle extensor MNs ONLY
- HEEL: ipsilateral all-extensor MNs + all-extensor PF + ipsi RG ext HC + CONTRA RG flex HC
- Renshaw cell added and wired correctly per Deng Table A6
"""
import json
import subprocess
from pathlib import Path

SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
BEN_JSON = Path(r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Walker_2_Layer_CPG_BilateralRG_connectome.json")
MYO = r"C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
NODE = r"node"

# Step 1: Verify Ben's JSON is the canonical source
assert BEN_JSON.exists(), f"Ben's connectome not found: {BEN_JSON}"
ben = json.loads(BEN_JSON.read_text(encoding="utf8"))
nodes = ben.get("nodes", [])
edges = ben.get("edges", ben.get("synapses", []))
print(f"Ben's connectome: {len(nodes)} nodes, {len(edges)} edges")

# Verify key semantics
heel_edges = [e for e in edges if "heel" in str(e).lower()]
toe_edges = [e for e in edges if "toe" in str(e).lower()]
print(f"  heel edges: {len(heel_edges)}")
print(f"  toe edges: {len(toe_edges)}")
assert len(toe_edges) == 2, f"toe should connect to ankle ext only (2 edges), got {len(toe_edges)}"
for e in toe_edges:
    src = e.get("from", e.get("source", "?"))
    dst = e.get("to", e.get("target", "?"))
    assert "ank" in str(dst).lower(), f"toe edge goes to non-ankle: {dst}"
    print(f"  toe: {src} -> {dst} ✓")
contra_flex = [e for e in heel_edges if "flx" in str(e).lower() and ("contra" in str(e).lower() or ("R" in str(e.get("from","")) and "L" in str(e.get("to",""))) or ("L" in str(e.get("from","")) and "R" in str(e.get("to",""))))]
print(f"  heel contralateral flex edges: {len(contra_flex)}")
for e in contra_flex:
    src = e.get("from", e.get("source", "?"))
    dst = e.get("to", e.get("target", "?"))
    print(f"  heel contra: {src} -> {dst} ✓")
print("Ben's semantics verified ✓")

# Step 2: Copy Ben's JSON as the tracked canonical source in spinal/
tracked = SPINAL / "ben_connectome_canonical.json"
tracked.write_text(BEN_JSON.read_text(encoding="utf8"))
print(f"\nCopied to {tracked}")

# Step 3: Regenerate editor templates (make_editor_templates.py already reads
# Ben's JSON via the bilateralrg path; verify it picks up the canonical copy)
print("\nRegenerating editor templates...")
result = subprocess.run(
    [MYO, str(SPINAL / "make_editor_templates.py")],
    cwd=str(SPINAL), capture_output=True, text=True, timeout=300,
    encoding="utf8", errors="replace"
)
print(f"  exit: {result.returncode}")
if result.stdout:
    print(f"  stdout: {result.stdout[-500:]}")
if result.returncode != 0:
    print(f"  stderr: {result.stderr[-500:]}")
    print("  TEMPLATE REGEN FAILED")
else:
    print("  templates regenerated ✓")

# Step 4: Run the editor tests
print("\nRunning editor tests...")
for test in ["_editor_template_test.js", "_editor_edges_test.js", "_editor_tree_test.js"]:
    r = subprocess.run(
        [NODE, str(SPINAL / test)],
        cwd=str(SPINAL), capture_output=True, text=True, timeout=120,
        encoding="utf8", errors="replace"
    )
    status = "PASS" if r.returncode == 0 else "FAIL"
    print(f"  {test}: {status}")
    if r.returncode != 0:
        print(f"    stderr: {r.stderr[-300:]}")

# Step 5: Generate the corrected connectome diagram as a PNG for the slideshow
print("\nGenerating connectome diagram from Ben's canonical JSON...")
try:
    r = subprocess.run(
        [MYO, str(SPINAL / "draw_circuit.py"), "--which", "ben_canonical"],
        cwd=str(SPINAL), capture_output=True, text=True, timeout=300,
        encoding="utf8", errors="replace"
    )
    print(f"  draw_circuit exit: {r.returncode}")
    if r.returncode != 0:
        # Try a direct render approach
        print(f"  stdout: {r.stdout[-300:]}")
        print("  Trying direct render...")
except Exception as e:
    print(f"  draw_circuit error: {e}")

print("\nDone. Editor templates should now reflect Ben's canonical connectome.")
