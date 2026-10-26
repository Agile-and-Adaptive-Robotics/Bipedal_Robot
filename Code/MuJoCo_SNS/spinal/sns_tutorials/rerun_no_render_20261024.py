# Re-run sns-toolbox tutorial runners 1/2/4/8 headless WITHOUT the graphviz
# diagram render (dot.exe not on PATH on EB475WS4). Simulation outputs are
# the point for the response-matching task; rendering is cosmetic.
import runpy
import sys

import sns_toolbox.renderer as _r

def _render_skip(*args, **kwargs):
    print("[render skipped: graphviz dot not on PATH (EB475WS4 2026-10-24)]")

_r.render = _render_skip

HERE = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\sns_tutorials"
for name in ("tutorial1_run.py", "tutorial2_run.py",
             "tutorial4_run.py", "tutorial8_run.py"):
    print("=" * 20, "RUNNING", name, "=" * 20, flush=True)
    try:
        runpy.run_path(HERE + "\\" + name, run_name="__main__")
        print("=" * 20, name, "DONE", "=" * 20, flush=True)
    except Exception as e:
        print("=" * 20, name, "FAILED:", e, "=" * 20, flush=True)
print("ALL_RERUNS_ATTEMPTED")
