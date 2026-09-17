"""Final LAMINATED curriculum run (2026-09-16): apply the stage-3 winner
EXACTLY as the curriculum objective evaluated it - via _curriculum.set_stage
(v10 multipliers base + stage-3 searched keys) - then the full 22 s ground
walk. Replaces _final_run.py for the laminated architecture (the old script
skipped several v10 multiplier overrides, so its config drifted from what
the study actually scored)."""
import io
import json
import shutil
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import _curriculum as C
import runner as R

st3 = json.loads(open("curriculum_stage3.json", encoding="utf-8").read())
v10 = json.loads(open("best_walk_params_v10.json", encoding="utf-8").read())
p = {**v10["multipliers"], **st3["params"]}
C.set_stage(3, p)
drive = float(st3["params"]["drive"])
print(f"LAMINATED curriculum final: study {st3.get('study')} "
      f"score {st3.get('score'):.3f} (trial {st3.get('trial')}), "
      f"drive {drive:.4f}, heel {p['heel_rge']:.2f}, toe {p['toe_rge']:.2f}, "
      f"ib {p['ib_rge']:.2f}, ia_in {p['ia_in']:.2f}, "
      f"pre {p['phase_reset_e']:.2f}/{p['phase_reset_f']:.2f}, "
      f"trim {p['ankle_post_walk_trim']:.2f}, renshaw 0.5, "
      f"kneext-inh {p.get('f1_kneext_inh', 0.0):.2f}", flush=True)

R.main(["--drive", repr(drive), "--time", "22"])
shutil.copy("spinal_run.npz", "curriculum_final.npz")
shutil.copy("spinal_run.npz", "ground_curriculum_final.npz")
print("saved curriculum_final.npz (+ ground_curriculum_final.npz for the "
      "hindlimb window auto-detect)")
