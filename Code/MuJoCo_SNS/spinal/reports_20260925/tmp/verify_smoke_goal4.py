"""goal4 wiring smoke verification (2026-09-25).

Checks, with evidence printed:
1. the npz-column-verification line appears in each smoke log;
2. each fresh db holds its fresh study with 2 COMPLETE trials, values
   printed (sentinel inversion check: seed trial must be a real air
   score > the -10 rhythm-gate floor);
3. the stage jsons exist and parse;
4. the chain-private npz files exist and their neuro_names / key_joints
   differ per variant (proves AARL_NET routed the builder AND that the
   chains did not share an npz);
5. the protected artifacts are untouched (spinal_run.npz mtime/size,
   optuna_walk.db mtime, 20 existing studies).
"""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import numpy as np
import optuna

os.chdir(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")

for v in ("w2lvar", "syn6"):
    print(f"===== {v} =====")
    log = open(rf"reports_20260925\logs\smoke_{v}_s1.log",
               encoding="utf-8", errors="replace").read()
    hit = [ln for ln in log.splitlines() if "npz columns verified" in ln]
    print("log column-verification line:", hit[0] if hit else "MISSING")

    study = optuna.load_study(
        study_name=f"curr_{v}_s1_air_deaff",
        storage=f"sqlite:///optuna_{v}.db")
    trials = study.trials
    print(f"study curr_{v}_s1_air_deaff: {len(trials)} trials, "
          f"states={[str(t.state).split('.')[-1] for t in trials]}")
    for t in trials:
        print(f"  trial {t.number}: value={t.value!r} params={t.params}")

    js = json.loads(open(f"curriculum_{v}_stage1.json",
                         encoding="utf-8").read())
    print("stage json:", {k: js[k] for k in
                          ("stage", "score", "trial", "study", "db",
                           "npz", "aarl_net")})

    z = np.load(f"spinal_run_{v}.npz", allow_pickle=True)
    names = [str(x) for x in z["neuro_names"]]
    joints = [str(x) for x in z["key_joints"]]
    print("npz neuro_names:", names)
    print("npz knee col idx:", joints.index("knee_angle_r"),
          "RG_E_r idx:", names.index("RG_E_r"),
          "cfg:", str(z["cfg"]))

print("===== protected artifacts =====")
for f in ("spinal_run.npz", "optuna_walk.db"):
    st = os.stat(f)
    print(f"{f}: size={st.st_size} mtime={st.st_mtime}")
studies = optuna.get_all_study_names(storage="sqlite:///optuna_walk.db")
print(f"optuna_walk.db studies: {len(studies)}")
print("variant study names leaked into optuna_walk.db:",
      [s for s in studies if s.startswith("curr_w2lvar")
       or s.startswith("curr_syn6")])
