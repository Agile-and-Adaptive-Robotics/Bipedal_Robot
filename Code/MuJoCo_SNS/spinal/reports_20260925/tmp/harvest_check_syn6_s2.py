"""Harvest state check for curr_syn6_s2_air_aff (stage 2, syn6 chain)."""
import json
import os
import sqlite3

SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
db = os.path.join(SP, "optuna_syn6.db")
con = sqlite3.connect(db)
cur = con.cursor()
cur.execute("SELECT study_id FROM studies WHERE study_name=?",
            ("curr_syn6_s2_air_aff",))
sid = cur.fetchone()[0]
cur.execute("SELECT state, COUNT(*) FROM trials WHERE study_id=? "
            "GROUP BY state", (sid,))
print("trial states:", cur.fetchall())
cur.execute("SELECT t.number, t.state, tv.value FROM trials t "
            "LEFT JOIN trial_values tv ON tv.trial_id = t.trial_id "
            "AND tv.objective = 0 WHERE t.study_id=? ORDER BY t.number",
            (sid,))
vals = cur.fetchall()
for n, st, v in vals:
    print(f"trial {n:3d} {st:9s} value={v}")
c_nan = c_gate = c_real = 0
for n, st, v in vals:
    if v is None:
        continue
    if abs(v + 200.0) < 1e-9:
        c_nan += 1
    elif v <= 25.0:
        c_gate += 1
    else:
        c_real += 1
print(f"floor census: -200={c_nan} <=25(gate/static)={c_gate} "
      f">25(real rhythm)={c_real}")
con.close()
j = json.load(open(os.path.join(SP, "curriculum_syn6_stage2.json"),
                   encoding="utf-8"))
print("winner json: score", j["score"], "trial", j["trial"],
      "params", json.dumps(j["params"]))
