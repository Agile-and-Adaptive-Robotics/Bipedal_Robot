"""Pre-launch state check for curr_syn6_s2 (stage-2 launch, syn6 chain).
Reports per-study trial states in optuna_syn6.db + confirms the stage-1
winner json that stage 2 chains from exists.
"""
import json
import os
import sqlite3

SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
db = os.path.join(SP, "optuna_syn6.db")
con = sqlite3.connect(db)
cur = con.cursor()
cur.execute("SELECT study_id, study_name FROM studies ORDER BY study_id")
for sid, name in cur.fetchall():
    cur.execute("SELECT state, COUNT(*) FROM trials WHERE study_id=? "
                "GROUP BY state", (sid,))
    print("study", sid, name, "->", cur.fetchall())
con.close()
s1 = os.path.join(SP, "curriculum_syn6_stage1.json")
print("stage1 json (chaining source) exists:", os.path.exists(s1))
if os.path.exists(s1):
    d = json.load(open(s1, encoding="utf-8"))
    print("stage1 winner: score", d["score"], "trial", d["trial"])
