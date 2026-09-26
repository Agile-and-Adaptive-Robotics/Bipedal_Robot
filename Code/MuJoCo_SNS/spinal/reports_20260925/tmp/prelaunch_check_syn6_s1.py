"""Pre-launch state check for curr_syn6_s1 (stage-1 relaunch, syn6 chain).
Reports: stage-1 json content, per-study trial states in optuna_syn6.db,
and whether any _curriculum_syn6 python process is already running.
"""
import json
import os
import sqlite3

SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
jp = os.path.join(SP, "curriculum_syn6_stage1.json")
print("stage1 json exists:", os.path.exists(jp))
if os.path.exists(jp):
    d = json.load(open(jp, encoding="utf-8"))
    print("json keys:", sorted(d.keys()))
    print("json stage/score/trial/study:", d.get("stage"), d.get("score"),
          d.get("trial"), d.get("study"))
    print("json params:", json.dumps(d.get("params")))

db = os.path.join(SP, "optuna_syn6.db")
con = sqlite3.connect(db)
cur = con.cursor()
cur.execute("SELECT study_id, study_name FROM studies ORDER BY study_id")
for sid, name in cur.fetchall():
    cur.execute("SELECT state, COUNT(*) FROM trials WHERE study_id=? "
                "GROUP BY state", (sid,))
    print("study", sid, name, "->", cur.fetchall())
con.close()
