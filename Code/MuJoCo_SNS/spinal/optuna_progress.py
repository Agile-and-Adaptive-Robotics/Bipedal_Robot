import sqlite3

con = sqlite3.connect("optuna_walk.db")
cur = con.cursor()
for study in ("ground_walk_v1", "ground_walk_v2"):
    sid = cur.execute("select study_id from studies where study_name=?",
                      (study,)).fetchone()
    if not sid:
        continue
    sid = sid[0]
    n = cur.execute(
        "select count(*) from trials where state='COMPLETE' and study_id=?",
        (sid,)).fetchone()[0]
    f = cur.execute(
        "select count(*) from trials where state='FAIL' and study_id=?",
        (sid,)).fetchone()[0]
    r = cur.execute(
        "select count(*) from trials where state='RUNNING' and study_id=?",
        (sid,)).fetchone()[0]
    rows = cur.execute(
        "select t.number, v.value from trials t "
        "join trial_values v on v.trial_id = t.trial_id "
        "where t.state='COMPLETE' and t.study_id=? order by t.number",
        (sid,)).fetchall()
    scores = [float(v) for _, v in rows]
    print(f"== {study}: complete={n} fail={f} running={r}")
    if scores:
        bi = scores.index(max(scores))
        print(f"   best {max(scores):.3f} (trial {rows[bi][0]})")
        print(f"   last 8: {[f'{s:.2f}' for s in scores[-8:]]}")
