import csv, sys
p = r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_export.csv"
rows = list(csv.DictReader(open(p, encoding="utf-8")))
total = len(rows)
arch = [r for r in rows if r["archived"].strip().lower() in ("true", "1", "yes")]
kept = [r for r in rows if r not in arch]
def n(rs, col):
    return sum(1 for r in rs if r[col] and r[col].strip())
print("export total:", total)
print("archived:", len(arch))
print("kept:", len(kept))
for label, rs in (("kept", kept), ("all", rows)):
    print(f"[{label}] notes={n(rs,'has_notes') if False else sum(1 for r in rs if r['has_notes'].strip().lower() in ('true','1','yes'))} animals={n(rs,'animals')} feedback={n(rs,'feedback')} afferents={n(rs,'afferents')} robot_sim={n(rs,'robot_sim')} animal_study={n(rs,'animal_study')} doi={n(rs,'doi')} pdf={sum(1 for r in rs if r['has_pdf'].strip().lower() in ('true','1','yes'))}")
