import csv
rows = list(csv.DictReader(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_export.csv",
                                encoding="utf-8-sig")))
ret = [r for r in rows if r["prune"] != "True" and r["archived"] != "True"]
T = lambda s: [t.strip() for t in s.replace(";", ",").split(",") if t.strip()]
ARTH = {"Insects", "Stick Insect", "Cockroach", "Arthropods"}
for r in ret:
    if T(r["afferents"]) and (set(T(r["animals"])) & ARTH):
        print(T(r["afferents"]), "|", T(r["animals"]), "|", r["primary"][:40])
