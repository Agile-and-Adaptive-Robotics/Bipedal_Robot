import csv
rows = list(csv.DictReader(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_export.csv",
                                encoding="utf-8-sig")))
ret = [r for r in rows if r["prune"] != "True" and r["archived"] != "True"]
T = lambda s: [t.strip() for t in s.replace(";", ",").split(",") if t.strip()]
AFF = {"Ia": ["Ia"], "Ib": ["Ib"], "II": ["II"], "Cut": ["Cutaneous"],
       "Other": ["Mechanosensory", "Flexor reflex afferents", "III/IV", "Nociceptive"]}
ANI = {"Cat": ["Cat"], "Human": ["Human"], "Rod": ["Mice", "Rat"],
       "Arth": ["Insects", "Stick Insect", "Cockroach", "Arthropods"],
       "LowV": ["Lamprey", "Salamander", "Frog", "Zebrafish"],
       "Other": ["Mammals", "Vertebrates", "Dog"], "NoTag": []}
g = {a: {n: 0 for n in ANI} for a in AFF}
for r in ret:
    at, an = T(r["afferents"]), T(r["animals"])
    if not at:
        continue
    for a, akeys in AFF.items():
        if set(at) & set(akeys):
            for n, nkeys in ANI.items():
                if (set(nkeys) & set(an)) if nkeys else (not an):
                    g[a][n] += 1
for a in g:
    print(a.ljust(5), " ".join(f"{n}={g[a][n]:<3}" for n in g[a]))
