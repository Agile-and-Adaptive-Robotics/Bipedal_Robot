# Supplementary author-axis check: seeded config + "Audet" + "sensory feedback".
# Audet (Frigon group) has a known in-window Crossref paper (10.1152/jn.00240.2026),
# so the author trigger must produce visible new rows.
import json

SRC = r"D:\Github\Bipedal_Robot\SADb_audit\watch\watch_config.json"
DST = r"D:\Github\Bipedal_Robot\tmp\watch_config_scratch2.json"

cfg = json.load(open(SRC, encoding="utf-8"))
assert "Audet" not in cfg["authors"] and "sensory feedback" not in cfg["keywords"]
cfg["authors"] = sorted(cfg["authors"] + ["Audet"])
cfg["keywords"] = cfg["keywords"] + ["sensory feedback"]
json.dump(cfg, open(DST, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
print("wrote", DST, "| authors:", len(cfg["authors"]), "(+Audet) | keywords:", len(cfg["keywords"]), "(+sensory feedback)")
