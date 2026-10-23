# Build the scratch config for the config-change proof: seeded config + ONE
# author ("Molkov", Rybak-lab, in the corpus at 3 retained papers so NOT in the
# seeded 20) + ONE keyword ("sensory feedback", corpus term deliberately left
# out of the seed). Written under tmp\ — the real config is untouched.
import json, shutil

SRC = r"D:\Github\Bipedal_Robot\SADb_audit\watch\watch_config.json"
DST = r"D:\Github\Bipedal_Robot\tmp\watch_config_scratch.json"

cfg = json.load(open(SRC, encoding="utf-8"))
assert "Molkov" not in cfg["authors"], "Molkov already seeded?"
assert "sensory feedback" not in cfg["keywords"], "sensory feedback already seeded?"
cfg["authors"] = sorted(cfg["authors"] + ["Molkov"])
cfg["keywords"] = cfg["keywords"] + ["sensory feedback"]
json.dump(cfg, open(DST, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
print("wrote", DST, "| authors:", len(cfg["authors"]), "(+Molkov) | keywords:", len(cfg["keywords"]), "(+sensory feedback)")
