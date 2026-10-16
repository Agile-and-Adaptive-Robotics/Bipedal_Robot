"""Final BPA winner summaries."""
import json

for f in ("curriculum_bpa_ns.json", "curriculum_bpa_spk.json"):
    d = json.load(open(f, encoding="utf-8"))
    print(f, "score", round(d["score"], 2), "trial", d["trial"])
    print("  ", {k: (round(v, 3) if isinstance(v, float) else v)
                 for k, v in d["params"].items()})
