"""Consolidate curation flags across batches into curation_flags_20260928.md (for Ben)."""
import glob, json, os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
by_rec = defaultdict(list)
for p in sorted(glob.glob(os.path.join(HERE, "curation_out", "batch_*.json"))):
    slug = os.path.basename(p)[6:8]
    for e in json.load(open(p, encoding="utf-8-sig")):
        fl = (e.get("flags") or "").strip()
        if fl:
            by_rec[(e["id"], slug)].append(fl)

lines = ["# Curation flags for Ben (2026-09-28 v2 campaign)", "",
         f"{len(by_rec)} papers carry flags from the subagent curators.",
         "Anything proposing NEW Feedback vocabulary needs Ben's naming decision",
         "before any Airtable change.", ""]
for (rid, slug), fls in sorted(by_rec.items()):
    lines.append(f"- `{rid}` (batch {slug}): " + " | ".join(fls))
open(os.path.join(HERE, "curation_flags_20260928.md"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
print(f"{len(by_rec)} flagged papers -> curation_flags_20260928.md")
