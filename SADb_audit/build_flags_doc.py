"""Regenerate curation_flags_20260928.md with hyperlinks (2026-10-01, LOCAL ONLY).

Every record ID becomes a clickable link straight to the Airtable record
(https://airtable.com/<base>/<table>/<recordId>) so Ben never has to search
for recXXX ids. Each entry also carries Author Year + short title + DOI link,
pulled from the cached export (no API calls — the workspace is over its
monthly API limit). Flags grouped by auto-tagged category.
"""
import glob, json, os, re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "appMQTnobUNRytIp7"
TABLE = "tblnnMrZszhboU4uD"
rec_url = lambda rid: f"https://airtable.com/{BASE}/{TABLE}/{rid}"

records = {r["id"]: r for r in json.load(
    open(os.path.join(HERE, "export", "sadb_export.json"), encoding="utf-8"))}

entries = []
for p in sorted(glob.glob(os.path.join(HERE, "curation_out", "batch_*.json"))):
    slug = os.path.basename(p)[6:8]
    for e in json.load(open(p, encoding="utf-8-sig")):
        fl = (e.get("flags") or "").strip()
        if not fl:
            continue
        r = records.get(e["id"], {})
        entries.append({
            "id": e["id"], "batch": slug, "flag": fl,
            "label": f"{r.get('primary') or '?'} {r.get('year') or ''}".strip(),
            "title": (r.get("title") or "")[:90], "doi": r.get("doi") or "",
        })


def category(fl):
    f = fl.lower()
    if any(k in f for k in ("duplicate", "abstract swap", "metadata", "mismatch",
                            "identical to", "wrong abstract")):
        return "Data concerns (possible duplicate / wrong metadata)"
    if any(k in f for k in ("candidate", "no exact feedback vocabulary", "fits no",
                            "no vocabulary entry", "vocabulary term",
                            "not fit", "left empty", "absent from")):
        return "New-pathway / vocabulary proposals"
    if any(k in f for k in ("animals vocabulary", "consider adding", "not in the animals",
                            "mapped to", "closest allowed")):
        return "Animal-vocabulary gaps"
    return "Judgment calls (tagging conservatism etc.)"


groups = defaultdict(list)
for e in entries:
    groups[category(e["flag"])].append(e)

ORDER = ["Data concerns (possible duplicate / wrong metadata)",
         "New-pathway / vocabulary proposals",
         "Animal-vocabulary gaps",
         "Judgment calls (tagging conservatism etc.)"]

lines = ["# Curation flags for Ben (2026-09-28 v2 campaign — hyperlinked edition)",
         "",
         f"{len(entries)} flagged papers. Every paper title below is a clickable link",
         "straight to its Airtable record (opens in your browser, logged in). DOI links",
         "open the publisher page. Categories are auto-tagged from the flag wording —",
         "skimming 'Data concerns' first is worthwhile. Regenerate after future batches:",
         "`myo python build_flags_doc.py`. No Airtable API calls are used.", ""]

for cat in ORDER:
    es = sorted(groups.get(cat, []), key=lambda e: (e["label"], e["batch"]))
    if not es:
        continue
    lines += [f"## {cat} ({len(es)})", ""]
    for e in es:
        doi = f" · [doi]({'https://doi.org/' + e['doi']})" if e["doi"] else ""
        lines.append(f"- [{e['label']} — {e['title']}]({rec_url(e['id'])})"
                     f"{doi}  \n  {e['flag']}")
    lines.append("")

out = os.path.join(HERE, "curation_flags_20260928.md")
open(out, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print(f"wrote {out}: {len(entries)} flags across {sum(1 for c in ORDER if groups.get(c))} categories")

# also dump a machine-readable copy (with category) for the docx builder
for e in entries:
    e["cat"] = category(e["flag"])
json.dump(entries, open(os.path.join(HERE, "curation_out", "_flags_flat.json"), "w"),
          ensure_ascii=False, indent=1)
