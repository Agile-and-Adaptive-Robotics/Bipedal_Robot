import json
import os

ROOT = r"D:\Github\Bipedal_Robot\SADb_audit"

# What record types does the export/app payload carry?
d = json.load(open(os.path.join(ROOT, "export", "sadb_export.json"), encoding="utf-8"))
print("export type:", type(d).__name__, "records:", len(d))
first = d[0]
print("paper-record fields:", sorted(first.keys()))

# does any non-paper table appear anywhere in the payload?
blob = open(os.path.join(ROOT, "export", "sadb_export.json"), encoding="utf-8").read()
for kw in ('"Feedback"', '"Models"', '"Review Papers"', '"Recorders"', '"pathway"', '"reviews"'):
    print(f"payload contains {kw!r}:", kw in blob)

# knowledge base: what page types does the app generate?
kb = os.path.join(ROOT, "knowledge_base")
if os.path.isdir(kb):
    subs = os.listdir(kb)
    print("knowledge_base entries:", len(subs), subs[:8])

# what does the app html embed beyond papers (tab names)?
html = open(os.path.join(ROOT, "app", "sadb_app.html"), encoding="utf-8", errors="replace").read()
for kw in ("Papers", "Feedback", "Models", "Review", "Recorder", "sadb-data", "CITES"):
    i = html.find(kw)
    print(f"app html has {kw!r}:", i >= 0)
