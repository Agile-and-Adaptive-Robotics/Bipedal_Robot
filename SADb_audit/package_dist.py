"""Package the advisor-ready SADb Explorer distribution (2026-09-28).

Assembles SADb_audit/dist/ :
  sadb_app.html          the app (copied from app/)
  README_START_HERE.md   advisor front door (patches N_PAPERS/DATE placeholders)
  CURATION_NOTES.md      curation state + open items (patched counts)
  sadb_export.csv        flat corpus
  knowledge_base/        markdown KB (copied)
then zips it to dist/SADb_Explorer_<date>.zip. Run AFTER the full rebuild
chain (export -> citation graph -> knowledge base -> app). Stdlib only.
"""
import datetime, json, os, shutil, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DIST = os.path.join(HERE, "dist")
os.makedirs(DIST, exist_ok=True)
records = json.load(open(os.path.join(HERE, "export", "sadb_export.json"), encoding="utf-8"))
n = len(records)
n_notes = sum(1 for r in records if r["has_notes"])
n_af = sum(1 for r in records if r.get("afferents"))
today = datetime.date.today().isoformat()

# fresh copies
shutil.copy2(os.path.join(HERE, "app", "sadb_app.html"), os.path.join(DIST, "sadb_app.html"))
shutil.copy2(os.path.join(HERE, "export", "sadb_export.csv"), os.path.join(DIST, "sadb_export.csv"))
kb_src = os.path.join(HERE, "knowledge_base")
kb_dst = os.path.join(DIST, "knowledge_base")
if os.path.isdir(kb_dst):
    shutil.rmtree(kb_dst)
if os.path.isdir(kb_src):
    shutil.copytree(kb_src, kb_dst)

readme = open(os.path.join(DIST, "README_START_HERE.md"), encoding="utf-8").read()
readme = readme.replace("N_PAPERS", str(n)).replace("DATE", today)
open(os.path.join(DIST, "README_START_HERE.md"), "w", encoding="utf-8").write(readme)

notes = f"""# Curation state (generated {today})

- {n} papers in the corpus; {n_notes} carry distilled curation notes.
- {n_af} papers have an afferent-type classification (Ia / Ib / II / III-IV /
  mechanosensory / cutaneous / heat / nociceptive / flexor-reflex afferents /
  legacy "type 1" terminology).
- Feedback pathways, animals, models, and review coverage are linked per paper
  in the Airtable base and reflected here.
- Papers whose abstract could not be located (no DOI, or no open abstract)
  are intentionally left without notes rather than guessed.

## Open items for the lab
- ~100 papers still lack notes (no findable open abstract; a library pass
  with institutional access would close most of these).
- Prune shortlist is flagged in Airtable (the base is over the free-plan
  1000-record cap); deletions await Ben's approval.
- Proposed new Feedback vocabulary entries (from curation flags) are logged
  in SADb_audit/curation_log.csv for review before any Airtable change.
"""
open(os.path.join(DIST, "CURATION_NOTES.md"), "w", encoding="utf-8").write(notes)

zip_path = os.path.join(DIST, f"SADb_Explorer_{today}.zip")
if os.path.exists(zip_path):
    os.remove(zip_path)
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
    for root, _, files in os.walk(DIST):
        for f in files:
            full = os.path.join(root, f)
            if f.endswith(".zip"):
                continue
            z.write(full, os.path.relpath(full, DIST))
print(f"packaged {zip_path} ({os.path.getsize(zip_path)/1e6:.1f} MB) — "
      f"{n} papers, {n_notes} curated, {n_af} afferent-classified")
