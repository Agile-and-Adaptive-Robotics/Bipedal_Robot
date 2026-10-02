"""Generate the SADb knowledge base (markdown) from the corpus export (2026-09-28).

Output: SADb_audit/knowledge_base/
  INDEX.md                    corpus overview + navigation
  clusters/<NN>-<slug>.md     per topic-cluster: stats, top authors, key papers
                              (top-cited, with curation notes)
  pathways/<slug>.md          per feedback pathway: papers + their notes
  afferents/<type>.md         per afferent type: papers + notes
  animals/<animal>.md         per animal: papers + notes

Regenerate after every curation batch:
  myo python export_corpus.py && myo python build_knowledge_base.py
Stdlib only. Slugs are filesystem-safe.
"""
import json, os, re
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
EXPORT = os.path.join(HERE, "export")
KB = os.path.join(HERE, "knowledge_base")
records = json.load(open(os.path.join(EXPORT, "sadb_export.json"), encoding="utf-8"))
layout = json.load(open(os.path.join(EXPORT, "sadb_layout.json"), encoding="utf-8"))
labels = {int(k): v for k, v in json.load(
    open(os.path.join(EXPORT, "cluster_labels.json"), encoding="utf-8")).items()}

for d in ("clusters", "pathways", "afferents", "animals"):
    os.makedirs(os.path.join(KB, d), exist_ok=True)


def slug(s):
    s = re.sub(r"[^a-zA-Z0-9]+", "-", s).strip("-").lower()
    return s[:60] or "x"


def paper_line(r):
    out = f"- **{r['primary'] or '?'} {r['year'] or ''}** — [{r['title'][:110]}]" \
          f"(https://doi.org/{r['doi']})" if r.get("doi") else \
          f"- **{r['primary'] or '?'} {r['year'] or ''}** — {r['title'][:110]}"
    bits = []
    if r.get("animals"):
        bits.append("animals: " + ", ".join(r["animals"]))
    if r.get("afferents"):
        bits.append("afferents: " + ", ".join(r["afferents"]))
    if r.get("feedback"):
        bits.append("pathways: " + "; ".join(r["feedback"]))
    if bits:
        out += "  \n  - " + " · ".join(bits)
    if r.get("notes"):
        out += "  \n  - " + r["notes"][:600]
    if r.get("robot_sim"):
        out += "  \n  - *Robot/sim:* " + r["robot_sim"][:400]
    return out


byid = {r["id"]: r for r in records}
for r in records:
    L = layout.get(r["id"], {})
    r["_cl"] = L.get("cl", -1)
    r["_c"] = L.get("c", 0)

n_notes = sum(1 for r in records if r["has_notes"])
n_af = sum(1 for r in records if r.get("afferents"))
n_rs = sum(1 for r in records if r.get("robot_sim"))

# ---- INDEX ----
idx = ["# SADb Knowledge Base", "",
       f"Generated from the Airtable \u201cSensory Feedback\u201d corpus: "
       f"**{len(records)} papers**, {n_notes} with curation notes, "
       f"{n_af} with afferent classification, {n_rs} with robot/sim translation notes.", "",
       "## Topic clusters", ""]
for cl in sorted(labels):
    rs = [r for r in records if r["_cl"] == cl]
    if not rs:
        continue
    idx.append(f"- [{labels[cl]}](clusters/{cl:02d}-{slug(labels[cl])}.md) — {len(rs)} papers")
idx += ["", "## Feedback pathways", ""]
fb_count = Counter(f for r in records for f in r["feedback"])
for f, n in fb_count.most_common():
    idx.append(f"- [{f}](pathways/{slug(f)}.md) — {n} papers")
idx += ["", "## Afferent types", ""]
af_count = Counter(a for r in records for a in r.get("afferents", []))
for a, n in af_count.most_common():
    idx.append(f"- [{a}](afferents/{slug(a)}.md) — {n} papers")
idx += ["", "## Animals", ""]
an_count = Counter(a for r in records for a in r["animals"])
for a, n in an_count.most_common():
    idx.append(f"- [{a}](animals/{slug(a)}.md) — {n} papers")
open(os.path.join(KB, "INDEX.md"), "w", encoding="utf-8").write("\n".join(idx) + "\n")

# ---- cluster pages ----
for cl, label in labels.items():
    rs = sorted((r for r in records if r["_cl"] == cl), key=lambda r: -r["_c"])
    if not rs:
        continue
    yrs = [r["year"] for r in rs if r["year"]]
    authors = Counter(r["primary"] for r in rs if r["primary"])
    fbs = Counter(f for r in rs for f in r["feedback"])
    h = [f"# Cluster {cl + 1}: {label}", "",
         f"{len(rs)} papers · years {min(yrs) if yrs else '?'}–{max(yrs) if yrs else '?'} · "
         f"{sum(1 for r in rs if r['has_notes'])} curated", "",
         "## Top authors", ""]
    h += [f"- {a} ({n})" for a, n in authors.most_common(10)]
    if fbs:
        h += ["", "## Pathways represented", ""]
        h += [f"- {f} ({n})" for f, n in fbs.most_common(12)]
    h += ["", "## Key papers (by citations)", ""]
    h += [paper_line(r) for r in rs[:25]]
    open(os.path.join(KB, "clusters", f"{cl:02d}-{slug(label)}.md"), "w",
         encoding="utf-8").write("\n".join(h) + "\n")

# ---- pathway / afferent / animal pages ----
def group_pages(key, subdir, title_fn):
    groups = defaultdict(list)
    for r in records:
        for v in (r.get(key) or []):
            groups[v].append(r)
    for v, rs in groups.items():
        rs = sorted(rs, key=lambda r: -(r["year"] or 0))
        h = [f"# {title_fn(v)}", "", f"{len(rs)} papers in the corpus.", ""]
        h += [paper_line(r) for r in rs]
        open(os.path.join(KB, subdir, slug(v) + ".md"), "w",
             encoding="utf-8").write("\n".join(h) + "\n")


group_pages("feedback", "pathways", lambda v: f"Feedback pathway: {v}")
group_pages("afferents", "afferents", lambda v: f"Afferent type: {v}")
group_pages("animals", "animals", lambda v: f"Animal: {v}")

n_pages = 1 + len(labels) + len(fb_count) + len(af_count) + len(an_count)
print(f"knowledge base: {n_pages} pages in {KB} "
      f"({len(records)} papers, {n_notes} with notes)")
