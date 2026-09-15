"""Build Primary Author / Secondary Authors payloads for all 552 SADb papers.

Rule (Ben, 2026-09-15): Primary Author = ONE surname, or 'A and B' when the paper
has exactly two authors. NO 'et al.'. All remaining surnames -> Secondary Authors
(multi-select bubbles). Author sources: personal Zotero items (creators) for the
99 originals + 50 demo + 383 rest-import; the RR-digest shortlist strings for the
20 digest papers.

Outputs (SADb_audit/author_fix/):
  payloads_XX.json   chunked update payloads (upsert on DOI; id-based for DOI-less)
  changes_report.csv record-by-record: old display author -> new primary (+ secondaries)
  unmatched.txt      records with no author source found
"""
import csv, json, os, re, unicodedata

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

def norm_title(t):
    t = unicodedata.normalize("NFKD", (t or "")).lower()
    t = re.sub(r"[^a-z0-9]+", " ", t)
    return " ".join(t.split())

def surname_of(creator):
    if "lastName" in creator and creator.get("lastName"):
        return str(creator["lastName"]).strip()
    if "name" in creator and creator.get("name"):
        n = str(creator["name"]).strip()
        parts = n.split()
        return parts[-1] if parts else n
    return ""

PERSONAL_JSON = os.environ.get("TEMP", SAD) + os.sep + "personal_fresh_authors.json"
if not os.path.exists(PERSONAL_JSON):
    PERSONAL_JSON = os.path.join(SAD, "personal_SAD_items_full.json")
with open(PERSONAL_JSON, encoding="utf-8-sig") as f:
    personal = json.load(f)

def clean_surnames(surnames):
    """Dedupe (accent/apostrophe-insensitive), drop initial-garbage and
    name-split artifacts ('Ray' when 'Le Ray' present)."""
    seen = set()
    dedup = []
    for s in surnames:
        k = unicodedata.normalize("NFD", s).encode("ascii", "ignore").decode().lower()
        k = k.replace("'", "").replace("\u2019", "")
        if not k or k in seen:
            continue
        if len(s) <= 2 and not any(v in s.lower() for v in "aeiou"):
            continue
        seen.add(k)
        dedup.append(s)
    # drop X when another surname ends with ' X' (split artifact)
    final = []
    for i, s in enumerate(dedup):
        ki = unicodedata.normalize("NFD", s).encode("ascii", "ignore").decode().lower()
        artifact = False
        for j, other in enumerate(dedup):
            if i == j:
                continue
            ko = unicodedata.normalize("NFD", other).encode("ascii", "ignore").decode().lower()
            if len(ko) > len(ki) and ko.endswith(" " + ki):
                artifact = True
                break
        if not artifact:
            final.append(s)
    return final

by_doi_cand, by_key_cand = {}, {}
for it in personal:
    surnames = [s for s in (it.get("surnames") or []) if s]
    if not surnames:
        continue
    doi = norm_doi(it.get("doi"))
    if doi:
        by_doi_cand.setdefault(doi, []).append(surnames)
    if it.get("key"):
        by_key_cand.setdefault(it["key"], []).append(surnames)

by_doi = {d: max(cands, key=len) for d, cands in by_doi_cand.items()}
by_key = {k: max(cands, key=len) for k, cands in by_key_cand.items()}

with open(os.path.join(SAD, "batch3", "rr_digest_shortlist.json"), encoding="utf-8") as f:
    shortlist = json.load(f)
digest_by_doi = {}
for s in shortlist:
    names = [p.strip().split(" ")[0] for p in s["authors"].split(";") if p.strip()]
    digest_by_doi[norm_doi(s["doi"])] = names
digest_by_doi[norm_doi("10.7554/elife.107480")] = ["Shevtsova", "Lockhart", "Rybak"]

targets = {}  # key -> (display_author, source)
DEMO_KEY_TO_ID = {}
with open(os.path.join(SAD, "airtable_created50_ids.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        DEMO_KEY_TO_ID[r["zotero_key"].strip()] = r["airtable_id"].strip()
with open(os.path.join(SAD, "airtable_papers_slim.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        raw_doi = (r.get("DOI") or "").strip()
        d = norm_doi(raw_doi)
        if d:
            targets[d] = (r.get("author", ""), "orig99", raw_doi)
with open(os.path.join(SAD, "airtable_created50_ids.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        targets["__key__" + r["zotero_key"].strip()] = (None, "demo50", "")
with open(os.path.join(SAD, "airtable_rest_import_clean.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        d = norm_doi(r.get("DOI"))
        k = r["zotero_key"].strip()
        targets[d if d else "__key__" + k] = (r.get("author", ""), "rest383", "")
for d in digest_by_doi:
        targets[d] = (None, "digest20", "")

records, unmatched = [], []
for key, (disp, source, raw_doi) in targets.items():
    surnames = None
    if source == "digest20":
        surnames = digest_by_doi.get(key)
    elif key.startswith("__key__"):
        surnames = by_key.get(key[7:])
    else:
        surnames = by_doi.get(key)
    if surnames:
        surnames = clean_surnames(surnames)
    if not surnames:
        unmatched.append((source, key))
        continue
    old = (disp or "").strip()
    old_n = re.sub(r"[^a-z]", "", old.lower())
    new_primary = surnames[0]
    if len(surnames) == 2:
        new_primary = surnames[0] + " and " + surnames[1]
    # guardrail: don't lose authors — if current display shows two authors
    # ("A and B") but the source only has one surname, keep the record out of
    # the primary-fix batch and log it.
    if old and " and " in old.lower() and len(surnames) < 2:
        unmatched.append(("skip-guard-2auth", key + " | old=" + old + " new=" + new_primary))
        continue
    new_n = re.sub(r"[^a-z]", "", new_primary.lower())
    primary_changed = bool(old) and old_n != new_n
    if len(surnames) > 2:
        secondary = surnames[1:]
    else:
        secondary = []
    rid = DEMO_KEY_TO_ID.get(key[7:]) if key.startswith("__key__") else None
    records.append({"match": key, "doi": "" if (key.startswith("__key__") or not key) else key,
                    "raw_doi": raw_doi if source == "orig99" else "",
                    "id": rid or "", "primary": new_primary, "secondary": secondary,
                    "old": disp, "source": source, "changed": primary_changed})

# manual overrides for records whose Zotero creators are mangled (Sept-11 imports)
OVERRIDES = {
    "recEfpqUr9aBwNnv4": ("Shik", ["Severin", "Orlovsky"]),  # Shik/Severin/Orlovsky 1966
}
for r in records:
    if r["id"] in OVERRIDES:
        r["primary"], r["secondary"] = OVERRIDES[r["id"]]
        r["changed"] = True

os.makedirs(HERE, exist_ok=True)
with open(os.path.join(HERE, "changes_report.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["match(doi_or_key)", "source", "old_display_author", "new_primary", "secondary_authors", "primary_changed"])
    for r in records:
        old = (r["old"] or "").strip()
        changed = "YES" if (old and old.lower() != r["primary"].lower()) else ("(was blank)" if not old else "no")
        w.writerow([r["match"], r["source"], old, r["primary"], "; ".join(r["secondary"]), changed])

with open(os.path.join(HERE, "unmatched.txt"), "w", encoding="utf-8") as f:
    for s, k in unmatched:
        f.write(f"{s}\t{k}\n")

KNOWN_IDS = {
    "YJVT4HJU": "recQwvJAHTH2EEnnZ",
    "YX6IC7SQ": "recY6tB6S2au4jFCY",
    "QRNCP58H": "recmmZTBWWeGV3Fvn",
}
payload = []
for r in records:
    if r["primary"] == "Ml":  # Shik M.L. mangled creator (any duplicate entry)
        r["primary"], r["secondary"] = "Shik", ["Severin", "Orlovsky"]
        r["changed"] = True
    if r["source"] == "orig99" and not r["changed"]:
        continue  # chunk-00 attempt 1 already wrote these (primary + secondaries)
    if r["id"]:
        payload.append({"id": r["id"], "fields": {"Primary Author": r["primary"], "Secondary Authors": r["secondary"]}})
    elif r["doi"]:
        use = r.get("raw_doi") or r["doi"]
        payload.append({"fields": {"DOI": use, "Primary Author": r["primary"], "Secondary Authors": r["secondary"]}})
    else:
        k = r["match"].replace("__key__", "")
        rid = KNOWN_IDS.get(k)
        if rid:
            payload.append({"id": rid, "fields": {"Primary Author": r["primary"], "Secondary Authors": r["secondary"]}})
        else:
            unmatched.append(("no-id-and-no-doi", r["match"]))

# ---------- DOI-less 'et al.' stragglers: fix by id; full data if title matches ----------
by_title_fresh = {}
for it in personal:
    t = norm_title(it.get("title"))
    if t and it.get("surnames"):
        by_title_fresh.setdefault(t, it["surnames"])
with open(os.path.join(HERE, "doi_less_etal.json"), encoding="utf-8") as f:
    stragglers = json.load(f)
for s in stragglers:
    sn = by_title_fresh.get(norm_title(s["title"]))
    if sn:
        sn = clean_surnames(sn)
    if sn and len(sn) >= 1:
        if len(sn) == 1:
            primary, secondary = sn[0], []
        elif len(sn) == 2:
            primary, secondary = sn[0] + " and " + sn[1], []
        else:
            primary, secondary = sn[0], sn[1:]
    else:
        primary = re.sub(r"\s*et al\.?\s*$", "", s["display"]).strip()
        secondary = []
    payload.append({"id": s["id"], "fields": {"Primary Author": primary, "Secondary Authors": secondary}})

chunks = [payload[i:i + 50] for i in range(0, len(payload), 50)]
for i, ch in enumerate(chunks):
    with open(os.path.join(HERE, f"payloads_{i:02d}.json"), "w", encoding="utf-8") as f:
        json.dump({"records": ch}, f, ensure_ascii=False)

print("records resolved:", len(records))
print("unmatched:", len(unmatched))
for u in unmatched[:25]:
    print("  UNMATCHED:", u)
print("chunks:", len(chunks))
changed = sum(1 for r in records if r["old"] and r["old"].strip().lower() != r["primary"].lower())
blank_old = sum(1 for r in records if not (r["old"] or "").strip())
print("primary changed:", changed, "| old blank:", blank_old)
