"""Re-check the 177 'no-PDF' papers against: (a) staged files by TITLE match,
(b) Zotero web API server-side files (personal + AARL group, via key),
(c) their Airtable record's Attachments field (via PAT).
Output: author_fix/paperhunt_status.csv
Keys via env: AT_PAT, ZOT_KEY. Not written to disk.
"""
import csv, os, json, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)
AT_PAT = os.environ["AT_PAT"]
ZK = os.environ["ZOT_KEY"]
BASE = "appMQTnobUNRytIp7"
UID = "631450"
GID = "735051"

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

def norm_title(t):
    import re, unicodedata
    t = unicodedata.normalize("NFKD", (t or "")).lower()
    t = t.replace("ı", "i")
    t = re.sub(r"[^a-z0-9]+", " ", t)
    return " ".join(t.split())

def get(url, headers=None):
    req = urllib.request.Request(url, headers=headers or {})
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.loads(r.read().decode())

# 1) load the 177 + their titles from corpus sources
missing = []
with open(os.path.join(HERE, "pdfs_missing_no_zotero.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        missing.append(r)
titles = {}
for f, col in ((os.path.join(SAD, "airtable_papers_slim.csv"), None),
               (os.path.join(SAD, "airtable_rest_import_clean.csv"), None)):
    with open(f, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            d = norm_doi(r.get("DOI"))
            if d:
                titles[d] = r.get("title", "")
with open(os.path.join(SAD, "batch3", "rr_digest_shortlist.json"), encoding="utf-8") as f:
    for s in json.load(f):
        titles[norm_doi(s["doi"])] = s["title"]

# 2) staged-by-title index
staged_by_title = {}
with open(os.path.join(SAD, "staged_with_folders.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        t = norm_title(r["title"])
        if t:
            staged_by_title.setdefault(t, []).append(r["att_key"])

def zot_children_file(scope, item_key):
    """Check server-side PDF file for an item via web API. Returns True/False/None."""
    try:
        kids = get(f"https://api.zotero.org/{scope}/items/{item_key}/children?format=json&key=" + ZK)
        for k in kids:
            if k["data"].get("contentType") == "application/pdf":
                att = k["key"]
                url = f"https://api.zotero.org/{scope}/items/{att}/file"
                req = urllib.request.Request(url, headers={"Authorization": "Bearer " + ZK}, method="GET")
                try:
                    with urllib.request.urlopen(req, timeout=30) as r:
                        head = r.read(8)
                    return head.startswith(b"%PDF")
                except urllib.error.HTTPError as e:
                    return False if e.code in (404, 403) else None
                except Exception:
                    return None
    except Exception:
        return None
    return False

def zot_search(scope, doi, title):
    """Find item for doi/title in a scope; returns (item_key, has_server_pdf)."""
    q = urllib.parse.quote(doi) if doi else urllib.parse.quote(title[:60])
    try:
        res = get(f"https://api.zotero.org/{scope}/items?format=json&key={ZK}&q={q}&itemType=-attachment&limit=5")
    except Exception:
        return None, None
    want_title = norm_title(title)
    for it in res:
        d = it["data"]
        if doi and norm_doi(d.get("DOI")) == doi:
            return it["key"], zot_children_file(scope, it["key"])
        if not doi and norm_title(d.get("title")) == want_title:
            return it["key"], zot_children_file(scope, it["key"])
    return None, False

rows = []
n = len(missing)
for i, m in enumerate(missing):
    d = norm_doi(m["doi"])
    title = titles.get(d, "")
    # (a) staging by title
    st = staged_by_title.get(norm_title(title), [])
    # (b) zotero web: personal then group
    uk, uf = zot_search("users/" + UID, d, title)
    time.sleep(0.1)
    gk, gf = zot_search("groups/" + GID, d, title)
    time.sleep(0.1)
    rows.append({"doi": d, "title": title[:80], "staged_title_match": ";".join(st[:2]),
                 "zot_user_key": uk or "", "zot_user_pdf": uf,
                 "zot_group_key": gk or "", "zot_group_pdf": gf})
    if (i + 1) % 25 == 0:
        print(f"web-checked {i+1}/{n}")

# (c) Airtable attachments for these 177 records
recs = {}
offset = ""
hdr = {"Authorization": "Bearer " + AT_PAT}
while True:
    q = BASE + "/Papers?pageSize=100&fields%5B%5D=DOI&fields%5B%5D=Attachments"
    if offset:
        q += "&offset=" + offset
    req = urllib.request.Request("https://api.airtable.com/v0/" + q, headers=hdr)
    with urllib.request.urlopen(req, timeout=40) as r:
        page = json.loads(r.read().decode())
    for r in page.get("records", []):
        recs[norm_doi((r.get("fields") or {}).get("DOI"))] = (r["id"], (r.get("fields") or {}).get("Attachments") or [])
    offset = page.get("offset")
    if not offset:
        break

final = []
for r in rows:
    rid, atts = recs.get(r["doi"], ("", []))
    r["airtable_id"] = rid
    r["airtable_has_file"] = "yes" if atts else "no"
    final.append(r)

with open(os.path.join(HERE, "paperhunt_status.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(final[0].keys()))
    w.writeheader()
    w.writerows(final)

n_st = sum(1 for r in final if r["staged_title_match"])
n_u = sum(1 for r in final if r["zot_user_pdf"] is True)
n_g = sum(1 for r in final if r["zot_group_pdf"] is True)
n_a = sum(1 for r in final if r["airtable_has_file"] == "yes")
print(f"177 recheck: staged-title {n_st} | zotero-user-pdf {n_u} | zotero-group-pdf {n_g} | airtable-already {n_a} | truly missing: {sum(1 for r in final if not (r['staged_title_match'] or r['zot_user_pdf'] or r['zot_group_pdf'] or r['airtable_has_file']=='yes'))}")
print("wrote paperhunt_status.csv")
