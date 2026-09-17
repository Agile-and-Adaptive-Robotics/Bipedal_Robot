"""Match Ben's found PDFs (D:\sadb_pdf_staging\pdfs) to Downloaded=x rows in
paperhunt_status.csv via Crossref surname+year, print proposed pairs.
"""
import csv, json, os, time, urllib.parse, urllib.request, re

HERE = os.path.dirname(os.path.abspath(__file__))
PDFS = r"D:\sadb_pdf_staging\pdfs"

def crossref(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(url, headers={"User-Agent": "sadb-match/0.1 (mailto:benjamin.bolen@pdx.edu)"})
    with urllib.request.urlopen(req, timeout=30) as r:
        m = json.loads(r.read().decode())["message"]
    fam = m.get("author", [{}])[0].get("family", "?") if m.get("author") else "?"
    year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
    return str(fam), str(year)

rows = []
with open(os.path.join(HERE, "paperhunt_status.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        rows.append(r)
marked = [r for r in rows if (r.get("Downloaded") or "").strip().lower() == "x"]
print("marked x:", len(marked))

files = [f for f in os.listdir(PDFS) if f.lower().endswith(".pdf")]

# parse filename -> (surname_token, year)
def parse_fname(f):
    stem = f[:-4]
    m = re.match(r"([A-Za-z\-]+)[_\- ]?(\d{4})", stem)
    if m:
        return m.group(1).lower(), m.group(2)
    m = re.search(r"(\d{4})", stem)
    if m:
        # "The Journal of Physiology - 2025 - Singh - ..." style
        parts = re.split(r"[-–]", stem)
        for p in parts:
            p = p.strip()
            if p and not p.isdigit() and len(p) > 2:
                return p.split()[0].lower(), m.group(1)
    return stem.lower(), ""

fidx = {}
for f in files:
    sn, yr = parse_fname(f)
    fidx.setdefault((sn, yr), []).append(f)

pairs, unmatched_rows, unmatched_files = [], [], list(files)
for r in marked:
    doi = r["doi"]
    try:
        fam, yr = crossref(doi)
    except Exception as e:
        fam, yr = "?", "?"
    key = (fam.lower(), yr)
    hit = fidx.get(key)
    title = r["title"][:55]
    if hit:
        f = hit[0]
        unmatched_files = [x for x in unmatched_files if x != f]
        pairs.append((doi, title, fam + " " + yr, f))
        print("MATCH  %-22s %-28s <- %s" % (fam + " " + yr, title, f))
    else:
        unmatched_rows.append((doi, title, fam + " " + yr))
        print("NOMATCH %-22s %-28s" % (fam + " " + yr, title))
    time.sleep(0.3)

print()
print("unmatched rows (%d):" % len(unmatched_rows))
for u in unmatched_rows:
    print("   ", u[2], "|", u[1])
print("unmatched files (%d):" % len(unmatched_files))
for u in unmatched_files:
    print("   ", u)

with open(os.path.join(HERE, "pdf_match_pairs.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["doi", "title", "crossref_author_year", "pdf_file"])
    for p in pairs:
        w.writerow(p)
print("wrote pdf_match_pairs.csv (%d pairs)" % len(pairs))
