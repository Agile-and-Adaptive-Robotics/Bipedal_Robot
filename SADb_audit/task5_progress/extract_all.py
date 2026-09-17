"""Task 5 main pipeline: match all 374 gap papers to staged PDFs, download any
missing from AARL Zotero web API, extract text, and write batch JSONs.
Env: ZOT_KEY. Run with anaconda base python + pypdf from D:\pdfx_tmp.
"""
import csv, json, os, re, subprocess, sys, time, urllib.parse, urllib.request

SAD = r"D:\Github\Bipedal_Robot\SADb_audit"
STAGING = r"D:\sadb_pdf_staging"
PROG = os.path.join(SAD, "task5_progress")
os.makedirs(PROG, exist_ok=True)
ZKEY = os.environ.get("ZOT_KEY", "")
sys.path.insert(0, r"D:\pdfx_tmp")
import pypdf

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

def norm_title(t):
    import unicodedata
    t = unicodedata.normalize("NFKD", (t or "")).lower()
    t = t.replace("ı", "i")
    t = re.sub(r"[^a-z0-9]+", " ", t)
    return " ".join(t.split())

def extract_pdf(path, max_chars=3500):
    try:
        r = pypdf.PdfReader(path)
        chunks = []
        n = 0
        for pg in r.pages[:3]:
            t = pg.extract_text() or ""
            chunks.append(t)
            n += len(t)
            if n > max_chars:
                break
        txt = " ".join(chunks)
        return re.sub(r"\s+", " ", txt).strip()[:max_chars]
    except Exception as e:
        return "EXTRACT-FAIL: " + str(e)[:100]

# 1) load all 374 in-scope rows
allrows = []
with open(os.path.join(SAD, "author_fix", "aarl_gap_by_folder.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if not r["folder_bucket"].startswith("OUT-OF-SCOPE"):
            allrows.append(r)
print("in-scope rows:", len(allrows))

# 2) staging index: att_key -> path; doi -> path; surname_year -> path
by_key, by_doi, by_sur_yr = {}, {}, {}
for f in os.listdir(STAGING):
    if not f.endswith(".pdf") or "__" not in f:
        continue
    path = os.path.join(STAGING, f)
    key = f.rsplit("__", 1)[1][:-4]
    by_key[key] = path
    # surname_year from filename prefix
    m = re.match(r"(.+)_(\d{4})__", f)
    if m:
        by_sur_yr.setdefault((m.group(1).lower(), m.group(2)), []).append(path)

# 3) inventory doi -> att_keys with files
inv_doi_keys = {}
with open(os.path.join(SAD, "pdf_inventory.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r["file_exists"].lower() == "true" and r["doi"]:
            k = r["att_key"]
            p = os.path.join(STAGING, "")
            # find the actual staged file for this att_key
            if k in by_key:
                inv_doi_keys.setdefault(norm_doi(r["doi"]), []).append(by_key[k])

# 4) match rows -> staged path (by att_key, then doi, then surname_year)
matched, need_download = 0, 0
for r in allrows:
    k = r["att_key"]
    path = by_key.get(k)
    if not path:
        d = norm_doi(r["doi"])
        cands = inv_doi_keys.get(d, [])
        if cands:
            path = cands[0]
    if not path:
        yr = ""
        if re.search(r"\d{4}", r["year"]):
            yr = re.search(r"\d{4}", r["year"]).group(0)
        sn = r["surname"].lower()
        cands = by_sur_yr.get((sn, yr), [])
        if cands:
            path = cands[0]
    if path:
        r["_path"] = path
        matched += 1
    else:
        r["_path"] = ""
        need_download += 1

print("matched to staged:", matched, "| need download:", need_download)

# 5) download missing from AARL Zotero web API (server-side file)
if need_download and ZKEY:
    got = 0
    for r in allrows:
        if r["_path"]:
            continue
        doi = r["doi"]
        att_key = r["att_key"]
        if not doi and not att_key:
            continue
        try:
            # try direct file download via att_key
            url = f"https://api.zotero.org/groups/735051/items/{att_key}/file"
            req = urllib.request.Request(url, headers={"Authorization": f"Bearer {ZKEY}"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
            if data[:4] == b"%PDF":
                path = os.path.join(STAGING, f"downloaded__{att_key}.pdf")
                with open(path, "wb") as f:
                    f.write(data)
                r["_path"] = path
                got += 1
                matched += 1
                need_download -= 1
                continue
        except urllib.error.HTTPError:
            pass
        except Exception:
            pass
        # fallback: search group by DOI, find any attachment child with file
        if doi:
            try:
                q = urllib.parse.quote(doi)
                surl = f"https://api.zotero.org/groups/735051/items?q={q}&format=json&limit=5&itemType=-attachment"
                sreq = urllib.request.Request(surl, headers={"Authorization": f"Bearer {ZKEY}"})
                with urllib.request.urlopen(sreq, timeout=30) as resp:
                    items = json.loads(resp.read().decode())
                for item in items:
                    if norm_doi(item["data"].get("DOI")) != doi:
                        continue
                    # get children
                    curl = f"https://api.zotero.org/groups/735051/items/{item['key']}/children?format=json"
                    creq = urllib.request.Request(curl, headers={"Authorization": f"Bearer {ZKEY}"})
                    with urllib.request.urlopen(creq, timeout=30) as resp:
                        kids = json.loads(resp.read().decode())
                    for kid in kids:
                        if kid["data"].get("contentType") != "application/pdf":
                            continue
                        furl = f"https://api.zotero.org/groups/735051/items/{kid['key']}/file"
                        freq = urllib.request.Request(furl, headers={"Authorization": f"Bearer {ZKEY}"})
                        with urllib.request.urlopen(freq, timeout=60) as resp:
                            data = resp.read()
                        if data[:4] == b"%PDF":
                            path = os.path.join(STAGING, f"downloaded__{kid['key']}.pdf")
                            with open(path, "wb") as f:
                                f.write(data)
                            r["_path"] = path
                            got += 1
                            matched += 1
                            need_download -= 1
                            break
                    if r["_path"]:
                        break
            except Exception:
                pass
        time.sleep(0.3)
    print("downloaded from AARL:", got, "| still missing:", need_download)

# 6) extract text
batch = []
for i, r in enumerate(allrows):
    text = extract_pdf(r["_path"]) if r["_path"] else ""
    batch.append({
        "att_key": r["att_key"], "doi": r["doi"], "title": r["title"],
        "surname": r["surname"], "year": r["year"], "folders": r["folders"],
        "folder_bucket": r["folder_bucket"], "path": r["_path"], "text": text,
    })
    if (i + 1) % 50 == 0:
        print("extracted", i + 1)

# 7) write batch JSONs (50 each)
CH = 50
for i in range(0, len(batch), CH):
    with open(os.path.join(PROG, f"all_extract_{i:03d}.json"), "w", encoding="utf-8") as f:
        json.dump(batch[i:i + CH], f, ensure_ascii=False)

n_good = sum(1 for b in batch if len(b["text"]) > 300 and not b["text"].startswith("EXTRACT-FAIL"))
n_none = sum(1 for b in batch if not b["path"])
print("total:", len(batch), "| good text:", n_good, "| no pdf:", n_none, "| thin:", len(batch) - n_good - n_none)
