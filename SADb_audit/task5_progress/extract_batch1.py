"""Task 5 batch extraction: for a list of att_keys, extract text from staged PDFs
and pull creators/year/DOI from the Zotero parent item via local API.
Writes per-paper extract files + one batch JSON; prints compact digests.
"""
import csv, json, os, re, subprocess, sys
sys.path.insert(0, r"D:\pdfx_tmp")
import pypdf

SAD = r"D:\Github\Bipedal_Robot\SADb_audit"
STAGING = r"D:\sadb_pdf_staging"
PROG = os.path.join(SAD, "task5_progress")
os.makedirs(PROG, exist_ok=True)

def zot_local(url):
    p = subprocess.run(["curl.exe", "-s", "-m", "30", url], capture_output=True)
    try:
        return json.loads(p.stdout.decode("utf-8", "replace"))
    except Exception:
        return None

def extract_pdf(path, max_pages=3, max_chars=3500):
    try:
        r = pypdf.PdfReader(path)
        chunks = []
        n = 0
        for i, pg in enumerate(r.pages[:max_pages]):
            t = pg.extract_text() or ""
            chunks.append(t)
            n += len(t)
            if n > max_chars:
                break
        txt = "\n".join(chunks)
        txt = re.sub(r"\s+", " ", txt).strip()
        return txt[:max_chars]
    except Exception as e:
        return "EXTRACT-FAIL: " + str(e)[:100]

def main(keys):
    rows = {}
    with open(os.path.join(SAD, "author_fix", "aarl_gap_by_folder.csv"), encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if r["att_key"] in keys:
                rows[r["att_key"]] = r
    batch = []
    for k in keys:
        r = rows.get(k)
        if not r:
            print("MISSING ROW", k)
            continue
        # find staged file
        cand = None
        for f in os.listdir(STAGING):
            if f.endswith("__" + k + ".pdf"):
                cand = os.path.join(STAGING, f)
                break
        # also check pdfs subfolder by att_key
        if not cand:
            sub = os.path.join(STAGING, "pdfs")
            if os.path.isdir(sub):
                for f in os.listdir(sub):
                    if f.endswith("__" + k + ".pdf"):
                        cand = os.path.join(sub, f)
                        break
        text = ""
        if cand:
            text = extract_pdf(cand)
        else:
            text = "NO-STAGED-PDF"
        # parent metadata from Zotero local API
        scope = r["library"]
        att = zot_local(f"http://localhost:23119/api/{scope}/items/{k}?format=json")
        parent = att["data"]["parentItem"] if att and att.get("data") else ""
        creators, doi, year, title = [], r["doi"], r["year"], r["title"]
        if parent:
            p = zot_local(f"http://localhost:23119/api/{scope}/items/{parent}?format=json")
            if p and p.get("data"):
                d = p["data"]
                doi = (d.get("DOI") or "").strip() or doi
                year = str(d.get("date") or "")[:4] or year
                title = " ".join((d.get("title") or "").split()) or title
                for c in d.get("creators") or []:
                    if c.get("lastName"):
                        creators.append(c["lastName"])
                    elif c.get("name"):
                        creators.append(c["name"].strip().split(" ")[-1])
        rec = {
            "att_key": k, "doi": doi, "title": title, "surname_csv": r["surname"],
            "year": year, "creators": creators, "folders": r["folders"],
            "folder_bucket": r["folder_bucket"], "zot_scope": scope,
            "zot_parent": parent, "staged_path": cand or "", "text": text,
        }
        batch.append(rec)
        print("=" * 100)
        print(k, "|", r["surname"], year, "|", (doi or "no-DOI")[:40])
        print("TITLE:", title[:110])
        print("AUTHORS:", creators[:12])
        print("TEXT:", text[:1500] if text else "(empty)")
    with open(os.path.join(PROG, "batch_sf_pf_extract.json"), "w", encoding="utf-8") as f:
        json.dump(batch, f, ensure_ascii=False, indent=1)
    print("WROTE batch_sf_pf_extract.json", len(batch))

if __name__ == "__main__":
    keys = ["D9NYJVGY","CPWAEC5K","4EMAIY84","FG6SXTAX","C5ZKWUUK","9BRRW9JI","RI68M9LY","YDL8YIDS",
            "TSTXZGW7","VZ9F38IZ","R4J2XKXT","7DF6JYJW","XI89DIYM","QD6TP8D5","IRX5WUBK","UDU4DVIN",
            "RIN8JZM4","NCSIFLYR","BNJN4GIA","22RM4EQK","TAZUXF69","ZMS66V28"]
    main(keys)
