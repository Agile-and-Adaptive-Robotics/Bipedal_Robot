"""Fetch Deng 2019 (Biomimetics 4(1):21) PDF from the AARL Zotero group.

Key is read from D:\Github\api_credentials_local.txt (outside the repo).
Outputs the PDF to tmp\bg15_figs\deng_2019.pdf.
"""
import io
import json
import re
import sys
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

CRED = r"D:\Github\api_credentials_local.txt"
txt = open(CRED, encoding="utf-8", errors="replace").read()
m = re.search(r"API\s*key\s*:\s*([A-Za-z0-9]+)", txt)
KEY = m.group(1)
GROUP = "735051"
OUT = r"D:\Github\Bipedal_Robot\tmp\bg15_figs\deng_2019.pdf"


def get_json(url):
    req = urllib.request.Request(url, headers={
        "Zotero-API-Version": "3",
        "Zotero-API-Key": KEY,
        "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode("utf-8"))


items = get_json(f"https://api.zotero.org/groups/{GROUP}/items"
                 f"?q=Deng%20neuromechanical&limit=25&format=json")
cands = []
for it in items:
    d = it["data"]
    if d.get("itemType") in ("attachment", "note"):
        continue
    title = d.get("title", "")
    doi = (d.get("DOI") or "").lower()
    if "neuromechanical model of rat hindlimb" in title.lower() \
            and "biomimetics4010021" in doi:
        cands.insert(0, it)
        print("DOI HIT:", title, "| key:", it["key"], "| doi:", doi)
    elif "neuromechanical model of rat hindlimb" in title.lower():
        cands.append(it)
        print("hit:", title, "| key:", it["key"], "| doi:", doi)

if not cands:
    print("no top-level hit; dumping titles")
    for it in items:
        d = it["data"]
        print("-", d.get("title", "?")[:90], d.get("itemType"))
    sys.exit(1)

item = cands[0]
kids = get_json(f"https://api.zotero.org/groups/{GROUP}/items"
                f"/{item['key']}/children?format=json&limit=25")
pdf_key = None
for k in kids:
    kd = k["data"]
    if kd.get("contentType") == "application/pdf":
        pdf_key = k["key"]
        print("pdf attachment:", kd.get("filename"), "key:", pdf_key)
        break

if not pdf_key:
    print("NO PDF ATTACHMENT on the Zotero item")
    sys.exit(1)

req = urllib.request.Request(
    f"https://api.zotero.org/groups/{GROUP}/items/{pdf_key}/file",
    headers={"Zotero-API-Version": "3", "Zotero-API-Key": KEY})
with urllib.request.urlopen(req, timeout=60) as r:
    data = r.read()
open(OUT, "wb").write(data)
print("wrote", OUT, len(data), "bytes")
