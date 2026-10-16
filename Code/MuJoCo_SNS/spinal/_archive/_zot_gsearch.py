"""Search AARL group library (735051) for items; list PDF attachments.
Usage: _zot_gsearch.py <query>       — top item hits
       _zot_gsearch.py <itemkey> pdf — list PDF attachments of an item
       _zot_gsearch.py <itemkey> get <out.pdf> — download first PDF
"""
import io
import json
import sys
import urllib.parse
import urllib.request

KEY = "34tQExoeiKS1yRe3Z1fbTF5z"
GID = "735051"
HDR = {"Zotero-API-Key": KEY, "Zotero-API-Version": "3"}
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")


def get(url):
    req = urllib.request.Request(url, headers=HDR)
    with urllib.request.urlopen(req, timeout=60) as r:
        raw = r.read()
    ct = r.headers.get("Content-Type", "")
    if "json" in ct:
        return json.loads(raw.decode("utf-8"))
    return raw


if len(sys.argv) >= 3 and sys.argv[2] == "pdf":
    kids = get(f"https://api.zotero.org/groups/{GID}/items/{sys.argv[1]}/"
               f"children?format=json")
    for c in kids:
        d = c["data"]
        print(f"  {c['key']} [{d['itemType']}] {d.get('contentType','')} "
              f"{d.get('title','')[:60]}")
elif len(sys.argv) >= 4 and sys.argv[2] == "get":
    out = sys.argv[3]
    raw = get(f"https://api.zotero.org/groups/{GID}/items/{sys.argv[1]}/file")
    with open(out, "wb") as f:
        f.write(raw)
    print(f"downloaded {out} ({len(raw)} bytes)")
else:
    q = urllib.parse.quote(" ".join(sys.argv[1:]))
    url = (f"https://api.zotero.org/groups/{GID}/items?format=json"
           f"&q={q}&limit=15&itemType=-attachment")
    items = get(url)
    print(f"{len(items)} hits for {sys.argv[1:]!r}")
    for it in items:
        d = it["data"]
        print(f"  {it['key']} [{d['itemType']}] "
              f"{it['meta'].get('creatorSummary','')} {d.get('date','')[:4]}"
              f"  {d.get('title','')[:85]}")
