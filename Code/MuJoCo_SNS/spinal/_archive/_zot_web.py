"""Zotero WEB API access test (AARL group library, read-only).
Usage: _zot_web.py profile|groups|gcoll <groupid>|gitems <groupid> <collkey>|gfile <groupid> <itemkey> <out>
Credentials in D:\\Github\\api_credentials_local.txt.
"""
import io
import json
import sys
import urllib.request

KEY = "34tQExoeiKS1yRe3Z1fbTF5z"
USER = "631450"
HDR = {"Zotero-API-Key": KEY, "Zotero-API-Version": "3"}
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")


def get_json(url):
    req = urllib.request.Request(url, headers=HDR)
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode("utf-8"))


mode = sys.argv[1]
if mode == "profile":
    d = get_json(f"https://api.zotero.org/users/{USER}/settings?limit=1")
    print("user endpoint OK:", d)
    gs = get_json(f"https://api.zotero.org/users/{USER}/groups?limit=50")
    for g in gs:
        print(f"group {g['data']['id']}: {g['data']['name']} "
              f"({g['data'].get('type','')})")
elif mode == "gcoll":
    gid = sys.argv[2]
    cols = get_json(f"https://api.zotero.org/groups/{gid}/collections?limit=100")
    for c in cols:
        print(f"{c['key']}  {c['data']['name']}  "
              f"(top-level={c['data']['parentCollection'] is False})")
elif mode == "gitems":
    gid, coll = sys.argv[2], sys.argv[3]
    items = get_json(f"https://api.zotero.org/groups/{gid}/items?format=json"
                     f"&collection={coll}&limit=25")
    print(f"{len(items)} items in {coll}")
    for it in items:
        d = it["data"]
        print(f"  {it['key']} {it['meta'].get('creatorSummary','')} "
              f"{d.get('date','')[:4]} [{d['itemType']}] "
              f"{d.get('title','')[:70]}")
elif mode == "gfile":
    gid, item, out = sys.argv[2], sys.argv[3], sys.argv[4]
    url = f"https://api.zotero.org/groups/{gid}/items/{item}/file"
    req = urllib.request.Request(url, headers=HDR)
    with urllib.request.urlopen(req, timeout=120) as r, open(out, "wb") as f:
        f.write(r.read())
    import os
    print(f"downloaded {out} ({os.path.getsize(out)} bytes)")
