"""Fetch Shevtsova + Shinohara items and PDFs from the AARL Zotero
group (735051). Query-string auth (header form gets stripped)."""
import io
import json
import sys
import urllib.parse
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
KEY = "34tQExoeiKS1yRe3Z1fbTF5z"


def get(url):
    sep = "&" if "?" in url else "?"
    req = urllib.request.Request(f"{url}{sep}key={KEY}&v=3",
                                 headers={"Zotero-API-Version": "3"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.loads(r.read().decode("utf-8"))


def find(q):
    uq = urllib.parse.quote(q)
    return [i for i in get(
        f"https://api.zotero.org/groups/735051/items?q={uq}&limit=15"
        f"&format=json")
        if i["data"].get("itemType") not in ("attachment", "note")]


for label, q in (("Shevtsova laminar", "Shevtsova laminar"),
                 ("Shevtsova RP107480", "RP107480"),
                 ("Shinohara", "Shinohara")):
    print(f"=== {label} ===")
    for it in find(q):
        d = it["data"]
        print(f"- {d.get('title', '?')[:95]}")
        print(f"  {d.get('itemType')} {d.get('date', '?')} "
              f"doi={d.get('DOI', '-')} key={it['key']}")
        for k in get(f"https://api.zotero.org/groups/735051/items/"
                     f"{it['key']}/children?format=json&limit=10"):
            if k["data"].get("itemType") == "attachment":
                print(f"    attach: {k['data'].get('contentType')} "
                      f"{k['data'].get('filename', '')[:55]} "
                      f"key={k['key']}")
