"""Search Ben's USER Zotero library (631450) for Shevtsova."""
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


for q in ("Shevtsova", "laminar organization locomotor", "RP107480"):
    uq = urllib.parse.quote(q)
    try:
        items = get("https://api.zotero.org/users/631450/items"
                    f"?q={uq}&limit=10&format=json")
    except Exception as e:
        print(q, "FAILED", e)
        continue
    print(f"=== {q}: {len(items)} hits ===")
    for it in items:
        d = it["data"]
        if d.get("itemType") in ("attachment", "note"):
            continue
        print(f"- {d.get('title', '?')[:95]}")
        print(f"  {d.get('itemType')} {d.get('date', '?')} "
              f"doi={d.get('DOI', '-')} key={it['key']}")
        for k in get(f"https://api.zotero.org/users/631450/items/"
                     f"{it['key']}/children?format=json&limit=10"):
            if k["data"].get("itemType") == "attachment":
                print(f"    attach: {k['data'].get('contentType')} "
                      f"{str(k['data'].get('filename', ''))[:50]} "
                      f"key={k['key']}")
