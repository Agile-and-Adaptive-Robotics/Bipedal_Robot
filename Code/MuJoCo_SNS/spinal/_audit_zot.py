"""Pull abstracts for specific Zotero keys + run topic queries (audit)."""
import json
import sys
import urllib.parse
import urllib.request

BASE = "http://localhost:23119/api/users/0"


def get(url):
    with urllib.request.urlopen(url, timeout=25) as r:
        return json.loads(r.read().decode("utf-8"))


def abstract(key):
    try:
        d = get(f"{BASE}/items/{key}?format=json")["data"]
        return d.get("abstractNote", "").replace("\n", " ")[:700], \
            [c.get("lastName") or "?" for c in d.get("creators", [])][:3], \
            d.get("date", "")
    except Exception as e:
        return f"(err {e})", [], ""


if sys.argv[1] == "abs":
    for key in sys.argv[2:]:
        a, cr, dt = abstract(key)
        print(f"== {key} [{', '.join(cr)} {dt}] ==")
        print(a)
        print()
elif sys.argv[1] == "q":
    for q in sys.argv[2:]:
        url = f"{BASE}/items?q={urllib.parse.quote(q)}&format=json&limit=8"
        try:
            items = get(url)
            print(f"== query {q!r}: {len(items)} ==")
            for it in items:
                d = it["data"]
                if d.get("itemType") in ("attachment", "note"):
                    continue
                yr = d.get("date", "")[:4]
                print(f"  {it['key']} {yr} {d.get('title','')[:95]}")
        except Exception as e:
            print(f"== query {q!r} ERR {e}")
