"""Find the Shevtsova 2026 paper in Zotero local API; print item metadata
and attachment keys (overnight session, Phase 1)."""
import json
import sys
import urllib.request

BASE = "http://localhost:23119/api/users/0"


def get(url):
    with urllib.request.urlopen(url, timeout=20) as r:
        return json.loads(r.read().decode("utf-8"))


def main():
    q = sys.argv[1] if len(sys.argv) > 1 else "Shevtsova"
    items = get(f"{BASE}/items?q={q}&format=json&limit=25")
    print(f"{len(items)} items for query {q!r}")
    for it in items:
        d = it["data"]
        key = it["key"]
        typ = d.get("itemType", "?")
        title = d.get("title", "")[:110]
        if typ in ("attachment", "note") and it.get("links", {}).get("up"):
            print(f"  [child] {typ} key={key} of={it['links']['up'][0]['href'].split('/')[-1]}")
            continue
        print(f"  {typ:12s} key={key}  {title}")
        if typ == "journalArticle":
            print(f"      DOI={d.get('DOI','')} date={d.get('date','')}")
        # children
        try:
            kids = get(f"{BASE}/items/{key}/children?format=json")
            for k in kids:
                kd = k["data"]
                if kd.get("itemType") == "attachment":
                    ct = kd.get("contentType", "")
                    print(f"      ATT key={k['key']} type={ct} "
                          f"title={kd.get('title','')[:60]}")
        except Exception as e:
            print(f"      (children: {e})")


if __name__ == "__main__":
    main()
