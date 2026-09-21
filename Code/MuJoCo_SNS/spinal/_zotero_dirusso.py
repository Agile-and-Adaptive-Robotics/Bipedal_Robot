"""Di Russo lookup: check key permissions, then search whichever
libraries the key can read (AARL group 735051 / user 631450)."""
import io
import json
import sys
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

KEY = "34tQExoeiKS1yRe3Z1fbTF5z"


def get(url, key=KEY):
    req = urllib.request.Request(url, headers={
        "X-Zotero-API-Key": key, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode("utf-8"))


try:
    perms = get("https://api.zotero.org/keys/current")
    print("key:", perms.get("userID"), "access:",
          perms.get("access", {}), "groups:", perms.get("access", {})
          .get("groups", {}))
except Exception as e:
    print("key check failed:", e)

LIBS = [f"groups/735051", f"users/631450"]
for lib in LIBS:
    try:
        items = get(f"https://api.zotero.org/{lib}/items"
                    f"?q=Di%20Russo&limit=25&format=json")
    except Exception as e:
        print(f"{lib}: FAILED {e}")
        continue
    print(f"\n=== {lib}: {len(items)} hits ===")
    for it in items:
        d = it["data"]
        if d.get("itemType") in ("attachment", "note"):
            continue
        print(f"- {d.get('title', '?')[:88]}")
        print(f"  {d.get('itemType')} {d.get('date', '?')} "
              f"doi={d.get('DOI', '-')} key={it['key']}")
        try:
            kids = get(f"https://api.zotero.org/{lib}/items"
                       f"/{it['key']}/children?format=json&limit=10")
            for k in kids:
                kd = k["data"]
                if kd.get("itemType") == "attachment":
                    print(f"    attach: {kd.get('contentType')} "
                          f"{kd.get('filename', '')[:55]} key={k['key']}")
        except Exception as e:
            print("    children failed:", e)
