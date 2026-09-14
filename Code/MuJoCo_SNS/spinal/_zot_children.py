"""List child attachments of a Zotero item (usage: key)."""
import json
import sys
import urllib.request

BASE = "http://localhost:23119/api/users/0"
key = sys.argv[1]
with urllib.request.urlopen(f"{BASE}/items/{key}/children?format=json",
                            timeout=20) as r:
    kids = json.loads(r.read().decode("utf-8"))
for k in kids:
    d = k["data"]
    if d.get("itemType") == "attachment":
        print(k["key"], d.get("contentType", ""), d.get("title", "")[:70])
