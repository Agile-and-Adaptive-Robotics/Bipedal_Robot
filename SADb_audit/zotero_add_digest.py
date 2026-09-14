#!/usr/bin/env python
# Build Zotero web-API payloads for the RR digest shortlist (SAD collection RBG8ZVYP).
# Usage (when Ben supplies his API key, userID 631450):
#   python zotero_add_digest.py --key XXXX            # dry run by default
#   python zotero_add_digest.py --key XXXX --post     # actually POST
# The key is passed on the command line only, never stored in a file (Ben's rule).
import json, argparse, subprocess, os, tempfile, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "batch3", "rr_digest_shortlist.json")

with open(SRC, encoding="utf-8") as f:
    rows = json.load(f)

# family, given from Crossref "authors" ("Family Given; ...")
def creators(s):
    out = []
    for part in (s or "").split(";"):
        part = part.strip()
        if not part:
            continue
        bits = part.rsplit(" ", 1)
        if len(bits) == 2:
            out.append({"creatorType": "author", "firstName": bits[1], "lastName": bits[0]})
        else:
            out.append({"creatorType": "author", "lastName": bits[0]})
    return out

items = []
for r in rows:
    if r.get("fetch") == "FAILED":
        continue
    items.append({
        "itemType": "journalArticle",
        "title": r.get("title", ""),
        "creators": creators(r.get("authors", ""))[:6],
        "publicationTitle": r.get("container", ""),
        "date": str(r.get("year") or ""),
        "DOI": r.get("doi", ""),
        "collections": ["RBG8ZVYP"],
        "tags": [],
        "relations": {},
    })

payload_path = os.path.join(HERE, "batch3", "zotero_digest_items.json")
with open(payload_path, "w", encoding="utf-8") as f:
    json.dump(items, f, indent=1, ensure_ascii=False)
print("wrote", payload_path, "with", len(items), "items")

ap = argparse.ArgumentParser()
ap.add_argument("--key")
ap.add_argument("--post", action="store_true")
a = ap.parse_args()
if a.key and a.post:
    # Zotero web API: POST /users/631450/items, JSON body, api-version header
    url = "https://api.zotero.org/users/631450/items"
    cmd = ["curl.exe", "-s", "-m", "60", "-X", "POST", url,
           "-H", "Authorization: Bearer " + a.key,
           "-H", "Zotero-API-Version: 3",
           "-H", "Content-Type: application/json",
           "--data-binary", "@" + payload_path]
    r = subprocess.run(cmd, capture_output=True)
    body = r.stdout.decode("utf-8", "replace")
    out_path = os.path.join(HERE, "batch3", "zotero_digest_post_result.json")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(body)
    print("POST result ->", out_path)
    print(body[:1500])
else:
    print("dry run: payload built; rerun with --key <APIKEY> --post to write")
