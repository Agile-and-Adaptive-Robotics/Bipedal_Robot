"""Execute Ben's 2026-10-02 flags-comment deletions + find 4 review orphans.

1. Snapshot + delete the flagged papers: Asif 2012, Dean 2013, Gao 2019,
   Gerstner 2009, Gomar 2014, Gubina 1974 (+ Blackburn 2004 if still live —
   Ben believes it's already gone).
2. Read the Review Papers table, list records whose Paper link is EMPTY
   (orphaned by the 290 cut) — candidates for Ben's '4 shitty reviews'.
"""
import json, os, re, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE, PAPERS = "appMQTnobUNRytIp7", "tblnnMrZszhboU4uD"
REVIEWS = "tblSEubKcRId4wYMK"
PAT = re.search(r"PAT\s*:\s*(pat[^\s]+)",
                open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()).group(1)


def req(path, params=None, method="GET"):
    q = "?" + urllib.parse.urlencode(params, doseq=True) if params else ""
    r = urllib.request.Request("https://api.airtable.com/v0/" + path + q,
                               headers={"Authorization": "Bearer " + PAT}, method=method)
    with urllib.request.urlopen(r, timeout=90) as resp:
        return json.loads(resp.read().decode())


TARGETS = {"reciUGAHPamTLce0k": "Asif 2012 (Ben: delete)",
           "recUWzWrumdQrLybs": "Dean 2013 (Ben: delete)",
           "recvwIzskqiUp9ikx": "Gao 2019 (Ben: delete)",
           "recKCncNuVAa7rZ4g": "Gerstner 2009 (Ben: delete)",
           "recwAaYQRw2Psco5n": "Gomar 2014 (Ben: delete)",
           "recZVayeFS2NnRjwz": "Gubina 1974 (Ben: delete)",
           "recxUQ8JuDHyBXIk0": "Blackburn 2004 (Ben: thought already deleted)"}

# live check + snapshot
records = {r["id"]: r for r in json.load(open(os.path.join(HERE, "export", "sadb_export.json"),
                                              encoding="utf-8"))}
live = []
for rid, why in TARGETS.items():
    try:
        got = req(f"{BASE}/{PAPERS}/{rid}")
        if got.get("fields"):
            live.append(rid)
            print(f"live, will delete: {why}")
            continue
    except urllib.error.HTTPError as e:
        if e.code in (403, 404):                 # Airtable masks not-found as 403
            print(f"already gone: {why}")
            continue
        raise
    print(f"already gone: {why}")
snap = [records[rid] for rid in live if rid in records]
json.dump(snap, open(os.path.join(HERE, "archive", "deleted_records_20261002b.json"), "w"),
          ensure_ascii=False, indent=1)

if live:
    q = "?" + urllib.parse.urlencode([("records[]", x) for x in live])
    req(f"{BASE}/{PAPERS}{q}", method="DELETE")
    print(f"deleted {len(live)} flagged papers")

# --- Review Papers orphans ---
orphans, offset = [], None
while True:
    p = [("pageSize", 100), ("fields[]", ["Name", "Paper"])] + ([("offset", offset)] if offset else [])
    d = req(f"{BASE}/{REVIEWS}", p)
    for r in d.get("records", []):
        if not r["fields"].get("Paper"):
            orphans.append((r["id"], r["fields"].get("Name", "?")))
    offset = d.get("offset")
    if not offset:
        break
print(f"\nReview Papers orphans (no Paper link after the cuts): {len(orphans)}")
for rid, name in orphans[:20]:
    print("  ", rid, name)
json.dump({rid: name for rid, name in orphans},
          open(os.path.join(HERE, "review_orphans_20261002.json"), "w"), indent=1)
