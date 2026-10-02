"""Delete the 47 orphaned Review Papers records + audit Models orphans (2026-10-02)."""
import json, os, re, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE, REVIEWS, MODELS = "appMQTnobUNRytIp7", "tblSEubKcRId4wYMK", "tblsBq9IEv7dZe6fn"
PAT = re.search(r"PAT\s*:\s*(pat[^\s]+)",
                open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()).group(1)


def req(path, params=None, method="GET"):
    q = "?" + urllib.parse.urlencode(params, doseq=True) if params else ""
    r = urllib.request.Request("https://api.airtable.com/v0/" + path + q,
                               headers={"Authorization": "Bearer " + PAT}, method=method)
    with urllib.request.urlopen(r, timeout=90) as resp:
        return json.loads(resp.read().decode())


orphans = json.load(open(os.path.join(HERE, "review_orphans_20261002.json"), encoding="utf-8"))
ids = list(orphans.keys())
print(f"deleting {len(ids)} orphaned Review Papers records")
ok = 0
for i in range(0, len(ids), 10):
    q = "?" + urllib.parse.urlencode([("records[]", x) for x in ids[i:i + 10]])
    req(f"{BASE}/{REVIEWS}{q}", method="DELETE")
    ok += len(ids[i:i + 10])
print(f"deleted {ok}/{len(ids)}")

# Models orphan audit (report only — Ben holds Models stubs)
m_orphans, offset = [], None
while True:
    p = [("pageSize", 100), ("fields[]", ["Name", "Paper"])] + ([("offset", offset)] if offset else [])
    d = req(f"{BASE}/{MODELS}", p)
    for r in d.get("records", []):
        f = r["fields"]
        has_data = bool((f.get("Notes") or "").strip() or f.get("Rules Used")
                        or f.get("Papers Cited") or f.get("Attachments"))
        if not f.get("Paper") and not has_data:
            m_orphans.append((r["id"], f.get("Name", "?")))
    offset = d.get("offset")
    if not offset:
        break
print(f"Models orphans (no Paper, no data): {len(m_orphans)} (report only)")
json.dump({rid: n for rid, n in m_orphans},
          open(os.path.join(HERE, "models_orphans_20261002.json"), "w"), indent=1)
