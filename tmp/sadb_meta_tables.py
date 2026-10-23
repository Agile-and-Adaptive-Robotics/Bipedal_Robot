"""One-off discovery: list tables of the SADb base via the Airtable meta API.
Same PAT source as export_corpus.py. Read-only metadata listing."""
import json, re, urllib.request

txt = open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()
pat = re.search(r"PAT\s*:\s*(pat[^\s]+)", txt).group(1)
req = urllib.request.Request(
    "https://api.airtable.com/v0/meta/bases/appMQTnobUNRytIp7/tables",
    headers={"Authorization": "Bearer " + pat})
with urllib.request.urlopen(req, timeout=60) as r:
    d = json.loads(r.read().decode())
for t in d["tables"]:
    fields = [f["name"] for f in t.get("fields", [])]
    print(t["id"], "|", t["name"], "| records:", len(fields), "fields:", ", ".join(fields))
