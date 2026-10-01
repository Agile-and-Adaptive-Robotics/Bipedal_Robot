"""Delete the four Ben-approved junk fields via the METADATA API (2026-10-01).

The plain data API 404s on field DELETE; field deletion lives under
/v0/meta/bases/{baseId}/tables/{tableId}/fields/{fieldId} and needs a
schema-scoped PAT. Falls back to a clear message for Ben's 2-click UI deletion.
"""
import re, urllib.request

PAT = re.search(r"PAT\s*:\s*(pat[^\s]+)",
                open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()).group(1)
targets = [
    ("tblnnMrZszhboU4uD", "fldnoM8tlsUo7SPJA", "Papers: Field 10 (unused)"),
    ("tblnnMrZszhboU4uD", "fldXcGzCtnfISOcEn", "Papers: Field 11 (unused)"),
    ("tblnnMrZszhboU4uD", "fldh983rtt2YtMZQX", "Papers: Models copy (unused, empty)"),
    ("tblot5mo4s5KgN5le", "fldjcTBtIC1xHzl7j", "Feedback: Models copy (unused on all 38)"),
]
for tid, fid, name in targets:
    req = urllib.request.Request(
        f"https://api.airtable.com/v0/meta/bases/appMQTnobUNRytIp7/tables/{tid}/fields/{fid}",
        headers={"Authorization": "Bearer " + PAT}, method="DELETE")
    try:
        urllib.request.urlopen(req, timeout=60)
        print(name, "-> DELETED")
    except Exception as e:
        print(name, "-> FAILED", str(e)[:110])
