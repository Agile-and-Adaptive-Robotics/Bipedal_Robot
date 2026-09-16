import os, json, urllib.request, urllib.error

p = os.environ.get("AT_PAT", "")
# write probe into the dead 'Field 10' column (harmless)
req = urllib.request.Request(
    "https://api.airtable.com/v0/appMQTnobUNRytIp7/Papers/recz4jiOLNrm1PACf",
    data=json.dumps({"fields": {"fldnoM8tlsUo7SPJA": "pat-write-probe"}}).encode(),
    headers={"Authorization": "Bearer " + p, "Content-Type": "application/json"},
    method="PATCH")
try:
    print("PATCH: HTTP", urllib.request.urlopen(req, timeout=30).status)
except urllib.error.HTTPError as e:
    print("PATCH: HTTP", e.code, e.read().decode()[:300])
