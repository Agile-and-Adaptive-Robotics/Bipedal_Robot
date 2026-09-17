import os, json, urllib.request, urllib.error

p = os.environ.get("AT_PAT", "")
base = "appMQTnobUNRytIp7"

# read test
try:
    req = urllib.request.Request(f"https://api.airtable.com/v0/{base}/Papers?maxRecords=1",
                                 headers={"Authorization": "Bearer " + p})
    print("READ: HTTP", urllib.request.urlopen(req, timeout=30).status)
except urllib.error.HTTPError as e:
    print("READ: HTTP", e.code, e.read().decode()[:300])

# upload test
try:
    req = urllib.request.Request(f"https://api.airtable.com/v0/{base}/attachments",
                                 data=json.dumps({"contentType": "application/pdf", "filename": "probe.pdf"}).encode(),
                                 headers={"Authorization": "Bearer " + p, "Content-Type": "application/json"},
                                 method="POST")
    print("UPLOAD: HTTP", urllib.request.urlopen(req, timeout=30).status)
except urllib.error.HTTPError as e:
    print("UPLOAD: HTTP", e.code, e.read().decode()[:400])
