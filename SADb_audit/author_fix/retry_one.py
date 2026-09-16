import os, json, urllib.parse, urllib.request

AT_PAT = os.environ["AT_PAT"]
TUNNEL = "https://richard-ooo-deutschland-vertical.trycloudflare.com"
STAGING = r"D:\sadb_pdf_staging"

rid = "recz4jiOLNrm1PACf"  # 10.1002/cne.20711 (Álvarez)
match = [f for f in os.listdir(STAGING) if f.startswith("Álvarez")]
print("candidates:", match)
f = match[0]
url = TUNNEL + "/" + urllib.parse.quote(f)
req = urllib.request.Request(
    f"https://api.airtable.com/v0/appMQTnobUNRytIp7/Papers/{rid}",
    data=json.dumps({"fields": {"Attachments": [{"url": url}]}}).encode(),
    headers={"Authorization": "Bearer " + AT_PAT, "Content-Type": "application/json"},
    method="PATCH")
with urllib.request.urlopen(req, timeout=90) as r:
    resp = json.loads(r.read().decode())
atts = (resp.get("fields") or {}).get("Attachments") or []
print("attached:", [(a.get("filename"), a.get("size")) for a in atts])
