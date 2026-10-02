"""Add 'Invertebrates' via typecast write (metadata API rejects option edits).

Sets Animals = [Vertebrates, Invertebrates] on Ting & Chiel 2017 — the paper
Ben's ruling was about. typecast=true makes Airtable create the new option.
"""
import json, re, urllib.request

PAT = re.search(r"PAT\s*:\s*(pat[^\s]+)",
                open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()).group(1)
URL = "https://api.airtable.com/v0/appMQTnobUNRytIp7/tblnnMrZszhboU4uD"
body = {"records": [{"id": "recMLbRa9efnWHXEx",
                     "fields": {"Animal": ["Vertebrates", "Invertebrates"]}}],
        "typecast": True}
r = urllib.request.Request(URL, data=json.dumps(body).encode(),
                           headers={"Authorization": "Bearer " + PAT,
                                    "Content-Type": "application/json"}, method="PATCH")
with urllib.request.urlopen(r, timeout=60) as resp:
    d = json.loads(resp.read().decode())
print("animals now:", d["records"][0]["fields"]["Animal"])
