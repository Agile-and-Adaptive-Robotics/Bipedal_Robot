import os, urllib.parse, urllib.request

STAGING = r"D:\sadb_pdf_staging"
TUNNEL = "https://richard-ooo-deutschland-vertical.trycloudflare.com"

files = sorted(os.listdir(STAGING))
ascii_file = next(f for f in files if all(ord(c) < 128 for c in f) and f.endswith(".pdf"))
uni_file = next((f for f in files if not all(ord(c) < 128 for c in f) and f.endswith(".pdf")), None)

for f in [ascii_file, uni_file]:
    if not f:
        continue
    url = TUNNEL + "/" + urllib.parse.quote(f)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "test"})
        with urllib.request.urlopen(req, timeout=30) as r:
            head = r.read(8)
            print(f, "-> HTTP", r.status, "head:", head[:5])
    except Exception as e:
        print(f, "-> FAIL:", str(e)[:120])
