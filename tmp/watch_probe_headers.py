import urllib.request, urllib.error, sys

url = ("https://api.openalex.org/works?filter=title_and_abstract.search:spinal%20locomotion,"
       "from_publication_date:2026-09-01,to_publication_date:2026-10-22&per-page=1")
req = urllib.request.Request(url, headers={"User-Agent": "SADbMonthlyWatch/1.0"})
try:
    with urllib.request.urlopen(req, timeout=30) as r:
        print("OK", r.status)
        print(dict(r.headers))
except urllib.error.HTTPError as e:
    print("HTTP", e.code)
    print("headers:", dict(e.headers))
    body = e.read().decode("utf-8", "replace")[:500]
    print("body:", body)
