# Probe the OpenAlex anonymous-tier budget: one cheap /works query per attempt.
import sys, time, json, urllib.request, urllib.error

url = ("https://api.openalex.org/works?filter=title_and_abstract.search:spinal%20locomotion,"
       "from_publication_date:2026-09-01,to_publication_date:2026-10-22&per-page=1")
tries = int(sys.argv[1]) if len(sys.argv) > 1 else 10
gap = int(sys.argv[2]) if len(sys.argv) > 2 else 120
for i in range(1, tries + 1):
    req = urllib.request.Request(url, headers={"User-Agent": "SADbMonthlyWatch/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.loads(r.read().decode("utf-8", "replace"))
        print("attempt %d/%d OK count=%s" % (i, tries, d["meta"]["count"]))
        sys.exit(0)
    except urllib.error.HTTPError as e:
        print("attempt %d/%d HTTP %d" % (i, tries, e.code))
    except Exception as e:
        print("attempt %d/%d %s %s" % (i, tries, type(e).__name__, e))
    if i < tries:
        time.sleep(gap)
sys.exit(1)
