"""Standalone echo verification for every curation_out batch (2026-09-28).

Per-record GETs (the list endpoint's records[] filter is silently ignored by
the REST API — it returns the unfiltered first 100; per-record GET cannot lie).
Checks Notes / Animal / Afferent Types / Animal Study Potential / Robot-Sim
Translation match the batch JSON exactly. Usage: myo python _verify_all.py [slugs...]
"""
import glob, json, os, re, sys, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "appMQTnobUNRytIp7"
PAPERS = "tblnnMrZszhboU4uD"


def pat():
    txt = open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()
    return re.search(r"PAT\s*:\s*(pat[^\s]+)", txt).group(1)


_PAT = pat()


def air(path, params=None):
    q = ""
    if params:
        q = "?" + urllib.parse.urlencode(params, doseq=True)
    req = urllib.request.Request("https://api.airtable.com/v0/" + path + q,
                                 headers={"Authorization": "Bearer " + _PAT})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


slugs = sys.argv[1:] or sorted(os.path.basename(p)[:2] for p in glob.glob(
    os.path.join(HERE, "curation_out", "batch_*.json")))
tot_ok = tot_bad = 0
for slug in slugs:
    batch = json.load(open(os.path.join(HERE, f"curation_out/batch_{slug}.json"),
                           encoding="utf-8-sig"))
    bad = 0
    for e in batch:
        d = air(f"{BASE}/{PAPERS}/{e['id']}")
        f = d.get("fields", {})
        an = set(x["name"] if isinstance(x, dict) else x for x in f.get("Animal", []))
        af = set(x["name"] if isinstance(x, dict) else x for x in f.get("Afferent Types", []))
        good = ((f.get("Notes") or "").strip()[:200] == (e.get("notes") or "").strip()[:200]
                and an == set(e.get("animals") or [])
                and af == set(e.get("afferents") or [])
                and (f.get("Animal Study Potential") or "").strip() == (e.get("animal_study") or "").strip()
                and (f.get("Robot/Sim Translation") or "").strip() == (e.get("robot_sim") or "").strip())
        if good:
            tot_ok += 1
        else:
            bad += 1
            tot_bad += 1
            print(f"  MISMATCH {slug} {e['id']}: notes_written={bool((f.get('Notes') or '').strip())} "
                  f"expected_notes={bool((e.get('notes') or '').strip())} "
                  f"an={sorted(an)}vs{sorted(set(e.get('animals') or []))} "
                  f"af={sorted(af)}vs{sorted(set(e.get('afferents') or []))}")
        time.sleep(0.22)
    print(f"batch {slug}: {len(batch) - bad}/{len(batch)} exact")
print(f"TOTAL echo: {tot_ok} exact, {tot_bad} mismatched")
