"""Supervised writer for curation batches (2026-09-28).

Usage: myo python write_curation.py curation_out/batch_XX.json [--dry]

Reads subagent-produced curation JSON (list of {id, notes, animals, afferents,
feedback, animal_study, robot_sim, is_model, model_name, is_review, review_name,
flags}), VALIDATES everything against the fixed vocabularies + live Feedback
names, writes Notes/Animal/Afferent Types/Animal Study Potential/Robot-Sim
Translation/Feedback to Airtable in batches of 10, re-GETs and echo-verifies,
creates Models / Review Papers twins where clearly warranted, and appends the
batch to curation_log.csv. Unknown feedback names are never written — they are
logged as proposals for Ben.

Validation gate (the fact-check):
  - notes 20..1500 chars; must not be an abstract paste (8-gram overlap < 40%)
  - animals/afferents restricted to the fixed vocabularies
  - animal_study / robot_sim <= 1000 chars
  - exact-duplicate notes across the batch are rejected
"""
import json, os, re, sys, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "appMQTnobUNRytIp7"
PAPERS = "tblnnMrZszhboU4uD"
F_NOTES = "fld3gPiUIKn26N6ji"
F_ANIMAL = "fld1x2BXLKIdA2dCw"
F_AFF = "fldcWSx0wyFJHjYZc"
F_ANST = "fldYtDbTshGMbRK6w"
F_ROBO = "fld5PGyWPjLYONLTD"
F_FB = "fldK2H6RfaSdLLLM7"

ANIMALS = {"Cat", "Dog", "Frog", "Hexapod", "Human", "Insects", "Lamprey", "Mammals",
           "Mice", "Rat", "Salamander", "Stick Insect", "Vertebrates", "Zebrafish",
           "Arthropods", "Cockroach", "Turtle"}
AFFERENTS = {"Type 1 (legacy)", "Ia", "Ib", "II", "III/IV", "Mechanosensory",
             "Cutaneous", "Heat", "Nociceptive", "Flexor reflex afferents"}


def pat():
    txt = open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()
    m = re.search(r"PAT\s*:\s*(pat[^\s]+)", txt)
    if not m:
        raise SystemExit("PAT not found")
    return m.group(1)


_PAT = pat()


def air(path, params=None):
    q = ""
    if params:
        q = "?" + urllib.parse.urlencode(params, doseq=True)
    req = urllib.request.Request("https://api.airtable.com/v0/" + path + q,
                                 headers={"Authorization": "Bearer " + _PAT})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def air_req(path, body, method):
    req = urllib.request.Request("https://api.airtable.com/v0/" + path,
                                 data=json.dumps(body).encode(),
                                 headers={"Authorization": "Bearer " + _PAT,
                                          "Content-Type": "application/json"},
                                 method=method)
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def ngrams(s, n=8):
    w = re.findall(r"[a-z]+", (s or "").lower())
    return {" ".join(w[i:i + n]) for i in range(max(0, len(w) - n + 1))}


def main(path, dry=False, verify_only=False):
    batch = json.load(open(os.path.join(HERE, path), encoding="utf-8-sig"))
    ground = json.load(open(os.path.join(HERE, "export", "grounding.json"), encoding="utf-8"))
    fb_names = {}
    for rec in air(f"{BASE}/tblot5mo4s5KgN5le", [("pageSize", 100)]).get("records", []):
        fb_names[rec["fields"].get("Name", "")] = rec["id"]

    errs, ok = [], []
    seen_notes = {}
    for i, e in enumerate(batch):
        rid = e.get("id", "")
        notes = (e.get("notes") or "").strip()
        if not (20 <= len(notes) <= 1500):
            errs.append(f"[{i}] {rid}: notes length {len(notes)} out of 20..1500")
            continue
        if notes in seen_notes:
            errs.append(f"[{i}] {rid}: duplicate notes of {seen_notes[notes]}")
            continue
        seen_notes[notes] = rid
        ab = (ground.get(rid, {}).get("abstract") or "")
        ga, gb = ngrams(notes), ngrams(ab)
        if ga and gb and len(ga & gb) / max(1, len(ga)) > 0.40:
            errs.append(f"[{i}] {rid}: notes look like an abstract paste "
                        f"({100*len(ga & gb)//max(1,len(ga))}% 8-gram overlap)")
            continue
        bad_an = set(e.get("animals") or []) - ANIMALS
        bad_af = set(e.get("afferents") or []) - AFFERENTS
        if bad_an or bad_af:
            errs.append(f"[{i}] {rid}: bad vocab animals={bad_an} afferents={bad_af}")
            continue
        if len(e.get("animal_study") or "") > 1000 or len(e.get("robot_sim") or "") > 1000:
            errs.append(f"[{i}] {rid}: animal_study/robot_sim > 1000 chars")
            continue
        fb_ids, proposals = [], []
        for nm in e.get("feedback") or []:
            nm = nm.strip()
            if nm in fb_names:
                fb_ids.append(fb_names[nm])
            else:
                proposals.append(nm)
        ok.append((e, fb_ids, proposals))

    print(f"validated: {len(ok)} ok, {len(errs)} rejected")
    for e in errs:
        print("  REJECT", e)
    if dry or not ok:
        return

    written, failed = 0, 0
    for i in range(0, len(ok), 10):
        chunk = ok[i:i + 10]
        body = {"records": [
            {"id": e["id"], "fields": {
                F_NOTES: e["notes"].strip(),
                F_ANIMAL: e.get("animals") or [],
                F_AFF: e.get("afferents") or [],
                F_ANST: (e.get("animal_study") or "").strip(),
                F_ROBO: (e.get("robot_sim") or "").strip(),
                F_FB: fbids,
            }} for e, fbids, _ in chunk], "typecast": False}
        try:
            air_req(f"{BASE}/{PAPERS}", body, "PATCH")
            written += len(chunk)
        except Exception as ex:
            failed += len(chunk)
            print("  WRITE FAIL", ex)
        time.sleep(0.4)

    # echo verification — PER-RECORD GETs: the list endpoint's records[] filter
    # is silently ignored by the REST API (returns the unfiltered first 100),
    # so batch verification must fetch each record individually.
    ver_ok, ver_bad = 0, []
    for e, _, _ in ok:
        try:
            got = air(f"{BASE}/{PAPERS}/{e['id']}")
        except Exception as ex:
            ver_bad.append(e["id"]); print("  verify GET fail", e["id"], ex); continue
        f = got.get("fields", {})
        an = set(x["name"] if isinstance(x, dict) else x for x in f.get("Animal", []))
        af = set(x["name"] if isinstance(x, dict) else x for x in f.get("Afferent Types", []))
        if ((f.get("Notes") or "").strip()[:200] == e["notes"].strip()[:200]
                and an == set(e.get("animals") or [])
                and af == set(e.get("afferents") or [])):
            ver_ok += 1
        else:
            ver_bad.append(e["id"])
        time.sleep(0.2)
    print(f"echo verify: {ver_ok}/{len(ok)} exact")

    # twin creation (Models / Review Papers)
    twins = []
    for e, _, _ in ok:
        if e.get("is_model") and e.get("model_name"):
            twins.append(("tblsBq9IEv7dZe6fn", "fldfZB5ZwOaZKWFF8", e["model_name"], e["id"]))
        if e.get("is_review") and e.get("review_name"):
            twins.append(("tblSEubKcRId4wYMK", "fldF2F714aLZPexcv", e["review_name"], e["id"]))
    made = 0
    existing = {}
    for tid, _, _, _ in twins:
        if tid in existing:
            continue
        existing[tid] = {r["fields"].get("Name", ""): r["id"]
                        for r in air(f"{BASE}/{tid}", [("pageSize", 100)]).get("records", [])}
        off = None
        while True:
            d = air(f"{BASE}/{tid}", [("pageSize", 100)] + ([("offset", off)] if off else []))
            existing[tid].update({r["fields"].get("Name", ""): r["id"]
                                  for r in d.get("records", [])})
            off = d.get("offset")
            if not off:
                break
    for tid, fld, name, paper_id in twins:
        try:
            if name in existing[tid]:
                # link existing twin to this paper too (additive)
                cur = air(f"{BASE}/{tid}", [("recordIds[]", [existing[tid][name]]),
                                            ("fields[]", ["Paper"])]).get("records", [{}])[0]
                linked = [x["id"] if isinstance(x, dict) else x
                          for x in cur["fields"].get("Paper", [])]
                if paper_id not in linked:
                    linked.append(paper_id)
                    air_req(f"{BASE}/{tid}", {"records": [{"id": existing[tid][name],
                                                           "fields": {fld: linked}}]}, "PATCH")
            else:
                air_req(f"{BASE}/{tid}", {"records": [{"fields": {**{"Name": name},
                                                                  fld: [paper_id]}}]}, "POST")
                made += 1
        except Exception as ex:
            print("  TWIN FAIL", name, ex)
        time.sleep(0.35)

    proposals = sorted({p for _, _, ps in ok for p in ps})
    flags = sorted({(e["id"], e.get("flags")) for e, _, _ in ok if e.get("flags")})
    with open(os.path.join(HERE, "curation_log.csv"), "a", encoding="utf-8", newline="") as fh:
        import csv as _c
        w = _c.writer(fh)
        w.writerow([time.strftime("%Y-%m-%d"), "CURV2", path,
                    f"{len(ok)} papers; written {written} (fail {failed}); echo {ver_ok}/{len(ok)}; "
                    f"twins created {made}; proposals: {'; '.join(proposals) if proposals else 'none'}; "
                    f"flags: {flags if flags else 'none'}"])
    print(f"log appended | twins created {made} | proposals: {proposals or 'none'}")
    if ver_bad:
        print("VERIFY FAILURES:", ver_bad)


if __name__ == "__main__":
    main(sys.argv[1], dry="--dry" in sys.argv, verify_only="--verify-only" in sys.argv)
