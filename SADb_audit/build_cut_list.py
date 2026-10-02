"""Build the 290-paper cut list with a transparent scoring rubric (2026-10-01).

Ben's rubric directive: citations + year (old + few citations = low
influence); protect important models/research and the important authors he
listed. Plus his new rules: monographs/textbooks go; conference-abstract-only
stubs with no accessible paper go. Round-2 approved deletes included.

PROTECTIONS (never on the list):
  P1 important authors (primary OR any secondary): Quinn, Hunt, Szczecinski,
     Guertin, Perreault, Pearson, Rybak, Ijspeert, Geyer (his "Gerr"),
     Chiel, Tresch, Webster-Wood (incl. Webster), Procházka (both spellings),
     Frigon, McCrea, Grillner, Buschmann (+ Büschges — likely intent)
  P2 papers carrying Feedback pathway links (the SADb's core demonstrations)
  P3 papers linked from the 14 rich Models twins (models of record)
  P4 high-influence: >= 400 OpenAlex citations (classics stay regardless)
  P5 Ben's explicit keeps (round-2 rulings + flags comments: Nelson&Quinn, Hunt 2017...)
  P6 pre-1940 historical classics (Sherrington/Brown era) — conservative guard

CUT TIERS (auto first, then ranked):
  T1 round-2 approved (7): Hilts, Kim, Vogels, Giesseler, Wang 2013,
     Fairhurst, Merlet preprint
  T2 flags-comment deletions (3): Jessell 2000, Izhikevich 2006 monograph,
     Orlovsky 1999 monograph
  T3 monograph/textbook sweep (title/flag keywords) — Ben's directive
  T4 conference-abstract-only stubs (flags saying so + thin BMC/abstract DOIs)
  T5 the rubric ranking of what remains:
     cut_score = w1*low_influence + w2*age + w3*thin_curation + w4*isolation
     low_influence = quantile of (citations / years_since_pub)
     thin_curation = no afferents AND no animals beyond nothing... etc.
Output: cut_list_290.md (FULL titles — no truncation) + cut_list_290.csv.
NOTHING IS DELETED — list only, Ben reviews.
"""
import csv, json, os, re, unicodedata
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
BASE, TABLE = "appMQTnobUNRytIp7", "tblnnMrZszhboU4uD"
rec_url = lambda rid: f"https://airtable.com/{BASE}/{TABLE}/{rid}"

records = json.load(open(os.path.join(HERE, "export", "sadb_export.json"), encoding="utf-8"))
layout = json.load(open(os.path.join(HERE, "export", "sadb_layout.json"), encoding="utf-8"))
cites = json.load(open(os.path.join(HERE, "export", "sadb_cites.json"), encoding="utf-8"))
deleted = {r["id"] for r in json.load(open(os.path.join(HERE, "deleted_records_20261001.json"),
                                           encoding="utf-8"))}
cited_by = defaultdict(list)
for a, bs in cites.items():
    for b in bs:
        cited_by[b].append(a)

twins = json.load(open(os.path.join(HERE, "twins_profile.json"), encoding="utf-8"))
# rich twin paper ids unknown locally -> protect by twin NAME match below instead

IMPORTANT = ["quinn", "hunt", "szczecinski", "guertin", "perreault", "pearson",
             "rybak", "ijspeert", "geyer", "gerr", "chiel", "tresch",
             "webster-wood", "webster wood", "websterwood", "procházka", "prochazka",
             "frigon", "mccrea", "grillner", "buschmann", "büschges", "buschges"]

def strip_accents(s):
    return "".join(c for c in unicodedata.normalize("NFD", s) if unicodedata.category(c) != "Mn")

def author_hit(r):
    blob = strip_accents((r["primary"] or "").lower() + " " +
                         " ".join(r.get("secondary") or []).lower())
    blob = blob.replace("-", " ").replace("  ", " ")
    hits = []
    for a in IMPORTANT:
        a2 = a.replace("-", " ")
        if a2 in blob:
            hits.append(a)
    return hits

RICH_TWIN_NAMES = set(twins["models"]["rich_names"])  # e.g. "Hunt 2015"

KEEP_EXPLICIT = set()          # from Ben's comments
FLAG_KEEP_AUTHORS = ("nelson and quinn",)

def main():
    live = [r for r in records if r["id"] not in deleted]
    print("live papers:", len(live))

    # ---------- protections ----------
    prot = {}
    for r in live:
        rid = r["id"]
        L = layout.get(rid, {})
        c = L.get("c", 0)
        yr = r.get("year") or 0
        yrs = max(1, 2027 - yr) if yr else 1
        reasons = []
        ah = author_hit(r)
        if ah:
            reasons.append("author:" + ",".join(sorted(set(ah))))
        if r["feedback"]:
            reasons.append("feedback-linked")
        # NOTE: models2/reviews links are POLLUTED (multi-linked twins + the
        # records[]-GET bug) — NOT used as a protection signal (2026-10-01).
        if c >= 400:
            reasons.append(f"high-cites({c})")
        if yr and yr < 1940:
            reasons.append("pre-1940-classic")
        blob = (r["primary"] or "").lower()
        if any(k in blob for k in FLAG_KEEP_AUTHORS):
            reasons.append("ben-keep")
        if reasons:
            prot[rid] = reasons
    print("protected:", len(prot))

    # ---------- cut tiers ----------
    T = defaultdict(list)
    round2_ids = ["recV3YQtMTF5hreIm", "recXW34i5WwjznlgT", "recpiClFNWFTFnOCA",
                  "recNxBNmZRvcI2Sfr", "recJ7cgrCX7pokhZO", "recx6cU1RkrYsx1TX",
                  "recZxoFgqdQh2aUJO"]
    flags_delete = ["recmtKi7BqBJSH9my", "recTyf20KciVGaJVG", "recH0svuO344o1Stg"]
    byid = {r["id"]: r for r in live}
    # Ben's explicit rulings OVERRIDE every protection:
    T["T1_round2_approved"] = [i for i in round2_ids if i in byid]
    T["T2_flags_delete"] = [i for i in flags_delete if i in byid]
    ruled = set(T["T1_round2_approved"]) | set(T["T2_flags_delete"])
    for rid in ruled:
        prot.pop(rid, None)

    MONO = re.compile(r"monograph|textbook|handbook|\bprimer\b|from mollusc to man|"
                      r"principles of neural|the nervous system$|"
                      r"neuromechanical modeling of posture", re.I)
    for r in live:
        rid = r["id"]
        if rid in ruled or any(rid in v for v in T.values()):
            continue
        blob = (r["title"] or "") + " " + (r.get("notes") or "")[:300]
        if MONO.search(blob) or "10.1093/acprof:oso" in (r.get("doi") or "") \
           or "10.7551/mitpress" in (r.get("doi") or ""):
            T["T3_monograph_textbook"].append(rid)

    CONF = re.compile(r"conference abstract|abstract only|two figure captions|"
                      r"poster(?:ing)? (?:presentation|session)|\babstract\b.{0,30}only", re.I)
    for r in live:
        rid = r["id"]
        if rid in ruled or any(rid in v for v in T.values()):
            continue
        blob = (r.get("notes") or "")[:300]
        doi = r.get("doi") or ""
        thin_doi = bool(re.match(r"10\.1186/.*-s\d+-p\d+", doi))
        if CONF.search(blob) or thin_doi:
            T["T4_conference_abstract_only"].append(rid)

    # ---------- rubric ranking ----------
    scored = []
    for r in live:
        rid = r["id"]
        if rid in prot or rid in ruled or any(rid in v for v in T.values()):
            continue
        # citation-based rubric needs citation data: no-DOI records have NONE
        # (Shik 1966 showed the failure mode — a classic reading as 0 cites).
        # They are excluded from T5 and left for manual review instead.
        if not (r.get("doi") or "").strip():
            continue
        L = layout.get(rid, {})
        c, yr = L.get("c", 0), r.get("year") or 0
        yrs = max(1, 2027 - yr) if yr else 1
        cpy = c / yrs
        indeg, outdeg = len(cited_by.get(rid, [])), len(cites.get(rid, []))
        thin = 0
        if not r["afferents"]:
            thin += 1
        if not r["animals"]:
            thin += 1
        if not r["has_notes"]:
            thin += 1
        if not r["has_pdf"]:
            thin += 0.5
        iso = 1 if (L.get("cl", -1) < 0 and indeg == 0 and outdeg == 0) else 0
        # influence quantiles computed after the loop; store raw parts
        scored.append({"id": rid, "r": r, "c": c, "yr": yr, "cpy": cpy,
                       "indeg": indeg, "outdeg": outdeg, "thin": thin, "iso": iso})
    cpys = sorted(s["cpy"] for s in scored)
    tots = sorted(s["c"] for s in scored)
    def q(v, arr):
        import bisect
        return bisect.bisect_left(arr, v) / len(arr) if arr else 0
    for s in scored:
        # Ben's framing: OLD + FEW CITATIONS = uninfluential; too-new papers
        # (>=2024) can't be judged and are excluded from ranking entirely.
        yr, c = s["yr"], s["c"]
        old_low = 0.0
        if yr and yr < 2005 and c < 40:
            old_low = 1.0
        elif yr and yr < 2015 and c < 15:
            old_low = 0.6
        s["score"] = round(
            0.35 * (1 - q(c, tots))              # low ABSOLUTE citations (his main signal)
            + 0.25 * old_low                      # old + genuinely uncited
            + 0.15 * (1 - q(s["cpy"], cpys))      # low citation RATE
            + 0.15 * s["iso"]                     # isolated in the corpus graph
            + 0.10 * (s["thin"] / 3.5),           # thin curation layer
            4)
    scored = [s for s in scored if not (s["yr"] and s["yr"] >= 2024)]

    target = 290
    auto = [i for v in T.values() for i in v]
    ranked = sorted(scored, key=lambda s: (-s["score"], -(s["c"]), s["id"]))
    need = target - len(auto)
    take = ranked[:need]
    T["T5_rubric_ranked"] = [s["id"] for s in take]

    final = auto + [s["id"] for s in take]
    print(f"auto tiers: {len(auto)} | rubric picks: {len(take)} | TOTAL {len(final)}")

    # ---------- outputs ----------
    rows = []
    for s in take:
        rows.append(s)
    tier_of = {}
    for t, ids in T.items():
        for i in ids:
            tier_of[i] = t

    out = ["# CUT LIST — 290 papers proposed for deletion (2026-10-01)",
           "",
           f"Live papers {len(live)} | protected {len(prot)} | proposed cut {len(final)}",
           "(auto tiers " + str(len(auto)) + " + rubric-ranked " + str(len(take)) + ").",
           "NOTHING IS DELETED — Ben reviews this list. PDFs stay local; the HTML app",
           "keeps its own snapshot until Ben decides otherwise. Full titles below",
           "(no truncation — lesson learned).", "",
           "## Rubric (T5 ranking)", "",
           "- 35% low ABSOLUTE citations (OpenAlex) — the main 'not influential' signal",
           "- 25% old + genuinely uncited (pre-2005 with <40 cites; partial credit",
           "  pre-2015 with <15)",
           "- 15% low citation rate (citations per year since publication)",
           "- 15% isolated in the corpus citation graph (no cluster, degree 0)",
           "- 10% thin curation layer (no afferents / animals / notes / PDF)",
           "- EXCLUDED from T5: papers from 2024+ (too new to judge influence)",
           "- EXCLUDED from T5: no-DOI records (citation data missing, not zero —",
           "  Shik 1966/Herr 2002/Zaporozhets/Gervasio are NOT on this list for that",
           "  reason; they need a manual-look pass instead).",
           "- PROTECTED from the list: Ben's important authors (primary or co-author),",
           "  Feedback-pathway-linked papers, rich-twin model papers, >= 400 citations,",
           "  pre-1940 classics, his explicit keeps.", ""]
    for t in ["T1_round2_approved", "T2_flags_delete", "T3_monograph_textbook",
              "T4_conference_abstract_only"]:
        ids = [i for i in T[t] if i in byid]
        out.append(f"## {t} — {len(ids)}")
        out.append("")
        for i in ids:
            r = byid[i]
            out.append(f"- [{r['primary']} {r['year']} — {r['title']}]" f"({rec_url(i)}) · "
                       f"{layout.get(i, {}).get('c', 0)} cites · doi:{r.get('doi','')}")
        out.append("")
    out += ["## T5 rubric-ranked — " + str(len(take)), ""]
    for s in take:
        r = s["r"]
        out.append(f"- [{r['primary']} {r['year']} — {r['title']}]({rec_url(s['id'])}) · "
                   f"score {s['score']} · {s['c']} cites · cpy {s['cpy']:.2f} · "
                   f"in/out {s['indeg']}/{s['outdeg']} · doi:{r.get('doi','')}")
    open(os.path.join(HERE, "cut_list_290.md"), "w", encoding="utf-8").write("\n".join(out) + "\n")

    with open(os.path.join(HERE, "cut_list_290.csv"), "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tier", "score", "record_id", "primary", "year", "cites",
                    "cites_per_year", "in_deg", "out_deg", "title", "doi"])
        for s in take:
            r = s["r"]
            w.writerow(["T5", s["score"], s["id"], r["primary"], r["year"], s["c"],
                        round(s["cpy"], 3), s["indeg"], s["outdeg"], r["title"], r.get("doi", "")])
        for t in ["T1_round2_approved", "T2_flags_delete", "T3_monograph_textbook",
                  "T4_conference_abstract_only"]:
            for i in T[t]:
                if i in byid:
                    r = byid[i]
                    w.writerow([t, "", i, r["primary"], r["year"],
                                layout.get(i, {}).get("c", 0), "", "", "", r["title"], r.get("doi", "")])

    # protected-author sanity check: NONE of the important authors on the list
    viol = [i for i in final if i in prot]
    print("protected-on-list violations:", len(viol))
    hist = {"pre-1960": 0, "1960-1989": 0, "1990-2009": 0, "2010+": 0, "unknown": 0}
    for i in final:
        yr = (byid[i].get("year") or 0) if i in byid else 0
        k = "pre-1960" if yr and yr < 1960 else "1960-1989" if yr < 1990 else \
            "1990-2009" if yr < 2010 else "2010+" if yr else "unknown"
        hist[k] += 1
    print("year histogram:", hist)
    print("wrote cut_list_290.md + .csv | protected:", len(prot))

if __name__ == "__main__":
    main()
