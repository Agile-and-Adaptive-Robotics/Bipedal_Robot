#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SADb monthly watch — the skimmer-curator sweep, generalizing the proven
rr-digest pattern (SADb_audit\\rr_digest_20260912.py + rr_digest_2023plus.md).

One run per month: sweep OpenAlex + Crossref (both open, no key; Scholar/WoS
scraping is deliberately NOT done — standing constraint) for NEW papers in the
lookback window that match Ben's watch authors or corpus-vocabulary keywords,
dedupe against the SADb corpus DOI set and every previously routed DOI, collapse
preprint/version twins, Crossref-verify every shortlist entry, then write:

  <out_root>\\reports\\YYYY-MM.md            the digest table (DOI/year/first author/short title/trigger/why)
  <out_root>\\reports\\YYYY-MM_shortlist.json  batch3 import pattern (title/authors/container/year/doi,
                                              first_author/trigger/why/fetch) — the same row shape as
                                              SADb_audit\\batch3\\rr_digest_shortlist.json, so
                                              SADb_audit\\zotero_add_digest.py and the Airtable
                                              Name/Author/Year/DOI import flow keep working unchanged.
  <out_root>\\state.json                     routed DOIs + per-month stored entries (idempotency).

IDEMPOTENT: state.json records every routed DOI and each month's stored report
rows. Re-running a month rebuilds the report from (stored rows + genuinely new
hits), so rows never duplicate, and a re-run that finds nothing new rewrites the
same bytes. Each month keeps the window of its first run unless --window is
given. An empty shortlist is a valid result and says so in the report.

NOTHING is ever written to Airtable or Zotero by this script — the machine
drafts, the human rules. Run with the myo env python by full path.
"""
import argparse
import json
import os
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from datetime import date, timedelta

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
UA = "SADbMonthlyWatch/1.0 (SADb_audit watch; personal research digest; python-urllib)"
REQUEST_TIMEOUT = 30
DELAY = 0.5          # seconds between API calls (be a good citizen)
OA_PAGE = 50         # small-ish page size
CR_ROWS = 50
OA_KW_PAGES = 4
KW_GROUP = 1         # keywords per OpenAlex query — larger pipe-OR groups get
                     # HTTP 401 on the anonymous (no-key) tier (observed 2026-10-23)

SCHEDULER_LINE = (
    'schtasks /Create /SC MONTHLY /D 1 /ST 09:00 /TN SADbMonthlyWatch '
    '/TR "\\"C:\\Users\\Ben Bolen\\.conda\\envs\\myo\\python.exe\\" '
    'D:\\Github\\Bipedal_Robot\\SADb_audit\\watch\\monthly_watch.py"'
)


# ---------------------------------------------------------------- helpers ---
def die(msg):
    print("ERROR: " + msg)
    sys.exit(1)


def norm_doi(doi):
    d = (doi or "").strip().lower()
    for pre in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "doi:"):
        if d.startswith(pre):
            d = d[len(pre):]
    return d


def norm_title(t):
    return re.sub(r"[^a-z0-9]+", "", (t or "").lower())


def strip_accents(s):
    return "".join(c for c in unicodedata.normalize("NFD", s) if not unicodedata.combining(c))


def family_of(name):
    """Family name from 'Family Given' or 'Given, Family' raw author strings."""
    n = (name or "").strip()
    if not n:
        return ""
    if "," in n:
        return n.split(",", 1)[0].strip()
    return n.split()[-1].strip() if n.split() else ""


def short_title(t, width=58):
    t = " ".join((t or "").split())
    return t if len(t) <= width else t[: width - 1].rstrip() + "…"


def strip_doi_version(d):
    """Version-suffix strip for twin detection: eLife's '.2'-style versioned
    DOIs and Research Square '/v1' DOIs collapse onto their base DOI. Used only
    with a title-compatibility guard, never as a bare merge key."""
    d = re.sub(r"/v\d+$", "", d)
    d = re.sub(r"\.\d+$", "", d)
    return d


def titles_compatible(a, b):
    a, b = a or "", b or ""
    if not a or not b:
        return False
    return a == b or a.startswith(b) or b.startswith(a)


# ------------------------------------------------------------------- HTTP ---
def http_json(url, mailto=""):
    """GET one JSON doc. OpenAlex/Crossref only; backoff retries on 429/5xx and
    on 401 — OpenAlex's anonymous (no-key) tier shares a small daily budget per
    IP and intermittently refuses with 401; a pause usually clears it."""
    sep = "&" if "?" in url else "?"
    if mailto:
        url = url + sep + "mailto=" + urllib.parse.quote(mailto)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    last = None
    for attempt, backoff in enumerate((0, 10, 30), start=1):
        try:
            with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as r:
                return json.loads(r.read().decode("utf-8", "replace"))
        except urllib.error.HTTPError as e:
            last = "HTTP %d %s" % (e.code, url[:120])
            if e.code not in (401, 429, 500, 502, 503) or attempt == 3:
                raise RuntimeError(last)
        except (urllib.error.URLError, OSError, ValueError) as e:
            last = "%s %s" % (type(e).__name__, url[:120])
            if attempt == 3:
                raise RuntimeError(last)
        time.sleep(backoff)
    raise RuntimeError(last or "request failed")


# ------------------------------------------------------------- corpus/state --
def load_corpus_dois(export_path):
    """DOIs of ALL corpus papers (retained + archived), papers only — never the
    non-paper records (models/reviews/feedback/recorders)."""
    with open(export_path, encoding="utf-8") as f:
        records = json.load(f)
    dois = set()
    for r in records:
        if (r.get("recordType") or "") != "paper":
            continue
        d = norm_doi(r.get("doi"))
        if d:
            dois.add(d)
    return dois, len([r for r in records if r.get("recordType") == "paper"])


def load_state(path):
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return {"schema": 1, "created": date.today().isoformat(), "routed_dois": {},
            "months": {}, "run_log": []}


# ------------------------------------------------------------- candidates ---
def oa_abstract(work):
    inv = work.get("abstract_inverted_index")
    if not inv:
        return ""
    pos = {}
    for word, idxs in inv.items():
        for i in idxs:
            pos[i] = word
    return " ".join(pos[i] for i in sorted(pos))


def from_oa(w):
    names = [((a.get("author") or {}).get("display_name") or "")
             for a in (w.get("authorships") or [])]
    return {
        "doi": norm_doi(w.get("doi")),
        "title": (w.get("display_name") or "").strip(),
        "year": w.get("publication_year") or "",
        "pubdate": w.get("publication_date") or "",
        "names": names,
        "families": [family_of(n) for n in names if family_of(n)],
        "abstract": oa_abstract(w),
        "container": (((w.get("primary_location") or {}).get("source") or {}).get("display_name") or ""),
        "type": w.get("type") or "",
    }


def from_cr(it):
    dp = ((it.get("issued") or {}).get("date-parts") or [[None]])[0]
    y = dp[0] if dp else ""
    names = it.get("author") or []
    pd = ""
    for key in ("published-online", "published-print", "issued", "created"):
        parts = (it.get(key) or {}).get("date-parts") or []
        if parts and parts[0] and len(parts[0]) >= 3:
            pd = "%04d-%02d-%02d" % tuple(parts[0][:3])
            break
    return {
        "doi": norm_doi(it.get("DOI")),
        "title": " ".join((it.get("title") or [""])[0].split()) if it.get("title") else "",
        "year": y,
        "pubdate": pd,
        "names": [" ".join(x for x in ((a.get("family") or ""), (a.get("given") or "")) if x).strip()
                  for a in names],
        "families": [(a.get("family") or "").strip() for a in names if (a.get("family") or "").strip()],
        "abstract": re.sub(r"<[^>]+>", " ", it.get("abstract") or ""),
        "container": (it.get("container-title") or [""])[0] if it.get("container-title") else "",
        "type": it.get("type") or "",
    }


def collect_candidates(cfg, w0, w1, log):
    """One OpenAlex author OR-query + grouped keyword OR-queries; per-term
    Crossref queries. Returns (candidates_by_doi, raw_hits, notes)."""
    cand = {}
    raw = [0]
    notes = []

    def add(rec, src):
        d = rec["doi"]
        if not d:
            return
        raw[0] += 1
        if d in cand:
            c = cand[d]
            for k in ("title", "abstract", "names", "families", "container", "pubdate"):
                if not c.get(k) and rec.get(k):
                    c[k] = rec[k]
            if src not in c["sources"]:
                c["sources"].append(src)
        else:
            rec["sources"] = [src]
            cand[d] = rec

    mailto = cfg.get("mailto") or ""
    oa_select = ("&select=doi,display_name,publication_year,publication_date,"
                 "authorships,primary_location,type,abstract_inverted_index")

    def safe(url):
        """One query; a persistently refused call is skipped and logged, never
        fatal — a partial sweep still produces an honest report."""
        try:
            return http_json(url, mailto)
        except RuntimeError as e:
            notes.append("query skipped after retries: %s" % str(e)[:140])
            log("query skipped after retries: %s" % str(e)[:140])
            return None

    # OpenAlex's anonymous (no-key) tier shares a small per-IP daily budget
    # (resets midnight UTC). Once it refuses with 429/401 after retries it
    # stays refused for hours — so the rest of the run's OpenAlex queries are
    # skipped at once (one logged note) and Crossref carries the sweep. The
    # state-dedupe design means a later full run simply folds in what this one
    # missed.
    oa = {"dead": False, "skipped": 0}

    def oa_safe(url):
        if oa["dead"]:
            oa["skipped"] += 1
            return None
        try:
            return http_json(url, mailto)
        except RuntimeError as e:
            msg = str(e)
            if "HTTP 429" in msg or "HTTP 401" in msg:
                oa["dead"] = True
                note = ("OpenAlex unavailable this run (HTTP 429/401 after retries — anonymous "
                        "per-IP daily budget exhausted; resets midnight UTC) — remaining OpenAlex "
                        "queries skipped, Crossref carried the sweep")
                notes.append(note)
                log(note)
            else:
                notes.append("query skipped after retries: %s" % msg[:140])
                log("query skipped after retries: %s" % msg[:140])
            return None

    # --- OpenAlex: one gated query per author -------------------------------
    # raw_author_name.search alone matches ANY name-holder (common surnames
    # like Zhang/Deng/Young return tens of thousands of works), so the domain
    # gate rides IN the filter: author surname AND a domain term in
    # title+abstract. The local trigger re-check still runs on every candidate.
    gate = "|".join(cfg.get("domain_gate") or [])
    authors = cfg["authors"]
    for a in authors:
        filt = "raw_author_name.search:" + urllib.parse.quote(a)
        if gate:
            filt += ",title_and_abstract.search:" + urllib.parse.quote(gate)
        url = ("https://api.openalex.org/works?filter=" + filt
               + ",from_publication_date:%s,to_publication_date:%s&per-page=%d" % (w0, w1, OA_PAGE)
               + oa_select)
        d = oa_safe(url)
        if d is None:
            continue
        total = (d.get("meta") or {}).get("count", 0)
        for w in d.get("results", []):
            add(from_oa(w), "openalex")
        if total > OA_PAGE:
            notes.append("OpenAlex author %s truncated at %d/%d" % (a, OA_PAGE, total))
        log("OpenAlex author %s: %d gated hits" % (a, total))
        time.sleep(DELAY)

    # --- OpenAlex: keywords in pipe-OR groups -------------------------------
    kws = cfg["keywords"]
    for gi in range(0, len(kws), KW_GROUP):
        group = kws[gi: gi + KW_GROUP]
        url = ("https://api.openalex.org/works?filter=title_and_abstract.search:"
               + urllib.parse.quote("|".join(group))
               + ",from_publication_date:%s,to_publication_date:%s&per-page=%d" % (w0, w1, OA_PAGE)
               + oa_select + "&cursor=")
        cursor = "*"
        pages = 0
        total = None
        while cursor and pages < OA_KW_PAGES:
            d = oa_safe(url + urllib.parse.quote(cursor))
            if d is None:
                break
            total = (d.get("meta") or {}).get("count", total)
            for w in d.get("results", []):
                add(from_oa(w), "openalex")
            cursor = (d.get("meta") or {}).get("next_cursor")
            pages += 1
            time.sleep(DELAY)
        if total and pages * OA_PAGE < total:
            notes.append("OpenAlex keyword group %d truncated at %d/%d" % (gi // KW_GROUP + 1, pages * OA_PAGE, total))
        log("OpenAlex keywords %s: %d hits" % ("|".join(group), total or 0))

    # --- Crossref: one query per author, one per keyword --------------------
    cr_base = "https://api.crossref.org/works?filter=from-pub-date:%s,until-pub-date:%s&rows=%d&select=DOI,title,author,issued,container-title,type,abstract" % (w0, w1, CR_ROWS)
    for a in authors:
        d = safe(cr_base + "&query.author=" + urllib.parse.quote(a))
        if d is not None:
            for it in (d.get("message") or {}).get("items", []):
                add(from_cr(it), "crossref")
            log("Crossref author %s: %d hits" % (a, (d.get("message") or {}).get("total-results", 0)))
        time.sleep(DELAY)
    for k in kws:
        d = safe(cr_base + "&query.bibliographic=" + urllib.parse.quote(k))
        if d is not None:
            for it in (d.get("message") or {}).get("items", []):
                add(from_cr(it), "crossref")
            log("Crossref keyword %s: %d hits" % (k, (d.get("message") or {}).get("total-results", 0)))
        time.sleep(DELAY)

    return cand, raw[0], notes


def in_window(rec, w0, w1):
    pd = rec.get("pubdate") or ""
    if len(pd) >= 10 and pd[:10] <= w1 and pd[:10] >= w0:
        return True
    y = rec.get("year")
    try:
        y = int(y)
    except (TypeError, ValueError):
        return False
    return w0[:4] <= str(y) <= w1[:4] if len(w0) >= 4 else False


# ----------------------------------------------------------- trigger logic --
def triggers_for(rec, cfg):
    """(authors_fired, keywords_fired). An author hit requires the family name
    to match a watch author AND a domain-gate term in title+abstract — surname
    matches alone bring in unrelated name-holders (a wireless 'J.-F. Frigon',
    Zenodo deposits by 'Marlène Frigon'), so the gate mechanizes the rr digest's
    scope judgment. Keyword hits are exact case-insensitive substrings."""
    text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
    text_low = strip_accents(text).lower()

    fams = set()
    for f in rec.get("families") or []:
        f = strip_accents(f).lower()
        if f:
            fams.add(f)
    if not fams:  # records without structured families: last-token fallback
        for n in rec.get("names") or []:
            f = strip_accents(family_of(n)).lower()
            if f:
                fams.add(f)
    a_fired = sorted({a for a in cfg["_authors_norm"] if a in fams})

    k_fired = []
    for k in cfg["_keywords_low"]:
        if k in text_low:
            k_fired.append(k)
    k_fired = sorted(set(k_fired))

    if a_fired and not any(g in text_low for g in cfg["_gate_low"]):
        a_fired = []  # scope gate: not the lab's domain
    return a_fired, k_fired


# ------------------------------------------------------------------- main ---
def main():
    t0 = time.time()
    ap = argparse.ArgumentParser(description="SADb monthly watch (OpenAlex + Crossref skimmer)")
    ap.add_argument("--config", default=os.path.join(SCRIPT_DIR, "watch_config.json"))
    ap.add_argument("--out", default=SCRIPT_DIR,
                    help="output ROOT override: reports\\ and state.json live under this dir")
    ap.add_argument("--month", default=None, help="report month YYYY-MM (window = previous month start through month end/today)")
    ap.add_argument("--window", nargs=2, metavar=("START", "END"), default=None,
                    help="explicit window YYYY-MM-DD YYYY-MM-DD (report month = END's month)")
    args = ap.parse_args()

    with open(args.config, encoding="utf-8") as f:
        cfg = json.load(f)
    cfg["_authors_norm"] = {strip_accents(a).lower() for a in cfg.get("authors") or []}
    cfg["_keywords_low"] = [k.lower() for k in cfg.get("keywords") or []]
    cfg["_gate_low"] = [g.lower() for g in cfg.get("domain_gate") or []]
    if not cfg["authors"] or not cfg["keywords"]:
        die("config must define non-empty 'authors' and 'keywords'")

    export_path = cfg.get("corpus_export") or os.path.join(SCRIPT_DIR, "..", "export", "sadb_export.json")
    if not os.path.isabs(export_path):
        export_path = os.path.normpath(os.path.join(SCRIPT_DIR, export_path))
    corpus_dois, n_papers = load_corpus_dois(export_path)

    reports_dir = os.path.join(args.out, "reports")
    os.makedirs(reports_dir, exist_ok=True)
    state_path = os.path.join(args.out, "state.json")
    state = load_state(state_path)
    routed = state.setdefault("routed_dois", {})
    # twin lookups over everything already routed (any month): a fresh hit that
    # is a version twin of a routed work must be dropped, or a re-run would add
    # eLife-v2-style rows under DOIs the first run never saw.
    routed_nt = set()
    routed_sv = {}
    for rdoi, rmeta in routed.items():
        nt = norm_title(rmeta.get("title") or "")
        if nt:
            routed_nt.add(nt)
        routed_sv.setdefault(strip_doi_version(rdoi), set()).add(nt)

    today = date.today()
    if args.window:
        w0, w1 = args.window
        month_key = w1[:7]
    else:
        month_key = args.month or today.strftime("%Y-%m")
        y, m = int(month_key[:4]), int(month_key[5:7])
        sy, sm = (y - 1, 12) if m == 1 else (y, m - 1)
        w0 = "%04d-%02d-01" % (sy, sm)
        ey = y + (1 if m == 12 else 0)
        em = 1 if m == 12 else m + 1
        last_day = date(ey, em, 1) - timedelta(days=1)
        w1 = str(min(last_day, today))
    # window-sticky months: a re-run reuses the month's first-run window
    stored_window = (state["months"].get(month_key) or {}).get("window")
    if stored_window and not args.window:
        w0, w1 = stored_window

    print("== SADb monthly watch %s ==" % month_key)
    print("window: %s .. %s   corpus papers: %d   corpus DOIs: %d" % (w0, w1, n_papers, len(corpus_dois)))

    def log(msg):
        print("  " + msg)

    cand, raw_hits, notes = collect_candidates(cfg, w0, w1, log)
    log("raw candidate records pulled (both sources, before dedupe): %d" % raw_hits)
    log("distinct candidate DOIs merged: %d" % len(cand))

    # window clamp
    cand = {d: r for d, r in cand.items() if in_window(r, w0, w1)}
    log("in-window candidates: %d" % len(cand))

    # corpus + routed dedupe, then twin-of-routed drop
    in_corpus = sorted(d for d in cand if d in corpus_dois)
    already = sorted(d for d in cand if d in routed)
    fresh = {d: r for d, r in cand.items() if d not in corpus_dois and d not in routed}
    twin_routed = []
    kept = {}
    for d, r in fresh.items():
        nt = norm_title(r.get("title"))
        sv = strip_doi_version(d)
        if (nt and nt in routed_nt) or (
                sv in routed_sv and any(titles_compatible(nt, t) for t in routed_sv[sv])):
            twin_routed.append(d)
            continue
        kept[d] = r
    fresh = kept
    log("dropped, already in corpus (%d papers): %d" % (n_papers, len(in_corpus)))
    log("dropped, routed in an earlier run/month: %d" % len(already))
    log("dropped, version twin of an already-routed work: %d" % len(twin_routed))

    # trigger filter
    hits = []
    gated_out = 0
    for d, r in sorted(fresh.items()):
        a, k = triggers_for(r, cfg)
        if k:
            hits.append((d, r, a, k))
        elif a:
            gated_out += 1  # author matched but no domain term: out of scope
    log("trigger hits (author/keyword): %d; author hits gated out of scope: %d" % (len(hits), gated_out))

    # twin collapse (preprint / version twins share a normalized title)
    groups = {}
    for d, r, a, k in hits:
        groups.setdefault(norm_title(r.get("title")), []).append((d, r, a, k))
    survivors = []
    twins_dropped = []
    for key, grp in groups.items():
        if len(grp) == 1:
            survivors.append(grp[0])
            continue

        def rank(t):
            d, r, a, k = t
            type_score = 4 if r.get("type") in ("journal-article", "article") else \
                         3 if r.get("type") in ("proceedings-article",) else 1
            return (-type_score, 0 if r.get("container") else 1, d)
        grp = sorted(grp, key=rank)
        survivors.append(grp[0])
        twins_dropped.extend((d, grp[0][0]) for d, _, _, _ in grp[1:])
    log("version twins collapsed (same title): %d DOI(s) dropped in favor of the published version" % len(twins_dropped))

    # twin collapse stage 2: version-suffix twins (e.g. eLife's '.2' DOIs) whose
    # title text drifted slightly — same version-stripped DOI + compatible title
    by_sv = {}
    for t in survivors:
        by_sv.setdefault(strip_doi_version(t[0]), []).append(t)
    survivors2 = []
    for sv, grp in by_sv.items():
        if len(grp) == 1:
            survivors2.append(grp[0])
            continue

        def rank(t):
            d, r, a_f, k_f = t
            type_score = 4 if r.get("type") in ("journal-article", "article") else \
                         3 if r.get("type") in ("proceedings-article",) else 1
            return (-type_score, 0 if r.get("container") else 1, d)
        acc = []
        for t in sorted(grp, key=rank):
            if any(titles_compatible(norm_title(t[1].get("title")), norm_title(a[1].get("title")))
                   for a in acc):
                twins_dropped.append((t[0], acc[0][0]))
            else:
                acc.append(t)
        survivors2.extend(acc)
    if len(survivors2) != len(survivors):
        log("version twins collapsed (version-suffixed DOIs, compatible titles): %d DOI(s) dropped"
            % (len(survivors) - len(survivors2)))
    survivors = survivors2
    log("twin DOIs dropped: %s" % (", ".join(sorted(d for d, _ in twins_dropped)) or "none"))

    # Crossref verification of EVERY shortlist entry
    entries = []
    verify_failed = []
    for d, r, a_f, k_f in sorted(survivors, key=lambda t: t[0]):
        try:
            cr = http_json("https://api.crossref.org/works/" + urllib.parse.quote(d), cfg.get("mailto") or "")
            m = cr.get("message") or {}
        except RuntimeError as e:
            verify_failed.append((d, str(e)[:100]))
            continue
        title = " ".join((m.get("title") or [r.get("title") or ""])[0].split())
        container = (m.get("container-title") or [r.get("container") or ""])[0]
        dp = ((m.get("issued") or {}).get("date-parts") or [[r.get("year")]])[0]
        year = dp[0] or r.get("year") or ""
        auth = m.get("author") or []
        authors_str = "; ".join(" ".join(x for x in ((au.get("family") or ""), (au.get("given") or "")) if x).strip()
                                for au in auth[:15])
        first_author = (auth[0].get("family") if auth and auth[0].get("family")
                        else family_of(r["names"][0]) if r.get("names") else "")
        trigger = []
        if a_f:
            trigger.append("author")
        if k_f:
            trigger.append("keyword")
        parts = []
        if a_f:
            parts.append("author hit: " + ", ".join(a_f))
        if k_f:
            parts.append("keyword hit: " + ", ".join(k_f))
        why = "; ".join(parts)
        entries.append({
            "doi": d, "title": title, "container": container, "year": year,
            "authors": authors_str, "first_author": first_author or "?",
            "trigger": "+".join(trigger) or "?", "why": why,
            "short_title": short_title(title),
        })
        time.sleep(DELAY)
    if verify_failed:
        log("Crossref verification FAILED for %d DOI(s) — dropped, not routed:" % len(verify_failed))
        for d, e in verify_failed:
            log("    %s  (%s)" % (d, e))

    # merge with the month's stored rows, sort deterministically
    month = state["months"].setdefault(month_key, {"window": [w0, w1], "entries": []})
    month["window"] = [w0, w1]
    by_doi = {e["doi"]: e for e in month["entries"]}
    for e in entries:
        by_doi[e["doi"]] = e
    all_rows = sorted(by_doi.values(), key=lambda e: (
        -int(e["year"]) if str(e["year"]).isdigit() else 0,
        (e["first_author"] or "").casefold(), e["doi"]))

    # ---- write outputs (report bytes depend only on the stored+new rows and
    # the run notes; same conditions → identical bytes) ----
    md = render_md(month_key, w0, w1, all_rows, cfg, n_papers, len(corpus_dois), notes)
    md_path = os.path.join(reports_dir, month_key + ".md")
    with open(md_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(md)
    rows_json = [{"doi": e["doi"], "why": e["why"], "title": e["title"], "container": e["container"],
                  "year": e["year"], "authors": e["authors"], "first_author": e["first_author"],
                  "trigger": e["trigger"], "fetch": "ok"} for e in all_rows]
    json_path = os.path.join(reports_dir, month_key + "_shortlist.json")
    with open(json_path, "w", encoding="utf-8", newline="\n") as f:
        json.dump(rows_json, f, indent=1, ensure_ascii=False)
        f.write("\n")

    # ---- commit state last ----
    for e in all_rows:
        routed[e["doi"]] = {"month": month_key, "title": e["title"]}
    month["entries"] = all_rows
    state.setdefault("run_log", []).append({
        "date": today.isoformat(), "month": month_key, "window": [w0, w1],
        "raw_hits": raw_hits, "candidates_distinct": len(cand), "in_corpus": len(in_corpus),
        "routed_already": len(already), "routed_twin_dropped": len(twin_routed),
        "trigger_hits": len(hits),
        "twins_dropped": len(twins_dropped), "verify_failed": len(verify_failed),
        "new_entries": len(entries), "shortlist_total": len(all_rows),
        "notes": notes,
    })
    state["run_log"] = state["run_log"][-60:]
    with open(state_path, "w", encoding="utf-8", newline="\n") as f:
        json.dump(state, f, indent=1, ensure_ascii=False)
        f.write("\n")

    print("report : %s (%d rows)" % (md_path, len(all_rows)))
    print("shortlist: %s" % json_path)
    for n in notes:
        print("note: " + n)
    print("elapsed: %.1f s" % (time.time() - t0))


def render_md(month_key, w0, w1, rows, cfg, n_papers, n_dois, notes=()):
    L = []
    L.append("# SADb monthly watch — %s" % month_key)
    L.append("")
    L.append("Window: **%s → %s** (previous-month-start through the month's run; one-month "
             "lookback plus small overlap; a month keeps the window of its first run)." % (w0, w1))
    L.append("")
    L.append("Method (the rr-digest pipeline, generalized — see `watch\\README.md`): OpenAlex + "
             "Crossref only (open APIs, no key; no Scholar/WoS scraping). Watch list: %d authors "
             "(top corpus primary authors by retained-paper count ∪ pinned core) and %d corpus-vocabulary "
             "keywords. Candidates deduped against the %d-paper corpus DOI set (%d DOIs, retained and "
             "archived) and every previously routed DOI (`watch\\state.json`); preprint/version twins "
             "collapsed; **every shortlist entry verified via Crossref**. Author hits carry a domain gate "
             "(`domain_gate` in `watch_config.json`) — surname matches alone pull in unrelated "
             "name-holders. Trigger fired per row in the table below."
             % (len(cfg.get("authors") or []), len(cfg.get("keywords") or []), n_papers, n_dois))
    L.append("")
    L.append("Shortlist JSON (Airtable Name/Author/Year/DOI import pattern): `reports\\%s_shortlist.json` — "
             "the batch3 row shape, so `SADb_audit\\zotero_add_digest.py` and the Airtable import flow keep "
             "working unchanged." % month_key)
    L.append("")
    L.append("**Nothing is written to Airtable or Zotero automatically. The machine drafts; the human rules.**")
    L.append("")
    L.append("## Shortlist (%d)" % len(rows))
    L.append("")
    if not rows:
        L.append("No papers met the watch triggers in this window. An empty shortlist is a valid result.")
        L.append("")
    else:
        L.append("| # | DOI | Year | First author | Title (short) | Trigger | Why relevant |")
        L.append("|---|-----|------|--------------|---------------|---------|--------------|")
        for i, e in enumerate(rows, 1):
            L.append("| %d | %s | %s | %s | %s | %s | %s |" % (
                i, e["doi"], e["year"], e["first_author"], e["short_title"], e["trigger"], e["why"]))
        L.append("")
    L.append("## Run notes")
    L.append("")
    if notes:
        for n in notes:
            L.append("- " + n)
    else:
        L.append("None — every query in the sweep completed.")
    L.append("")
    L.append("## Scheduling (provided, not installed)")
    L.append("")
    L.append(r'Paste in **cmd** (the `\"` escaping is required so the inner quotes survive):')
    L.append("")
    L.append("```")
    L.append(SCHEDULER_LINE)
    L.append("```")
    L.append("")
    L.append("Or run it by hand any time: `\"C:\\Users\\Ben Bolen\\.conda\\envs\\myo\\python.exe\" "
             "D:\\Github\\Bipedal_Robot\\SADb_audit\\watch\\monthly_watch.py`")
    L.append("")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
