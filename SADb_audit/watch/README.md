# SADb monthly watch — the skimmer–curator workflow

Ben's monthly literature sweep, generalizing the proven rr-digest pattern
(`SADb_audit\rr_digest_2023plus.md`, `rr_digest_20260912.py`). One run per month
drafts a shortlist of NEW papers (OpenAlex + Crossref only — both open, no key;
Google Scholar / Web of Science scraping is deliberately NOT done, standing
constraint). **The machine drafts, the human rules: nothing is ever written to
Airtable or Zotero automatically.**

## Files

- `monthly_watch.py` — the sweep (myo env, stdlib only).
- `watch_config.json` — Ben's watch list. `authors` = union of
  `authors_derived` (every corpus primary author with ≥4 *retained* papers,
  computed 2026-10-22 from `export\sadb_export.json`, papers only — never the
  non-paper records) and `authors_pinned` (Rybak, Shevtsova, McCrea, Deng,
  Grillner, Pearson, Prochazka, Jankowska), deduped — 20 names. `keywords`
  seeded from the corpus vocabulary (Ia/Ib/group II afferents, Renshaw, spinal
  locomotion, central pattern generator, locomotor CPG, proprioceptive
  feedback, fictive locomotion, muscle spindle, Golgi tendon organ, cutaneous
  reflex, locomotor network, hindlimb locomotion, locomotor recovery).
  Deliberately NOT seeded even though the corpus mentions them: "sensory
  feedback" and "spinal cord injury" — the first real run showed those two
  alone produce ~40 of 53 rows, nearly all clinical/haptics papers outside the
  corpus's lane; edit the config to add them back if wanted.
  `mailto` is for the OpenAlex/Crossref polite pool — LEFT EMPTY on purpose;
  never hardcode a personal email in this public repo. `domain_gate` scopes
  author hits (surname matches alone pull in unrelated name-holders, e.g. the
  wireless "Jean-François Frigon").
- `state.json` — created by runs: routed DOIs + each month's stored report
  rows. This is what makes re-runs idempotent.
- `reports\YYYY-MM.md` + `reports\YYYY-MM_shortlist.json` — the digest table
  and the batch3-pattern import file (`title` / `authors` "Family Given; …" /
  `container` / `year` / `doi` / `first_author` / `trigger` / `why` / `fetch` —
  the same row shape as `batch3\rr_digest_shortlist.json`, so
  `SADb_audit\zotero_add_digest.py` and the Airtable Name/Author/Year/DOI
  import flow keep working unchanged).

## Run

```
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" D:\Github\Bipedal_Robot\SADb_audit\watch\monthly_watch.py
```

Defaults (all match the spec): window = first day of the PREVIOUS month
through today (one-month lookback plus small overlap); report month = current
month; output root = this folder (`reports\` + `state.json` here). A full run
is ~70 modest API calls with a 0.5 s delay — a couple of minutes.

CLI overrides for testing:

- `--config PATH` — alternate config (e.g. a scratch copy under `tmp\`).
- `--out PATH` — output ROOT override; `reports\` and `state.json` go under
  it, so scratch runs never touch the real state.
- `--month YYYY-MM` — report month; window = previous month's first day
  through that month's last day (or today if current).
- `--window YYYY-MM-DD YYYY-MM-DD` — explicit window; report month = END's
  month.

## Idempotency

`state.json` tracks every DOI already routed (corpus members excluded
upstream) and stores each month's report rows. Re-running a month rebuilds the
report from stored rows + genuinely new hits — rows never duplicate, and a
re-run that finds nothing new rewrites the identical bytes. Each month keeps
the window of its first run (a later re-run doesn't stretch it) unless
`--window` is passed explicitly. An empty shortlist is a valid result and says
so in the report.

## Pipeline (exactly the established rr-digest one)

1. Pull candidates: OpenAlex (one domain-gated query per watch author —
   `raw_author_name.search:<surname>` AND a `domain_gate` term in
   `title_and_abstract.search`, because surname-only matching returns tens of
   thousands of works for common names — plus grouped pipe-OR keyword queries;
   all inside the date window) and Crossref (one query per author / keyword,
   pub-date window). Small pages, 0.5 s delay, one backoff retry on 429/5xx.
2. Dedupe against the corpus DOI set (all 943 papers — retained AND archived —
   from `export\sadb_export.json`), against every previously routed DOI, and
   against version twins of already-routed works (eLife publishes versioned
   `.2` DOIs whose title text drifts; without this a re-run would re-add the
   same paper under the version DOI).
3. Collapse preprint/version twins: same normalized title (keep the published
   version, prefer the one with a journal container), then same
   version-stripped DOI with a compatible title (same guard both directions).
4. Trigger: author hit (family name matches a watch author AND a
   `domain_gate` term appears in title/abstract) or keyword hit (exact
   case-insensitive substring). The "why relevant" line names what fired.
5. Crossref-verify EVERY shortlist entry (`GET /works/{DOI}`); failures are
   dropped and logged, never routed.
6. Write `reports\YYYY-MM.md` + `reports\YYYY-MM_shortlist.json`, then commit
   state. Run log (counts per run) is kept in `state.json` only.

## OpenAlex budget reality (no-key tier)

OpenAlex's anonymous tier (no API key, per the standing no-credentials rule)
shares a small per-IP daily budget — observed 2026-10-23: 1000 credits/day
(10 per query, $0.10), **reset at midnight UTC**. When it is exhausted every
OpenAlex call answers HTTP 429 with a `Retry-After` of hours. The script
handles this honestly: after retries fail it marks OpenAlex unavailable for
the rest of that run, skips the remaining OpenAlex queries (one run-note,
also printed to the report's *Run notes* section), and Crossref carries the
sweep. A later run with budget simply folds in whatever the starved run
missed — the routed-DOI dedupe makes that merge safe. The first full run of
2026-10-22/23 exhausted the pool with diagnosing + test sweeps, so the
2026-10 report of record was swept Crossref-complete but OpenAlex-starved
(disclosed in its Run notes); re-run after midnight UTC for the full sweep.
Crossref has no such budget.

## Scheduling (provided here, NOT installed)

Monthly on the 1st at 09:00. Paste into **cmd** — the `\"` escaping is
required so the inner quotes survive cmd's argument parsing (the python path
has spaces):

```
schtasks /Create /SC MONTHLY /D 1 /ST 09:00 /TN SADbMonthlyWatch /TR "\"C:\Users\Ben Bolen\.conda\envs\myo\python.exe\" D:\Github\Bipedal_Robot\SADb_audit\watch\monthly_watch.py"
```

Check / delete with `schtasks /Query /TN SADbMonthlyWatch` /
`schtasks /Delete /TN SADbMonthlyWatch /F`. A ZCode-side cron is an acceptable
alternative, but the script is standalone (no imports from the repo, stdlib
only, myo python by full path).
