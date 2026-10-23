"""Final read-back of export/sadb_export.json (read-only)."""
import json, os
p = os.path.join(r"D:\Github\Bipedal_Robot\SADb_audit", "export", "sadb_export.json")
recs = json.load(open(p, encoding="utf-8"))
by_type, by_type_arch = {}, {}
for r in recs:
    rt = r.get("recordType", "paper")
    by_type[rt] = by_type.get(rt, 0) + 1
    if r.get("archived"):
        by_type_arch[rt] = by_type_arch.get(rt, 0) + 1
arch = sum(1 for r in recs if r.get("archived"))
print("total:", len(recs))
print("by type:", by_type)
print("archived by type:", by_type_arch)
print("archived total:", arch, "| live:", len(recs) - arch)
papers = [r for r in recs if r.get("recordType", "paper") == "paper"]
print("papers:", len(papers), "(archived", sum(1 for r in papers if r.get("archived")), ")")
nonp = [r for r in recs if r.get("recordType", "paper") != "paper"]
print("non-paper:", len(nonp), "(archived", sum(1 for r in nonp if r.get("archived")), ")")
# spot-check orphan reconstruction sources of record
byid = {r["id"]: r for r in recs}
for rid, want in [("recX0lfqPjQtMSR9q", ("model", "Shevtsova 2016")),
                  ("recnvMUdPsOHvevrx", ("model", "Geyer and Herr 2010")),
                  ("recze0Jzfzci1xDVA", ("model", "Rybak 2006")),
                  ("rec4p9PNw5BZoSDSP", ("review", "Prochazka 2007")),
                  ("rec2uSKIvzoT6tszs", ("feedback", "Ipsilateral excitatory - swimming edge cell")),
                  ("recIXfJZUYWQ0HNyJ", ("feedback", "Cross inhibitory - swimming edge cell")),
                  ("recmmxPL2uZaMON2N", ("feedback", ""))]:
    r = byid.get(rid)
    assert r and r["recordType"] == want[0] and r["title"] == want[1] and r.get("archived"), \
        f"orphan spot-check failed for {rid}: {r}"
print("orphan spot-check: 7/7 manifest satellite records reconstructed with names from prune_shortlist.csv")
one = next(r for r in recs if r["id"] == next(iter(__import__("json").load(
    open(os.path.join(r"D:\Github\Bipedal_Robot\SADb_audit", "review_orphans_20261002.json"),
         encoding="utf-8")))))
print("review-orphan sample:", one["id"], "|", one["recordType"], "|", one["title"], "| archived:", one["archived"])
