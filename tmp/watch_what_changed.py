import json

rows = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\watch\reports\2026-10_shortlist.json", encoding="utf-8"))
st = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\watch\state.json", encoding="utf-8"))
log = st["run_log"]
print("run_log:")
for e in log:
    print("  ", e["date"], "| new_entries:", e["new_entries"], "| shortlist_total:", e["shortlist_total"],
          "| raw:", e["raw_hits"], "| trigger_hits:", e["trigger_hits"], "| twins:", e["twins_dropped"])

# The added row = the one whose routing happened in run 2; find rows present in
# shortlist and diff against what run 1 stored... state now holds all 17, so
# instead show which row is NOT among the first 16 of the previous report order.
prev = set(json.load(open(r"D:\Github\Bipedal_Robot\tmp\watch_prev16.json", encoding="utf-8"))) if False else None
# Simpler: report rows whose DOI routing month is 2026-10 (all) — print the row
# added by run 2 by comparing to the earlier saved copy if we have one; we do not,
# so list the last row of the current sort (the most likely insert at the tail).
print("\ncurrent rows (order):")
for r in rows:
    print("  %-34s | %s | %s | %s" % (r["doi"], r["year"], r["first_author"], r["why"][:60]))
