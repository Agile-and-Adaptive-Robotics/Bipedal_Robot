note = """

## NOTE (2026-10-22 round 2): measurements above are ROUND-1 numbers

Round 2 applied the principal's selection pass after this table was written.
Current numbers of record are in REPORT.md (SNS 6,794 / SADb 2,063 narrative;
280 pages). Round-2 additions since the round-2 tables: the SADb prune paragraph
now states the verified pre-prune base size of 1,393 records (957 post-prune +
436 removed; papers 561 + 382 = 943 audited; consistent with deletion_manifest.json
88+290+4 papers and 3+1+3+47 orphaned links) per Ben's recollection ("like 1300").
"""
with open(
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\COMPLIANCE_TABLE.md",
    "a",
    encoding="utf-8",
) as f:
    f.write(note)
print("appended")
