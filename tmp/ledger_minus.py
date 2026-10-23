note = """
- 2026-10-23 (ZCode, round-4 review of Revision C): one self-scored **minus score**
  added — round 2 presented the Excel "Torque actual" rows (28.8–50.6 N·m span,
  259.9 mm) as final adjoint-reduced presentation data; Ben clarified rows 16–17
  of the workbook are sanity checks, and Revision C removed the paragraph.
  Lesson: a number matching its source file is not provenance — verify the row's
  intended status before presenting it as a result.
"""
with open(
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\summary.md",
    "a",
    encoding="utf-8",
) as f:
    f.write(note)
print("ledger updated")
