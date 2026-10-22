import os

os.chdir(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022")
for f in ("REPORT.md", os.path.join("patches", "README.md")):
    s = open(f, encoding="utf-8").read()
    n = s.count("Notes/advisor_revision_20261022")
    s = s.replace(
        "Notes/advisor_revision_20261022",
        "Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022",
    )
    # avoid doubled path if some occurrences were already full
    s = s.replace(
        "Documentation/Reports and Papers/Dissertation/Documentation/Reports and Papers/Dissertation/",
        "Documentation/Reports and Papers/Dissertation/",
    )
    with open(f, "w", encoding="utf-8", newline="") as fh:
        fh.write(s)
    print(f, "replaced", n)
