import os
p = r"D:\Github\Bipedal_Robot\SADb_audit\README.md"
t = open(p, encoding="utf-8").read()
anchor = "(one row per paper + WORKFLOW rows; field meanings in its header row)"
add = """

### How to read curation_log.csv (for Ben)

Each row is one paper. Left to right: date, batch number (1-5 = curation batches; RR = Research-Rabbit discovery; VV = VOSviewer; PDFS/AUTHORS/RULINGS = maintenance passes), the Zotero key, the Airtable record id, a short title, then yes/no "note written", the Animal tags, the Feedback links, and a free-text column saying where the paper's text came from and anything Ben must rule on. Rows whose zotero_key is WORKFLOW are not papers - they are milestone/verification notes and your rulings. For the current backlog state, read only the WORKFLOW rows."""
if anchor in t and "How to read curation_log.csv" not in t:
    t = t.replace(anchor, anchor + add, 1)
    open(p, "w", encoding="utf-8").write(t)
    print("inserted")
else:
    print("anchor missing or already present")
