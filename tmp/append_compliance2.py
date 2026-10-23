base = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\COMPLIANCE_TABLE.md"
add = open(r"D:\Github\Bipedal_Robot\tmp\compliance_append2.md", encoding="utf-8").read()
with open(base, "a", encoding="utf-8") as f:
    f.write(add)
print("appended", len(add), "chars")
