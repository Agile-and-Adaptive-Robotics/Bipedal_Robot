p = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\APP_1393_TASK.md"
s = open(p, encoding="utf-8").read()
old = "# Side task: carry the full 1,393-record pre-prune SADb base in the packaged app"
new = "# Side task (2 deliverables): (1) full 1,393-record pre-prune base in the packaged app; (2) the monthly skimmer"
if old in s:
    s = s.replace(old, new)
    with open(p, "w", encoding="utf-8", newline="") as f:
        f.write(s)
    print("title updated")
else:
    print("title already:", s.splitlines()[0])
