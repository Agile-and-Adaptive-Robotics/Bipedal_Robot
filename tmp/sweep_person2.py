import os

os.chdir(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\proposals")

subs = {
    "70-overall.tex": [
        ("encountered in the author's master's routing work", "encountered in my master's routing work"),
    ],
    "80-conclusion.tex": [
        ("the optimized muscle routing of the author's earlier work", "the optimized muscle routing of my earlier work"),
    ],
    "94-AppendixC.tex": [
        ("set from CAD by the author", "set from CAD by me"),
        ("the author's fcec-ray construction", "my fcec-ray construction"),
        ("its axis assignment was revised by the author", "its axis assignment was revised by me"),
    ],
}

for fn, pairs in subs.items():
    s = open(fn, encoding="utf-8").read()
    for a, b in pairs:
        n = s.count(a)
        s = s.replace(a, b)
        print(fn, n, "x", a[:45])
    with open(fn, "w", encoding="utf-8", newline="") as f:
        f.write(s)

for fn in ("70-overall.tex", "80-conclusion.tex", "94-AppendixC.tex"):
    s = open(fn, encoding="utf-8").read()
    print(fn, "remaining 'the author':", s.count("the author"))
