import os

os.chdir(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\proposals")

subs = {
    "15-background.tex": [
        ("Work from the author's own program", "Work from my own program"),
        ("In the author's master's work", "In my master's work"),
        ("From the author's master's study.", "From my master's study."),
        ("the author's own included", "my own included"),
        ("designed at the author's laboratory", "designed at my laboratory"),
        ("the author's program first applied", "my program first applied"),
    ],
    "60-sns.tex": [
        ("the author and Morrow constructed", "Morrow and I constructed"),
        ("The toe joint, the author's addition to the Li design", "The toe joint, my addition to the Li design"),
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

# verify none remain in any proposal
import glob

for fn in glob.glob("*.tex"):
    s = open(fn, encoding="utf-8").read()
    if "the author" in s:
        i = s.find("the author")
        print("REMAINS in", fn, ":", s[i - 60 : i + 60])
print("sweep done")
