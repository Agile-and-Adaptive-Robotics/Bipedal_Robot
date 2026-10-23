import collections
import re


def widths(path):
    log = open(path, encoding="utf-8", errors="replace").read()
    i = log.rfind("=== PASS 4 ===")
    tail = log[i:]
    return collections.Counter(re.findall(r"Overfull \\hbox \(([\d.]+)pt", tail))


old = widths(r"D:\Github\Bipedal_Robot\tmp\build_live.log")
new = widths(r"D:\Github\Bipedal_Robot\tmp\build_live3.log")
print("old total:", sum(old.values()), "new total:", sum(new.values()))
print("removed:", old - new)
print("added:", new - old)
