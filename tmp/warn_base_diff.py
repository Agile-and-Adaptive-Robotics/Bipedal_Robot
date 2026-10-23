import collections
import re


def widths(path):
    log = open(path, encoding="utf-8", errors="replace").read()
    i = log.rfind("=== PASS 4 ===")
    tail = log[i:]
    return collections.Counter(re.findall(r"Overfull \\hbox \(([\d.]+)pt", tail))


base = widths(r"D:\Github\Bipedal_Robot\tmp\build_base.log")
new = widths(r"D:\Github\Bipedal_Robot\tmp\build_live3.log")
print("base total:", sum(base.values()), "new total:", sum(new.values()))
print("in base not in new:", base - new)
print("in new not in base:", new - base)
