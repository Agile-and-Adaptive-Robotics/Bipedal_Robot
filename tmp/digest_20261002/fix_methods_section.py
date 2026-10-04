"""Repair the methods margin/section region (v2): keep one audit comment, restore the
balance-section header as a standalone line, drop the commented duplicate fragment."""
from pathlib import Path

p = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\20-methods.tex")
t = p.read_text(encoding="utf8")

marker = "% [ZCODE 2026-09-26 Xi-audit: margin sentence updated."
i0 = t.find(marker)
assert i0 > 0
i1 = t.find(marker, i0 + 1)
assert i1 > 0, "second audit copy not found"
eol = t.find("\n", i1)
tail = t[i1:eol]
assert "\\par\\section{Balance Platform and Inverted-Pendulum Testbeds}" in tail, "section glue missing"
assert tail.count("With two BPAs per joint") == 1

audit1 = t[i0:i1].rstrip("\n")          # first audit comment, complete line
sec = "\\section{Balance Platform and Inverted-Pendulum Testbeds}\\label{sec:balance_testbeds}"
t = t[:i0] + audit1 + "\n\n" + sec + "\n" + t[eol + 1:]

assert t.count(marker) == 1, t.count(marker)
assert t.count(sec) == 1
assert t.count("\\label{sec:balance_testbeds}") == 1
p.write_text(t, encoding="utf8", newline="\n")
print("repaired: audit copies =", t.count(marker), "| sections =", t.count(sec))
