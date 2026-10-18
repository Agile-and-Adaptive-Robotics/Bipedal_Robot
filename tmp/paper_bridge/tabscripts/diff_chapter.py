# -*- coding: utf-8 -*-
"""Numeric diff: fragment tables T1/T2 vs the chapter tables of record
(30-torque.tex tab:FlxPinGof / tab:ExtPinGof, committed state).
Compares the six numeric cells of every row."""
import re

CH = (r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
      r"\ProofFinal\chapters\30-torque.tex")
FR = r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge\sections\tables_main.tex"


def rows_from(path, begin_marker, end_marker):
    src = open(path, encoding="utf-8").read()
    block = src[src.find(begin_marker):src.find(end_marker, src.find(begin_marker))]
    out = []
    for line in block.splitlines():
        line = line.strip()
        if "&" not in line or not line.endswith("\\\\"):
            continue
        cells = line[:-2].split("&")
        nums = []
        for cell in cells[1:]:  # drop label cell
            nums.extend(float(x) for x in re.findall(r"-?\d+\.\d+", cell))
        if nums:
            out.append((cells[0].strip()[:34], nums))
    return out


def strip_qty(label):
    return re.sub(r"\$l_\{0\}=\\qty\{|\}\{\\cm\}\$|\$|\\,", "", label)


pairs = [
    ("T1 tab:FlxPinGof",
     (CH, r"\label{tab:FlxPinGof}", r"\end{table}"),
     (FR, r"\label{tab:FlxPinGof}", r"\end{table}")),
    ("T2 tab:ExtPinGof",
     (CH, r"\label{tab:ExtPinGof}", r"\end{table}"),
     (FR, r"\label{tab:ExtPinGof}", r"\end{table}")),
]

fails = 0
for name, chap, frag in pairs:
    c_rows = rows_from(*chap)
    f_rows = rows_from(*frag)
    print("%s: chapter %d rows, fragment %d rows" % (name, len(c_rows), len(f_rows)))
    if len(c_rows) != len(f_rows):
        fails += 1
        continue
    for (cl, cn), (fl, fn) in zip(c_rows, f_rows):
        ok = cn == fn and strip_qty(cl) == strip_qty(fl)
        print("  %-36s %s  %s" % (strip_qty(fl),
                                  "labels-match" if strip_qty(cl) == strip_qty(fl)
                                  else "LABEL-DIFF [%s vs %s]" % (strip_qty(cl), strip_qty(fl)),
                                  "numbers-match" if cn == fn else
                                  "NUMBER-DIFF %s vs %s" % (cn, fn)))
        if not ok:
            fails += 1
print("\n" + ("DIFF CLEAN" if fails == 0 else "%d ROW DIFFERENCES" % fails))
