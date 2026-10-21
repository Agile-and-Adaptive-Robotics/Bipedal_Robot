r"""Stock-connectome 30 s reference: Ben's OLD Ground_Standalone.asim (the
battle-tested 2026-10-16 V0 control source), UNTOUCHED wiring, only the sim
window extended (SimEndTime 5.1 -> 30.1, chart windows [0,30], 1 ms collect)
plus the same Root Track chart. Purpose: does the OLD connectome free-drop
stand/march for 30 s, or did the old 5.1 s runs simply end before the fall?
Writes the asim ONLY into this campaign folder (run_stock30\).
"""
from __future__ import annotations

import os
import re
import sys

REPO = r"D:\GitHub\Bipedal_Robot"
SRC = os.path.join(REPO, "Neuromechanical_Models", "Walker_2_Layer_CPG_BilateralRG",
                   "Walker_2_Layer_CPG_BilateralRG_Ground_Standalone.asim")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "run_stock30")
OUT = os.path.join(OUT_DIR, "Walker_2_Layer_CPG_BilateralRG_stock30.asim")

SIM_END, CHART_END = 30.1, 30
NEWPFX = "d73a"


def g(n: int) -> str:
    return f"{NEWPFX}{n:04d}-0000-4000-8000-{900000000000 + n}"


def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(SRC, encoding="utf-8", errors="surrogateescape", newline="") as f:
        s = f.read()
    if NEWPFX in s:
        raise RuntimeError("prefix present")

    s, n = re.subn(r"<SimEndTime>5.1</SimEndTime>", f"<SimEndTime>{SIM_END}</SimEndTime>", s)
    assert n == 1
    s, n = re.subn(r"<EndTime>5</EndTime>", f"<EndTime>{CHART_END}</EndTime>", s)
    s, n = re.subn(r"<StartTime>0.5</StartTime>", "<StartTime>0</StartTime>", s)
    s, n = re.subn(r"<CollectInterval>0.0002</CollectInterval>",
                   "<CollectInterval>0.001</CollectInterval>", s)

    i_root = s.find("<Name>Root</Name>")
    root_id = re.search(r"<ID>([^<]+)</ID>", s[i_root:i_root + 200]).group(1)
    joints = {}
    for m in re.finditer(r"<Joint>(.*?)</Joint>", s, re.S):
        b = m.group(1)
        nm = re.search(r"<Name>([^<]*)</Name>", b)
        gid = re.search(r"<ID>([^<]+)</ID>", b)
        if nm and gid:
            joints[nm.group(1)] = gid.group(1)
    cols = [("RootX", "WorldPositionX", root_id),
            ("RootY", "WorldPositionY", root_id),
            ("RootZ", "WorldPositionZ", root_id),
            ("hip_L", "JointRotation", joints["hip_L"]),
            ("hip_R", "JointRotation", joints["hip_R"])]
    chart = [f"<DataChart>\n<ID>{g(41)}</ID>\n<Type>FileChart</Type>\n"
             f"<OutputFilename>Root Track.txt</OutputFilename>\n"
             f"<Name>Root Track</Name>\n<SetStartEndTime>True</SetStartEndTime>\n"
             f"<CollectTimeWindow>-1</CollectTimeWindow>\n"
             f"<AlwaysActive>False</AlwaysActive>\n<StartTime>0</StartTime>\n"
             f"<EndTime>{CHART_END}</EndTime>\n<CollectInterval>0.001</CollectInterval>\n"
             f"<DataColumns>"]
    for n2, (cn, dt, tid) in enumerate(cols):
        chart.append(f"<DataColumn>\n<ModuleName/>\n<Type>DataColumn</Type>\n"
                     f"<ID>{g(42 + n2)}</ID>\n<ColumnName>{cn}</ColumnName>\n"
                     f"<DataType>{dt}</DataType>\n<TargetID>{tid}</TargetID>\n"
                     f"</DataColumn>")
    chart.append("</DataColumns>\n</DataChart>")
    a2 = s.rfind("</DataChart>") + len("</DataChart>")
    s = s[:a2] + "\n" + "\n".join(chart) + s[a2:]

    import xml.etree.ElementTree as ET
    ET.fromstring(s.encode("utf-8", "surrogateescape"))
    with open(OUT, "w", encoding="utf-8", errors="surrogateescape", newline="") as f:
        f.write(s)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
