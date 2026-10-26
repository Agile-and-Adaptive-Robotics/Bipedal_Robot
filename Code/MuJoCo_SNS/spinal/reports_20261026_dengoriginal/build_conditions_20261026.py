r"""Build Deng-2019-original perturbation conditions from the exported asim.

Order B1 (Ben, 2026-10-26). Source of truth: deng2019_original_export.asim
(the AnimatLab GUI's own serializer run headless via DengExport.exe on the
ModelDB 264844 project; 90 neurons / 196 connexions / 42 adapters verified
GUID-identical to the mine json).

Conditions (PREREG_B1_deng_original_20261026.md, written BEFORE any run):
  B0_baseline    - stimuli exactly as delivered (RG_Stin/PF_Stin disabled,
                   Stimulus_1 torque-z 0.05 N.m on, 2.0-2.4 s)
  P1_rgflexor    - RG_Stin retargeted to LH_RG F (e8db7e8f), +2 nA, 2.0-2.5 s
  P2_kneepfflexor- PF_Stin retargeted to LH_K&A PF F (48ad42ce), +2 nA, 2.0-2.5 s
  P3_authorsblock- RG_Stin enabled AS DELIVERED (-10 nA on LH_RG E, 2.0-2.1 s)

Common edits on every condition copy:
  SetSimEnd True, SimEndTime 7.01 -> 10.1 s (window at 2.0-2.5 + ~4 cycles)
  existing charts: StartTime 0, EndTime 10, CollectInterval 1 ms
  + one new FileChart "deng_rg_mn" (1 ms) recording the prereg trace set:
    RG E/F both sides, knee MN E/F both sides, LH K&A PF E/F, LH_HipZ PF E

Every edit is logged to <run>/edits.json. Deterministic new GUIDs b1d1%04d.
Usage: D:\Anaconda\envs\myo\python.exe build_conditions_20261026.py
"""
from __future__ import annotations
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "deng2019_original_export.asim")

SIM_END = 10.1
CHART_END = 10
DT_CHART = 0.001

# key GUIDs (from deng2019_mine_20261025.json, cross-checked in the export)
LH_RG_E  = "a888264c-bb2d-4bd6-ac77-6a4b8aa01656"
LH_RG_F  = "e8db7e8f-3dc7-43af-b1fa-25875605031b" if False else "e8db7e8f"  # full below
LH_RG_F  = "e8db7e8f-8d99-4a29-9b2d-6bf0f95cf630" if False else None
# NOTE: fill full GUIDs from the mine at runtime instead of hardcoding halves.
MINE = os.path.join(HERE, "originals", "deng2019_mine_20261025.json")

RG_STIN_ID = "51a2ee1b-6b96-44a0-8567-b92cdad012b4"
PF_STIN_ID = "38bda618-a888-4293-b0d2-557ead8fc3db"

def g(n: int) -> str:
    return "b1d1%04d-0000-4000-8000-%012d" % (n, 900000000000 + n)

def load_nodes() -> dict:
    m = json.load(open(MINE, encoding="utf-8"))
    out = {}
    for n in m["nodes"]:
        if n["class"] != "NonSpiking" or not n["name"]:
            continue   # OffPage connectors and non-neurons share names - exclude
        if n["name"] in out and out[n["name"]] != n["guid"]:
            raise RuntimeError("duplicate NonSpiking name: " + n["name"])
        out[n["name"]] = n["guid"]
    return out

def patch_stimulus(x: str, stim_id: str, edits: dict, what: str, log: list) -> str:
    i = x.find("<ID>%s</ID>" % stim_id)
    if i < 0:
        raise RuntimeError("stimulus %s not found" % stim_id)
    s = x.rfind("<Stimulus>", 0, i)
    e = x.find("</Stimulus>", i) + len("</Stimulus>")
    blk = x[s:e]
    for tag, val in edits.items():
        pat = r"<%s>[^<]*</%s>" % (tag, tag)
        new = "<%s>%s</%s>" % (tag, val, tag)
        blk2, n = re.subn(pat, new, blk, count=1)
        if n != 1:
            raise RuntimeError("%s: tag %s patch count %d" % (what, tag, n))
        blk = blk2
        log.append({"stim": what, "tag": tag, "value": str(val)})
    return x[:s] + blk + x[e:]

def common_edits(x: str, log: list) -> str:
    x2, n = re.subn(r"<SetSimEnd>False</SetSimEnd>", "<SetSimEnd>True</SetSimEnd>", x, count=1)
    if n != 1: raise RuntimeError("SetSimEnd patch %d" % n)
    x = x2; log.append({"edit": "SetSimEnd", "to": "True"})
    x2, n = re.subn(r"<SimEndTime>7\.01</SimEndTime>", "<SimEndTime>%s</SimEndTime>" % SIM_END, x, count=1)
    if n != 1: raise RuntimeError("SimEndTime patch %d" % n)
    x = x2; log.append({"edit": "SimEndTime", "to": SIM_END})
    x2, n = re.subn(r"<StartTime>0\.5</StartTime>", "<StartTime>0</StartTime>", x)
    log.append({"edit": "chart StartTime 0.5->0", "count": n}); x = x2
    x2, n = re.subn(r"<StartTime>4</StartTime>", "<StartTime>0</StartTime>", x)
    log.append({"edit": "chart StartTime 4->0", "count": n}); x = x2
    x2, n = re.subn(r"<EndTime>5</EndTime>", "<EndTime>%d</EndTime>" % CHART_END, x)
    log.append({"edit": "chart EndTime 5->10", "count": n}); x = x2
    x2, n = re.subn(r"<EndTime>9</EndTime>", "<EndTime>%d</EndTime>" % CHART_END, x)
    log.append({"edit": "chart EndTime 9->10", "count": n}); x = x2
    x2, n = re.subn(r"<EndTime>7</EndTime>", "<EndTime>%d</EndTime>" % CHART_END, x)
    log.append({"edit": "chart EndTime 7->10", "count": n}); x = x2
    x2, n = re.subn(r"<CollectInterval>5e-005</CollectInterval>",
                    "<CollectInterval>%s</CollectInterval>" % DT_CHART, x)
    log.append({"edit": "CollectInterval 5e-5->1ms", "count": n}); x = x2
    return x

def add_neural_chart(x: str, nodes: dict, log: list) -> str:
    cols = [
        ("LH_RG_E",   nodes["LH_RG E"]),
        ("LH_RG_F",   nodes["LH_RG F"]),
        ("RH_RG_E",   nodes["RH_RG E"]),
        ("RH_RG_F",   nodes["RH_RG F"]),
        ("LH_KneeMN_E", nodes["LH_Knee MN E"]),
        ("LH_KneeMN_F", nodes["LH_Knee MN F"]),
        ("RH_KneeMN_E", nodes["RH_Knee MN E"]),
        ("RH_KneeMN_F", nodes["RH_Knee MN F"]),
        ("LH_KAPF_E", nodes["LH_K&A PF E"]),
        ("LH_KAPF_F", nodes["LH_K&A PF F"]),
        ("LH_HipPF_E", nodes["LH_HipZ PF E"]),
    ]
    chart = ['<DataChart>', '<ID>%s</ID>' % g(1), '<Type>FileChart</Type>',
             '<OutputFilename>deng_rg_mn.txt</OutputFilename>', '<Name>deng_rg_mn</Name>',
             '<SetStartEndTime>True</SetStartEndTime>', '<CollectTimeWindow>-1</CollectTimeWindow>',
             '<AlwaysActive>False</AlwaysActive>', '<StartTime>0</StartTime>',
             '<EndTime>%d</EndTime>' % CHART_END, '<CollectInterval>%s</CollectInterval>' % DT_CHART,
             '<DataColumns>']
    for k, (cn, tid) in enumerate(cols):
        chart += ['<DataColumn>', '<ModuleName/>', '<Type>DataColumn</Type>',
                  '<ID>%s</ID>' % g(10 + k), '<ColumnName>%s</ColumnName>' % cn,
                  '<DataType>MembraneVoltage</DataType>', '<TargetID>%s</TargetID>' % tid,
                  '</DataColumn>']
        log.append({"chart_col": cn, "target": tid})
    chart += ['</DataColumns>', '</DataChart>']
    anchor = x.rfind("</DataChart>") + len("</DataChart>")
    return x[:anchor] + "\n" + "\n".join(chart) + x[anchor:]

CONDITIONS = {
    "B0_baseline": lambda x, nodes, log: x,
    "P1_rgflexor": lambda x, nodes, log: patch_stimulus(
        x, RG_STIN_ID,
        {"Enabled": "True", "TargetNodeID": nodes["LH_RG F"],
         "CurrentOn": "2e-009", "EndTime": "2.5"},
        "P1 RG_Stin->LH_RG F +2nA 2.0-2.5", log),
    "P2_kneepfflexor": lambda x, nodes, log: patch_stimulus(
        x, PF_STIN_ID,
        {"Enabled": "True", "TargetNodeID": nodes["LH_K&A PF F"],
         "CurrentOn": "2e-009", "EndTime": "2.5"},
        "P2 PF_Stin->LH_K&A PF F +2nA 2.0-2.5", log),
    "P3_authorsblock": lambda x, nodes, log: patch_stimulus(
        x, RG_STIN_ID, {"Enabled": "True"},
        "P3 RG_Stin as delivered (-10nA LH_RG E 2.0-2.1)", log),
    # paper-exact windows (Deng 2023: 2 nA, 2.0-2.1 s - the 2.5 in Ben's order
    # was an approximation; the paper and the project's own blocks say 2.1)
    "P1b_rgflexor_2p1": lambda x, nodes, log: patch_stimulus(
        x, RG_STIN_ID,
        {"Enabled": "True", "TargetNodeID": nodes["LH_RG F"],
         "CurrentOn": "2e-009", "EndTime": "2.1"},
        "P1b RG_Stin->LH_RG F +2nA 2.0-2.1 (paper-exact)", log),
    "P2b_kneepfflexor_2p1": lambda x, nodes, log: patch_stimulus(
        x, PF_STIN_ID,
        {"Enabled": "True", "TargetNodeID": nodes["LH_K&A PF F"],
         "CurrentOn": "2e-009", "EndTime": "2.1"},
        "P2b PF_Stin->LH_K&A PF F +2nA 2.0-2.1 (paper-exact)", log),
}

def main() -> int:
    nodes = load_nodes()
    # sanity: the two stin targets must match the mine
    assert nodes["LH_RG E"].startswith("a888264c"), nodes["LH_RG E"]
    assert nodes["LH_K&A PF F"].startswith("48ad42ce"), nodes["LH_K&A PF F"]
    src = open(SRC, encoding="utf-8", errors="surrogateescape").read()
    for cond, fn in CONDITIONS.items():
        log = []
        x = common_edits(src, log)
        x = fn(x, nodes, log)
        x = add_neural_chart(x, nodes, log)
        # validate XML well-formedness of the whole file
        import xml.etree.ElementTree as ET
        ET.fromstring(x.encode("utf-8", "surrogateescape"))
        rd = os.path.join(HERE, "run_" + cond)
        os.makedirs(rd, exist_ok=True)
        with open(os.path.join(rd, cond + ".asim"), "w", encoding="utf-8",
                  errors="surrogateescape", newline="") as f:
            f.write(x)
        json.dump({"condition": cond, "edits": log},
                  open(os.path.join(rd, "edits.json"), "w", encoding="utf-8"), indent=2)
        print("built run_%s/%s.asim (%d bytes, %d edits)" % (cond, cond, len(x), len(log)))
    return 0

if __name__ == "__main__":
    sys.exit(main())
