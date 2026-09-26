# M4 inspector: bilateralrg template commissural edges + w2laproj pf_drive edges
import io, json, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
d = json.load(open(SP + r"\connectome_templates.json", encoding="utf-8"))

WANT = {"comm_c1", "c1_SynAmp2.749", "comm_v3", "v3_SynAmp0.1_weak",
        "v3_to_contra_InE", "kickoff_antiphase", "pf_drive", "rg_laminate"}
t = d["bilateralrg"]
print("=== bilateralrg: edges of interest ===")
for e in t["edges"]:
    if e["tag"] in WANT:
        print(json.dumps(e, ensure_ascii=False))
print()
print("=== bilateralrg: RG/comm/C1/V3 nodes ===")
for n in t["nodes"]:
    lab = n["label"]
    if any(k in lab for k in ("RG", "comm", "Comm", "C1", "V3", "Stim")):
        print(json.dumps(n, ensure_ascii=False))

print()
print("=== w2laproj: pf_drive + rg_laminate edges ===")
t2 = d["w2laproj"]
for e in t2["edges"]:
    if e["tag"] in ("pf_drive", "rg_laminate"):
        print(json.dumps(e, ensure_ascii=False))
print()
print("=== w2laproj nodes containing 'PF' or 'RG' or Stim ===")
for n in t2["nodes"]:
    if any(k in n["label"] for k in ("RG", "PF", "Stim")):
        print(json.dumps(n, ensure_ascii=False))
