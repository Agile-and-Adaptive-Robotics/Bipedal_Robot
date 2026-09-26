# M4 inspector 2: c1/V3 node types + bilateralrg RG-included node list
import io, json, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
d = json.load(open(SP + r"\connectome_templates.json", encoding="utf-8"))
t = d["bilateralrg"]
print("=== bilateralrg: all node labels/types (RG-side families) ===")
for n in t["nodes"]:
    if n["type"] in ("IN-C", "IN-V3", "HC-RG-E", "HC-RG-F", "IN-InE", "IN-InF", "PORT-load"):
        print(f"  {n['type']:<10} {n['label']}")
