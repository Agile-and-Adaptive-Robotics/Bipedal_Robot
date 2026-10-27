import json

m = json.load(open("Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json"))
try:
    b = json.load(open("Solid_Models/OpenSim/Gait2392_Robotbody/bpa_actuators_27.json"))
    dias = {}
    for a in b.get("actuators", b if isinstance(b, list) else []):
        name = a.get("name", a.get("master_name", ""))
        d = a.get("diameter_mm", a.get("diameter", a.get("bpa_diameter_mm", "?")))
        dias[name] = d
    src = "bpa_actuators_27.json"
except Exception as e:
    dias, src = {}, f"bpa_actuators_27.json unreadable ({e})"

print(f"| # | Muscle | Primary DoF | Secondary DoF(s) | BPA dia (mm) |")
print(f"|---|---|---|---|---|")
for a in m["actuators"]:
    dia = dias.get(a["master_name"], dias.get(str(a["id"]), "?"))
    sec = ", ".join(a["secondary_dofs"]).replace("_r", "") or "--"
    print(f"| {a['id']} | {a['master_name'].replace(',', ',')} | {a['primary_dof'].replace('_r','')} | {sec} | {dia} |")
print(f"\n(diameter source: {src}; ? = not in that file — fall back to PAM_Prop.mat)")
