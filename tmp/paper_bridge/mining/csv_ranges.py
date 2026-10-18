import csv

BASE = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
for f in ["Ext_20mm_pressure_for_human_torque.csv", "Flx_20mm_pressure_for_human_torque.csv"]:
    rows = list(csv.DictReader(open(BASE + "\\" + f)))
    print(f, len(rows), "rows; cols", list(rows[0].keys()))
    for col in rows[0]:
        if "pressure" in col:
            vals = [float(r[col]) for r in rows if r[col] not in ("", "nan")]
            print("  ", col, "min %.1f max %.1f" % (min(vals), max(vals)))
    ht = [float(r["human_target_Nm"]) for r in rows]
    ang = [float(r["knee_angle_deg"]) for r in rows]
    print("   human_target min %.3f max %.3f; angle %.1f..%.1f" % (min(ht), max(ht), min(ang), max(ang)))
