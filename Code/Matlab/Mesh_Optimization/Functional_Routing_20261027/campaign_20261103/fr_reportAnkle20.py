"""Report 20 mm distal optimization, dense pose checks, and assembly fit."""
import csv
import json
import sys
from pathlib import Path

from fr_exportPilotOpenSim import export_edl_pair
from fr_plotDistalStages import plot


def report(folder, geometry):
    plot(folder, geometry)
    records = [json.loads(path.read_text()) for path in sorted(folder.glob("*.json"))]
    assert len(records) == 20 and all(row["diameter_mm"] == 20 for row in records)
    rows = []
    for r in records:
        name = f"{r['group']}_{r['branch']}"
        e = r.get("coupled")
        peak = None
        joint, second = "", None
        if e is None:
            ankle = 100*r["fractionMet"]
            no_force = 100*r["noForceFraction"]
            failed = [c["name"] for c in r["constraints"] if c["violated"]]
            bone = -1000*next(c["cDim"] for c in r["constraints"] if c["name"] == "bone_exclusion")
        else:
            with (folder/f"{name}_coupled_torque.csv").open(newline="") as stream:
                grid = list(csv.DictReader(stream))

            def met_fraction(axis):
                hit = 0
                for p in grid:
                    human = float(p[f"human_{axis}_Nm"])
                    robot = float(p[f"robot_{axis}_Nm"])
                    signed = robot*((human > 0)-(human < 0))
                    hit += p["no_force"] == "False" and signed >= abs(human)-1e-6
                return 100*hit/len(grid)

            ankle, second = met_fraction("ankle"), met_fraction("second")
            joint = "MTP" if e["second"] == "mtp_angle_r" else "subtalar"
            no_force = 100*sum(e["noForce"])/len(e["noForce"])
            bone = 1000*e["boneMargin_m"]
            peak = max((x for x in e["tendonForce_N"] if x is not None), default=None)
            margins = dict(boneMargin_m="bone", leadMargin_m="tibial tendon lead",
                           footMargin_m="vessel above ankle", distalGuideMargin_m="sleeve bay",
                           wrapMargin_m="vessel bend", shaftCableMargin_m="shaft tendon",
                           bundleMargin_m="within-group sleeve", bundleHardwareMargin_m="within-group fitting")
            failed = [label for key, label in margins.items()
                      if e.get(key) is not None and e[key] < (-1e-9 if key == "bundleMargin_m" else -1e-4)]
            if no_force > 0:
                failed.append("no-force poses")
            if r["tendon_mm"]/1000 < e["tendonRequired_m"]-1e-6:
                failed.append("minimum tendon span")
            if r.get("layer") and not r["layer"]["layerPass"]:
                failed.append("EHL layer")
        rows.append(dict(group=r["group"], branch=r["branch"], BPA_count=r["count"],
                         sheaves=r["topologySheaves"], gain=r["gain"], rest_mm=r["rest_mm"],
                         tendon_mm=r["tendon_mm"], evaluations=r["trialEvaluations"],
                         ankle_met_pct=ankle, second_joint=joint, second_met_pct=second,
                         no_force_pct=no_force, bone_margin_mm=bone, peak_per_tendon_N=peak,
                         local_pass=r["localConstraintsPass"], failures=", ".join(failed)))

    with (folder/"ankle20_group_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    lines = ["# Tibia-down 20 mm routing verification", "",
             "Each group and branch screened independent 20 mm BPAs. Identical load paths share the "
             "geometry calculation and sum their wrenches. The general packaging screen spaces the "
             "vessels about the representative route. The fixed-length EDL trial solves and exports "
             "two distinct tendons and mount positions. The EDL rows below refer to the general "
             "optimization search. The separately routed pair is used in the compliant assembly "
             "and its OpenSim preview.", "",
             "The table uses 21 × 21 coupled poses, with a 31-point ankle-only check for soleus. "
             "Sheave counts above zero describe ideal topologies.", "",
             "| Group | Branch | BPAs | Sheaves | Ankle met | Second joint met | No force | Bone (mm) | Local screen |",
             "|---|---|---:|---:|---:|---:|---:|---:|---|"]
    for r in rows:
        second = f"{r['second_joint']} {r['second_met_pct']:.1f}%" if r["second_joint"] else "—"
        lines.append(f"| {r['group']} | {r['branch']} | {r['BPA_count']} | {r['sheaves']} | "
                     f"{r['ankle_met_pct']:.1f}% | {second} | {r['no_force_pct']:.1f}% | "
                     f"{r['bone_margin_mm']:.2f} | {r['failures'] or ('pass' if r['local_pass'] else 'fail')} |")
    lines.extend(["", "## Simultaneous rigid hardware check", "",
                  "The compliant assembly uses the separately routed fixed-length EDL pair.", ""])
    for branch in ("rigid", "compliant"):
        file = folder/f"assembled_{branch}.csv"
        if not file.is_file():
            lines.append(f"- {branch.capitalize()} check is pending.")
            continue
        with file.open(newline="") as stream:
            pairs = list(csv.DictReader(stream))
        failures = sum(p["rigid_pass"] not in {"1", "true", "True"} for p in pairs)
        worst = min(pairs, key=lambda p: float(p["rigid_margin_mm"]))
        lines.append(f"- {branch.capitalize()} has {failures} of {len(pairs)} group pairs below the "
                     f"rigid hardware requirement. The smallest margin is {float(worst['rigid_margin_mm']):.2f} mm "
                     f"between {worst['group_a']} and {worst['group_b']}.")
    (folder/"ANKLE20_VERIFICATION.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    pair = folder/"edl_fixed_pair_separate_guides"
    recheck = pair/"final_ehl_recheck.json"
    if recheck.is_file():
        comparison = json.loads(recheck.read_text())
        if comparison["pairedPass"]:
            selected = folder/"edl_with_selected_ehl"
            selected.mkdir(exist_ok=True)
            r = json.loads((pair/"ext_dig_compliant.json").read_text())
            ehl = json.loads((folder/"ext_hal_compliant.json").read_text())
            r["innerRoute"] = {key: ehl[key] for key in ("points_m","bodies")}
            r["ehlLocalPass"] = comparison["ehlLocalPass"]
            r["assembledLocalPass"] = comparison["pairedPass"]
            for key in ("rigidMargin_m","orderingMargin_m"):
                r["layer"][key] = min(check[key] for check in comparison["clearance"])
            r["layer"]["sleeveOverlap_m"] = max(check["sleeveOverlap_m"] for check in comparison["clearance"])
            route = selected/"ext_dig_compliant.json"
            route.write_text(json.dumps(r), encoding="utf-8")
            plot(selected, geometry)
            export_edl_pair(route, geometry)


if __name__ == "__main__":
    report(Path(sys.argv[1]), Path(sys.argv[2]))
