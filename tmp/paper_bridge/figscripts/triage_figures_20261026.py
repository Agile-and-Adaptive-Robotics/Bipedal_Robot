import os
import shutil

root = "Documentation/Reports and Papers/Knee_Torque_Test/Figures"
comp = os.path.join(root, "Figure components")
dep = os.path.join(root, "deprecated")
os.makedirs(dep, exist_ok=True)

keep = {
    "testJigs.pdf", "KneeICR.eps", "KneeICR-eps-converted-to.pdf",
    "xiFrameTransform.pdf", "bktFrame.pdf",
    "ExtPin_group.eps", "ExtPin_group-eps-converted-to.pdf",
    "ExtensorRoute_redesigned.eps", "ExtensorRoute_redesigned-eps-converted-to.pdf", "ExtensorRoute_redesigned.pdf",
    "RouteVariants_20261020.eps", "RouteVariants_20261020-eps-converted-to.pdf", "RouteVariants_20261020.pdf",
    "FlxPin_allmodels_20261023.eps", "FlxPin_allmodels_20261023-eps-converted-to.pdf",
    "Ext10mm_52cm_20261023.eps", "Ext10mm_52cm_20261023-eps-converted-to.pdf",
    "FlxGrp_20261023.eps", "FlxGrp_20261023-eps-converted-to.pdf",
    "Ext20mm_cam_20261023.eps", "Ext20mm_cam_20261023-eps-converted-to.pdf",
    "LoadCellCalibration.png", "RobotSim_frames.pdf", "extensor_misalign.png",
}
components = {
    "KneeExt_10mm_passthru.png",
    "Flx20mm_56cm.eps", "Flx20mm_56cm.pdf", "Flx20mm_56cm-eps-converted-to.pdf",
}

present = set(os.listdir(root))
missing_keep = keep - present
missing_comp = components - present
unclassified = sorted(present - keep - components - {"Figure components", "deprecated"})

moved_dep, moved_comp = [], []
for f in sorted(present):
    src = os.path.join(root, f)
    if os.path.isdir(src):
        continue
    if f in components:
        shutil.move(src, os.path.join(comp, f)); moved_comp.append(f)
    elif f not in keep:
        shutil.move(src, os.path.join(dep, f)); moved_dep.append(f)

print(f"kept at root: {len(keep & present)}")
print(f"moved to Figure components: {moved_comp}")
print(f"moved to deprecated: {len(moved_dep)} files")
print(f"missing-but-expected keep: {sorted(missing_keep)}")
print(f"missing-but-expected comp: {sorted(missing_comp)}")
print(f"unclassified leftovers: {unclassified}")
