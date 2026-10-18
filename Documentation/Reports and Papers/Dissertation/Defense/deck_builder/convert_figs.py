import os, subprocess

GS = r"D:\Anaconda\envs\gs\Library\bin\gswin64c.exe"
ROOT = os.path.dirname(os.path.abspath(__file__))
FIGS = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures"
OUT  = os.path.join(ROOT, "deck_assets", "figs")
os.makedirs(OUT, exist_ok=True)

JOBS = [
    ("20-muscles/01_testJigs1.pdf",              "testJigs_force.png"),
    ("20-muscles/02a_MaxForce10-eps-converted-to.pdf", "MaxForce10.png"),
    ("20-muscles/02b_MaxForce20-eps-converted-to.pdf", "MaxForce20.png"),
    ("20-muscles/03a_MaxForce10_all-eps-converted-to.pdf", "MaxForce10_all.png"),
    ("20-muscles/03b_kmax10-eps-converted-to.pdf",     "kmax10.png"),
    ("20-muscles/04a_FStar10-eps-converted-to.pdf",    "FStar10.png"),
    ("20-muscles/04b_FStar20-eps-converted-to.pdf",    "FStar20.png"),
    ("20-muscles/05_Comparison-eps-converted-to.pdf",  "Comparison.png"),
    ("30-torque/FlxPin_group-eps-converted-to.pdf",  "FlxPin_group.png"),
    ("30-torque/ExtPin_group-eps-converted-to.pdf",  "ExtPin_group.png"),
    ("30-torque/Ext10mm_52cm-eps-converted-to.pdf",  "Ext10mm_52cm.png"),
    ("30-torque/FlxGrp-eps-converted-to.pdf",        "FlxGrp.png"),
    ("30-torque/KneeICR-eps-converted-to.pdf",       "KneeICR.png"),
    ("30-torque/testJigs.pdf",      "testJigs_knee.png"),
    ("30-torque/bktFrame.pdf",      "bktFrame.png"),
    ("60-sns/animatlab_phase1_joint_angles.pdf", "animatlab_phase1_joint_angles.png"),
    # batch C (Ch.5/6 + close): c_-prefixed rasters for slides_c.py
    ("60-sns/animatlab_phase1_joint_angles.pdf", "c_animatlab_phase1_joint_angles.png"),
]
fails = []
for src, dst in JOBS:
    s = os.path.join(FIGS, src)
    d = os.path.join(OUT, dst)
    r = subprocess.run([GS, "-dBATCH", "-dNOPAUSE", "-dQUIET", "-sDEVICE=png16m",
                        "-r220", "-dEPSCrop", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                        f"-sOutputFile={d}", s], capture_output=True, text=True)
    ok = os.path.exists(d) and os.path.getsize(d) > 5000
    print(("OK  " if ok else "FAIL"), dst, "" if ok else r.stderr[-200:])
    if not ok: fails.append(src)
print("FAILS:", fails)
