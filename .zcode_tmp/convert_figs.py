import os, subprocess

GS = r"D:\Anaconda\envs\gs\Library\bin\gswin64c.exe"
FIGS = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\figs"
OUT  = r"D:\GitHub\Bipedal_Robot\.zcode_tmp\deck_assets\figs"
os.makedirs(OUT, exist_ok=True)

JOBS = [
    ("Aim1/01_testJigs1.pdf",              "testJigs_force.png"),
    ("Aim1/02a_MaxForce10-eps-converted-to.pdf", "MaxForce10.png"),
    ("Aim1/02b_MaxForce20-eps-converted-to.pdf", "MaxForce20.png"),
    ("Aim1/03a_MaxForce10_all-eps-converted-to.pdf", "MaxForce10_all.png"),
    ("Aim1/03b_kmax10-eps-converted-to.pdf",     "kmax10.png"),
    ("Aim1/04a_FStar10-eps-converted-to.pdf",    "FStar10.png"),
    ("Aim1/04b_FStar20-eps-converted-to.pdf",    "FStar20.png"),
    ("Aim1/05_Comparison-eps-converted-to.pdf",  "Comparison.png"),
    ("Aim2/FlxPin_group.eps",  "FlxPin_group.png"),
    ("Aim2/ExtPin_group.eps",  "ExtPin_group.png"),
    ("Aim2/Ext10mm_52cm.eps",  "Ext10mm_52cm.png"),
    ("Aim2/FlxGrp.eps",        "FlxGrp.png"),
    ("Aim2/KneeICR.eps",       "KneeICR.png"),
    ("Aim2/testJigs.pdf",      "testJigs_knee.png"),
    ("Aim2/bktFrame.pdf",      "bktFrame.png"),
    ("Preliminary/animatlab_phase1_joint_angles.pdf", "animatlab_phase1_joint_angles.png"),
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
