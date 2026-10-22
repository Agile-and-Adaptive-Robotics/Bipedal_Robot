import glob
import os
import shutil
import subprocess

GD = r"C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.7\resources\app\git\cmd\git.exe"
ROOT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022"
TMP = os.path.join(os.environ["TEMP"], "patchchk")

for chapter in ("10-introduction", "15-background", "30-torque"):
    if os.path.isdir(TMP):
        shutil.rmtree(TMP)
    os.makedirs(os.path.join(TMP, "chapters"))
    shutil.copy(
        os.path.join(ROOT, "baseline", "chapters", chapter + ".tex"),
        os.path.join(TMP, "chapters", chapter + ".tex"),
    )
    subs = sorted(glob.glob(os.path.join(ROOT, "patches", chapter, "*.patch")))
    ok = True
    for sp in subs:
        r = subprocess.run(
            [GD, "apply", "--unidiff-zero", sp], cwd=TMP, capture_output=True, text=True
        )
        if r.returncode != 0:
            print(chapter, "FAILED applying", os.path.basename(sp), r.stderr[:150])
            ok = False
            break
    if ok:
        got = open(os.path.join(TMP, "chapters", chapter + ".tex"), "rb").read()
        want = open(os.path.join(ROOT, "proposals", chapter + ".tex"), "rb").read()
        print(chapter, len(subs), "sub-patches compose to proposal:", got == want)
