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
        r = subprocess.run([GD, "apply", "--unidiff-zero", sp], cwd=TMP, capture_output=True, text=True)
        if r.returncode != 0:
            print(chapter, "FAILED", os.path.basename(sp), r.stderr[:150])
            ok = False
            break
    if not ok:
        continue
    got = open(os.path.join(TMP, "chapters", chapter + ".tex"), "rb").read()
    want = open(os.path.join(ROOT, "proposals", chapter + ".tex"), "rb").read()
    same = got == want
    print(chapter, len(subs), "sub-patches; byte-identical:", same)
    if not same:
        print("  got head:", got[:20], "want head:", want[:20])
        print("  got len", len(got), "want len", len(want))
        # first differing byte
        for i in range(min(len(got), len(want))):
            if got[i] != want[i]:
                print("  first diff at", i, got[i - 40 : i + 20], b"|", want[i - 40 : i + 20])
                break
