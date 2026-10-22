import glob
import os
import shutil
import subprocess
import difflib

GD = r"C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.7\resources\app\git\cmd\git.exe"
ROOT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022"
TMP = os.path.join(os.environ["TEMP"], "patchchk")

chapter = "30-torque"  # smallest case
if os.path.isdir(TMP):
    shutil.rmtree(TMP)
os.makedirs(os.path.join(TMP, "chapters"))
shutil.copy(
    os.path.join(ROOT, "baseline", "chapters", chapter + ".tex"),
    os.path.join(TMP, "chapters", chapter + ".tex"),
)
for sp in sorted(glob.glob(os.path.join(ROOT, "patches", chapter, "*.patch"))):
    r = subprocess.run([GD, "apply", "--unidiff-zero", sp], cwd=TMP, capture_output=True, text=True)
    print(os.path.basename(sp), "rc", r.returncode, r.stderr[:100])

got = open(os.path.join(TMP, "chapters", chapter + ".tex"), encoding="utf-8").read().splitlines()
want = open(os.path.join(ROOT, "proposals", chapter + ".tex"), encoding="utf-8").read().splitlines()
d = list(difflib.unified_diff(got, want, lineterm="", n=0))
print("diff lines:", len(d))
print("\n".join(d[:40]))
