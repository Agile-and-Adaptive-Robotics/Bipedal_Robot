import os
import subprocess

GD = r"C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.7\resources\app\git\cmd\git.exe"
ROOT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022"
p = os.path.join(ROOT, r"patches\30-torque\02-publication-note.patch")
TMP = os.path.join(os.environ["TEMP"], "patchchk")

data = open(p, "rb").read()
crlf = data.replace(b"\n", b"\r\n")
open(p, "wb").write(crlf)
r = subprocess.run(
    [GD, "apply", "--check", p], cwd=TMP, capture_output=True, text=True
)
print("CRLF retry:", "OK" if r.returncode == 0 else r.stderr.strip()[:200])
