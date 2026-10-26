# Apply the three proven bridge-build preconditions to the desktop's
# blockset copy (per BRIDGE_REPORT.md + laptop_build_bridge.m):
# 1. install.m MJ_VER -> '2.3.7'
# 2. tools/setupBuild.m -> selectedCompilerWin="MINGW"
# 3. AARL LOCAL PATCH in src/mj.cpp: after mj_makeData in initData(),
#    reset from keyframe 0 (nkey>0) else mj_resetData(m, d).
from pathlib import Path
import re

R = Path(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\mujoco_bridge"
         r"\blockset\mujoco-simulink-blockset-main")

# 1
p = R / "install.m"
s = p.read_text(encoding="utf-8", errors="replace")
s2 = re.sub(r"MJ_VER\s*=\s*'[^']*'", "MJ_VER = '2.3.7'", s)
print("install.m MJ_VER set:", s2 != s)
p.write_text(s2, encoding="utf-8")

# 2
p = R / "tools" / "setupBuild.m"
s = p.read_text(encoding="utf-8", errors="replace")
s2 = re.sub(r'selectedCompilerWin\s*=\s*"[^"]*"',
            'selectedCompilerWin="MINGW"', s)
print("setupBuild MINGW set:", s2 != s)
p.write_text(s2, encoding="utf-8")

# 3
p = R / "src" / "mj.cpp"
s = p.read_text(encoding="utf-8", errors="replace")
if "AARL LOCAL PATCH" in s:
    print("mj.cpp already patched")
else:
    m = re.search(r"(mj_makeData\s*\([^;]*;)", s)
    assert m, "mj_makeData call not found in mj.cpp"
    patch = (m.group(1) +
             "\n    // AARL LOCAL PATCH: bare mj_makeData zero-fills qpos;"
             "\n    // start from keyframe 0 (the standing pose) when the"
             "\n    // model defines one, else the model default state."
             "\n    if (m->nkey > 0) { mj_resetDataKeyframe(m, d, 0); }"
             "\n    else { mj_resetData(m, d); }")
    s2 = s[:m.start(1)] + patch + s[m.end(1):]
    p.write_text(s2, encoding="utf-8")
    print("mj.cpp patched after mj_makeData")
