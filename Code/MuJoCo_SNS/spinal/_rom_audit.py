"""Inspect driver-joint limit attributes in the converted MJCF."""
import io
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
from pathlib import Path

import runner as R

txt = Path(R.MODEL).read_text(encoding="utf-8")
names = ["knee_angle", "hip_flexion", "ankle_angle", "subtalar_angle",
         "mtp_angle", "hip_adduction", "hip_rotation", "lumbar_extension",
         "pelvis_tilt"]
for base in names:
    for side in ("r", "l"):
        m = re.search(rf'<joint name="{base}_{side}"[^/>]*/>', txt)
        if m:
            tag = m.group(0)
            lim = re.search(r'limited="([^"]+)"', tag)
            rng = re.search(r'range="([^"]+)"', tag)
            print(f"{base}_{side:1s}  limited={lim.group(1) if lim else '-':5s}"
                  f"  range={rng.group(1) if rng else 'NONE'}")
