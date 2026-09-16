"""Parse stock gait2392_simbody.osim coordinate ranges (radians)."""
import io
import re
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

OSIM = Path(r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim"
            r"\Gait2392_Robotbody\gait2392_simbody.osim")
txt = OSIM.read_text(encoding="utf-8", errors="replace")
for m in re.finditer(
        r"<Coordinate\s+name=\"([^\"]+)\"[^>]*>.*?<range>\s*([-0-9.eE]+)\s+"
        r"([-0-9.eE]+)\s*</range>", txt, flags=re.S):
    print(f"{m.group(1):24s} {float(m.group(2)):+.4f} {float(m.group(3)):+.4f}"
          f"  ({float(m.group(2)) * 57.2958:+.1f} {float(m.group(3)) * 57.2958:+.1f} deg)")
