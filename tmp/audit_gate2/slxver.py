"""Audit: read the Simulik save version from .slx (zip) files."""
import re
import sys
import zipfile

FILES = [
    (r"D:\GitHub\Bipedal_Robot\tmp\audit_gate2\SNS_Library_HEAD.slx", "committed HEAD SNS_Library"),
    (r"D:\GitHub\Bipedal_Robot\tmp\audit_gate2\KneeReflexDemo_HEAD.slx", "committed HEAD KneeReflexDemo"),
    (r"D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\SNS_Library.slx", "working SNS_Library"),
    (r"D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\demos\KneeReflexDemo.slx", "working KneeReflexDemo"),
]

pat = re.compile(rb'(SavedInVersion|Release|SimulinkVersion)="([^"]*)"')

for path, label in FILES:
    try:
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist() if n.endswith(".xml")][:8]
            found = set()
            for n in names:
                data = z.read(n)
                for m in pat.finditer(data):
                    found.add(m.group(0).decode())
            print(f"{label}:")
            for f in sorted(found):
                print(f"   {f}")
            if not found:
                print("   (no version attrs in first xml entries)", names)
    except Exception as e:
        print(f"{label}: ERROR {e}")
