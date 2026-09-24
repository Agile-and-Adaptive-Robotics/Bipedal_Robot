"""Dump mtimes (ISO, UTC) of chart .txt files in the Walker_2_Layer_CPG folder.
Usage: _w2l_mtimes.py  (prints JSON to stdout; caller redirects to a file)
"""
import json
import os
import datetime

FOLDER = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG"
out = {}
for name in sorted(os.listdir(FOLDER)):
    if name.lower().endswith(".txt"):
        p = os.path.join(FOLDER, name)
        st = os.stat(p)
        out[name] = {
            "mtime_local": datetime.datetime.fromtimestamp(st.st_mtime).isoformat(),
            "size": st.st_size,
        }
print(json.dumps(out, indent=1))
