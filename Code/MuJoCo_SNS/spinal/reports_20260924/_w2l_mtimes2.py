"""Dump mtimes (ISO, local) of chart .txt files in a model folder.
Usage: _w2l_mtimes2.py <folder>  (prints JSON to stdout; caller redirects)
"""
import json
import os
import sys
import datetime

FOLDER = sys.argv[1]
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
