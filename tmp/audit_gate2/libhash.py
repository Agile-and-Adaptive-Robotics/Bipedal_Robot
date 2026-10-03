"""Hash the Simulink XML payload of SNS_Library.slx (ignoring metadata/time)."""
import hashlib
import sys
import zipfile

path = r"D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_SNS.slsx"  # placeholder guard
path = r"D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\SNS_Library.slx"
h = hashlib.sha256()
block_count = 0
with zipfile.ZipFile(path) as z:
    for n in sorted(z.namelist()):
        if n.startswith("simulink/") and n.endswith(".xml"):
            h.update(n.encode())
            h.update(z.read(n))
            data = z.read(n)
            block_count += data.count(b"<Block ")
print(f"payload sha256: {h.hexdigest()}")
print(f"<Block occurrences across simulink/*.xml: {block_count}")
