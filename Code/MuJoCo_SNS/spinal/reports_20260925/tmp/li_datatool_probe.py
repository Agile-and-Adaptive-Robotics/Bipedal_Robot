import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
p = r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\DataTool_7.txt'
with open(p, 'r', errors='replace') as fh:
    lines = fh.readlines()
print("lines:", len(lines))
for i, ln in enumerate(lines[:8]):
    print(f"{i:3d}: {ln.rstrip()[:250]}")
print("...")
for i, ln in enumerate(lines[-4:], start=len(lines) - 4):
    print(f"{i:3d}: {ln.rstrip()[:250]}")
