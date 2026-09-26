import io, sys, glob, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

f = sorted(glob.glob(r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\Trace_2020*.txt'))[0]
with open(f, 'r', errors='replace') as fh:
    lines = fh.readlines()
print("total lines:", len(lines))
for i, ln in enumerate(lines[:40]):
    print(f"{i:5d}: {ln.rstrip()[:160]}")
print("...")
for i, ln in enumerate(lines[-15:], start=len(lines) - 15):
    print(f"{i:5d}: {ln.rstrip()[:160]}")
