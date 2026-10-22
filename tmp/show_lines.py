import sys

path = sys.argv[1]
lo = int(sys.argv[2])
hi = int(sys.argv[3])
lines = open(path, encoding="utf-8-sig", errors="replace").read().split("\n")
for i in range(lo - 1, min(hi, len(lines))):
    t = lines[i].strip()
    if not t or t.startswith("%"):
        continue
    print(f"--- line {i+1} ---")
    print(t[:2200])
