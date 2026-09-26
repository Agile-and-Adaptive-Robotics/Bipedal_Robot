import glob
import os

os.chdir(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
         r"\reports_20260925\logs")
print("== .done files (first line = exit code) ==")
for f in sorted(glob.glob("curr_*.done")):
    first = open(f, encoding="utf-8", errors="replace").readline().strip()
    print(f"{f}: {first}")

print("\n== stage-log key lines ==")
for f in sorted(glob.glob("curr_*_s?.log")):
    print(f"--- {f}")
    for line in open(f, encoding="utf-8", errors="replace"):
        s = line.strip()
        if (s.startswith("seeded") or "columns verified" in s
                or s.startswith("== stage") or "pinned AARL" in s
                or s.startswith("== best")):
            print("   " + s[:400])

print("\n== smoke logs key lines ==")
for f in ("smoke_w2lvar_s1.log", "smoke_syn6_s1.log",
          "smoke_w2lvar_s4.log", "smoke_syn6_s4.log"):
    print(f"--- {f}")
    for line in open(f, encoding="utf-8", errors="replace"):
        s = line.strip()
        if ("columns verified" in s or s.startswith("seeded")
                or s.startswith("== stage")):
            print("   " + s[:400])
