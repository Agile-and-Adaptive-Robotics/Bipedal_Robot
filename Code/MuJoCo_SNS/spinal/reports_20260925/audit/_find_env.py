import os

os.chdir(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
for f in sorted(os.listdir(".")):
    if f.endswith(".py"):
        try:
            lines = open(f, encoding="utf-8", errors="replace").read().splitlines()
        except Exception as e:
            print(f, "ERR", e)
            continue
        for i, l in enumerate(lines, start=1):
            if "AARL_NPZ" in l or "AARL_NET" in l:
                print(f"{f}:{i}: {l.strip()}")
