# Recon part 3: utf-8-safe score grep (reports tail) + provenance of decomp winners.
import os, re, glob, subprocess, sys

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

pat = re.compile(r"(-?18[5-9]\.\d|-?19\d\.\d|-?20[0-9]\.\d|of-record|gait_lib|s3k|syn6|pruned|winner)", re.I)
for d in ["reports_20260925", "reports_spiking_20261002", "campaigns/20260930"]:
    for p in sorted(glob.glob(d + "/**/*.md", recursive=True)):
        try:
            txt = open(p, encoding="utf-8", errors="replace").read()
        except Exception:
            continue
        hits = [ln.strip() for ln in txt.splitlines() if pat.search(ln)]
        if hits:
            print(f"\n--- {p} ({len(hits)} hit lines) ---")
            for h in hits[:20]:
                print("   ", h[:200].encode("ascii", "replace").decode("ascii"))

print("\n=== GIT PROVENANCE ===")
for f in ["decomp_winner_s3k_mutfix.json", "decomp_winner_syn6.json",
          "decomp_winner_pruned.json", "GAIT_DECOMPOSITION_20261006.json",
          "gait_refs/ong_speed_050.npz"]:
    r = subprocess.run(["git", "log", "--oneline", "-3", "--", f],
                       capture_output=True, text=True)
    tracked = subprocess.run(["git", "ls-files", "--", f],
                             capture_output=True, text=True).stdout.strip()
    print(f"{f}: tracked={bool(tracked)} log: {r.stdout.strip()[:200] or '(no commits)'}")
