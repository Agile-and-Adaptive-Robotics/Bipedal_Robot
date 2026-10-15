"""Render skeletal walking gifs for all model winners using the MuJoCo offscreen renderer."""
import subprocess
import sys
from pathlib import Path

SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
MYO = r"C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
SLIDESHOW = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\slideshow")
MODEL_XML = r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml"
W2L_XML = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_mjcf_fixed.xml"

MODELS = [
    ("s3k_mutfix", "spinal_run_s3k_mutfix_best_20261003.npz", MODEL_XML),
    ("syn6", "spinal_run_syn6_best_20261003.npz", MODEL_XML),
    ("pruned", "spinal_run_pruned_best_20261003.npz", MODEL_XML),
    ("w2lvar", "spinal_run_w2lvar_s1_best_20261003.npz", MODEL_XML),
    ("spiking_mirror", "spinal_run_spk_s1_best_20261003.npz", MODEL_XML),
]

for label, npz_name, xml in MODELS:
    npz = SPINAL / npz_name
    gif = SLIDESHOW / f"{label}_skeletal.gif"
    if not npz.exists():
        print(f"SKIP {label}: {npz_name} not found")
        continue
    print(f"Rendering {label}...")
    r = subprocess.run(
        [MYO, str(SPINAL / "render_walking_gif.py"), str(npz), str(gif), xml],
        cwd=str(SPINAL), capture_output=True, text=True, timeout=300,
        encoding="utf8", errors="replace"
    )
    if r.returncode == 0:
        print(f"  ✓ {gif.name}")
    else:
        print(f"  ✗ {r.stdout[-200:]} {r.stderr[-200:]}")

# BPA robot gifs (different model)
for label, npz_name in [("bpa_nonspiking", "bpa_ns_best_20261003.npz"),
                          ("bpa_spiking", "bpa_spk_best_20261003.npz")]:
    npz = SPINAL / npz_name
    gif = SLIDESHOW / f"{label}_skeletal.gif"
    if not npz.exists():
        print(f"SKIP {label}: {npz_name} not found")
        continue
    # BPA walker uses its own MJCF
    bpa_xml = SPINAL / "w2l_mujoco" / "w2l_mjcf_bpa.xml"
    if not bpa_xml.exists():
        print(f"SKIP {label}: BPA model xml not found")
        continue
    print(f"Rendering {label}...")
    r = subprocess.run(
        [MYO, str(SPINAL / "render_walking_gif.py"), str(npz), str(gif), str(bpa_xml)],
        cwd=str(SPINAL), capture_output=True, text=True, timeout=300,
        encoding="utf8", errors="replace"
    )
    if r.returncode == 0:
        print(f"  ✓ {gif.name}")
    else:
        print(f"  ✗ {r.stdout[-200:]} {r.stderr[-200:]}")

# Remove the old stick-figure gifs
removed = 0
for old in SLIDESHOW.glob("*_walk.gif"):
    old.unlink()
    removed += 1
    print(f"  Removed old stick-figure: {old.name}")
print(f"\nDone. Removed {removed} stick-figure gifs; skeletal gifs in {SLIDESHOW}")
