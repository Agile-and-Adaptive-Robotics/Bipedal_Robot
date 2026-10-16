"""Extract the best-trial params from each study's curriculum json and save as replayable winner files."""
import json
from pathlib import Path

SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")

SOURCES = {
    "s3k_mutfix": "curriculum_stage4.json",  # generic stage-4 params (the runner knows the study)
    "syn6": "curriculum_syn6_stage5_reseed_20261003.json",
    "pruned": "curriculum_pruned_walk_2_fix_20261003.json",
    "w2lvar": "curriculum_w2lvar_stage3_fix_20261003.json",
    "bpa_ns": "curriculum_bpa_ns.json",
    "bpa_spk": "curriculum_bpa_spk.json",
}

for label, fname in SOURCES.items():
    p = SPINAL / fname
    if not p.exists():
        print(f"{label}: {fname} NOT FOUND")
        continue
    d = json.loads(p.read_text())
    print(f"{label}: keys={list(d.keys())[:6]}")
    # Save as a winner file the runner can use
    out = SPINAL / f"decomp_winner_{label}.json"
    out.write_text(json.dumps(d, indent=2))
    print(f"  -> wrote {out.name}")
