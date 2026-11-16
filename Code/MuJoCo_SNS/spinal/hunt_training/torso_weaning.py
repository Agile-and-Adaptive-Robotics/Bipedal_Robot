"""Run early torso-support curriculum rungs and stop at the first failed gate."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def run(base, output, scales):
    output.mkdir(parents=True, exist_ok=True)
    records = []
    for stage in (1, 2):
        for scale in scales:
            config = dict(base, stage=stage, support_attachment="torso",
                          torso_guide_scale=scale)
            name = f"stage{stage}_torso_{scale:g}".replace(".", "p")
            path = output / f"{name}.json"
            result = output / f"{name}.result.json"
            if result.exists():
                row = json.loads(result.read_text(encoding="utf-8"))
                if not all(row["configuration"].get(k) == value for k, value in config.items()):
                    raise ValueError(f"Checkpoint settings differ from this run: {result}")
            else:
                if path.exists() and json.loads(path.read_text(encoding="utf-8")) != config:
                    raise ValueError(f"Pending configuration differs from this run: {path}")
                path.write_text(json.dumps(config, indent=2), encoding="utf-8")
                with (output / f"{name}.log").open("a", encoding="utf-8") as log:
                    completed = subprocess.run([sys.executable, str(HERE / "physical_eval.py"), str(path)],
                                               stdout=log, stderr=subprocess.STDOUT, check=False)
                if completed.returncode or not result.exists():
                    return records + [dict(stage=stage, scale=scale, accepted=False,
                                          reason="evaluator error", exit_code=completed.returncode)]
                row = json.loads(result.read_text(encoding="utf-8"))
            assessment = row["assessment"]
            accepted = bool(assessment.get("accepted") and assessment.get("trunk_support_gate"))
            records.append(dict(stage=stage, scale=scale, accepted=accepted, result=str(result),
                                trunk_gate=assessment.get("trunk_support_gate")))
            if not accepted:
                return records
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scales", type=float, nargs="+", default=(1., .8, .6))
    args = parser.parse_args()
    if not 0 <= args.scales[-1] <= args.scales[0] <= 1 or any(
            a <= b for a, b in zip(args.scales, args.scales[1:])):
        parser.error("Use strictly descending torso guide scales in [0, 1]")
    base = json.loads(args.config.read_text(encoding="utf-8-sig"))
    records = run(base, args.output.resolve(), args.scales)
    (args.output / "WEANING_RESULT.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
    print(json.dumps(records, indent=2))
    return 0 if all(row["accepted"] for row in records) else 2


if __name__ == "__main__":
    raise SystemExit(main())
