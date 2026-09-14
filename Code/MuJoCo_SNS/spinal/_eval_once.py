"""Run one runner eval and print the metrics dict as JSON.

Usage: python _eval_once.py --fitted --best [--phase-reset E F] ...
Passes all args to runner.main; prints one JSON line 'METRICS {...}'.
"""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import runner as R

if __name__ == "__main__":
    m = R.main(sys.argv[1:])
    m["kine_score"] = float(m["kine_score"])
    k = m.pop("kine", None)
    print("METRICS " + json.dumps({"metrics": m, "kine": k}))
