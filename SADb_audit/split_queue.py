"""Split curation_queue/queue.json into per-batch files of N papers (2026-09-28)."""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
Q = os.path.join(HERE, "curation_queue")
queue = json.load(open(os.path.join(Q, "queue.json"), encoding="utf-8"))
N = 25
n = 0
for i in range(0, len(queue), N):
    slug = f"{i//N:02d}"
    json.dump(queue[i:i + N], open(os.path.join(Q, f"queue_{slug}.json"), "w"),
              ensure_ascii=False, indent=1)
    n += 1
print(f"{len(queue)} papers -> {n} batch files of <= {N}")
