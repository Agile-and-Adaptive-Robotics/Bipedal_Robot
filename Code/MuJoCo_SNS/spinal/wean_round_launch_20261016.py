"""Launch N parallel optuna workers for one wean round + extract best.

Usage: python wean_round_launch_20261016.py <r> <n_per_worker> <nworkers>
Runs from schtasks (survives ssh drops — the drive-sweep pattern). Each
worker = subprocess `wean_tune_20261016.py round <r> <n> <w>` with its own
log round{r}_w{w}.log; after all finish, runs `best <r>`. Writes master
log round{r}.log in reports_20261016/wean_human/.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reports_20261016", "wean_human")
PY = r"D:\Anaconda\envs\myo\python.exe"

r = sys.argv[1]
n = sys.argv[2]
nw = int(sys.argv[3])
log = open(os.path.join(OUT, f"round{r}.log"), "a", encoding="utf-8",
           buffering=1)


def say(msg):
    log.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n")


say(f"LAUNCH round {r} x{n} trials, {nw} workers")
procs = []
for w in range(nw):
    wl = open(os.path.join(OUT, f"round{r}_w{w}.log"), "a",
              encoding="utf-8", buffering=1)
    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    p = subprocess.Popen([PY, "-u", os.path.join(HERE,
                                                 "wean_tune_20261016.py"),
                          "round", r, n, str(w)], stdout=wl, stderr=wl,
                         cwd=HERE, env=env)
    procs.append((w, p, wl))
    say(f"worker {w} pid {p.pid}")
    time.sleep(20)          # stagger: sqlite create-study + seed race guard
rcs = []
for w, p, wl in procs:
    rc = p.wait()
    wl.close()
    rcs.append(rc)
    say(f"worker {w} EXIT {rc}")
say(f"all workers done rcs={rcs}")
if all(rc == 0 for rc in rcs):
    b = subprocess.run([PY, "-u", os.path.join(HERE,
                                               "wean_tune_20261016.py"),
                        "best", r], cwd=HERE, capture_output=True,
                       text=True, env={**os.environ,
                                       "PYTHONIOENCODING": "utf-8"})
    say("BEST stdout: " + b.stdout.strip())
    if b.returncode != 0:
        say("BEST stderr: " + b.stderr[-2000:])
say(f"ROUND{r}_DONE")
log.close()
