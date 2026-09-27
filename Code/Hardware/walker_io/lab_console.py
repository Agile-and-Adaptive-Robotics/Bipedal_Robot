"""SKELETON: lab-side operator console — computer-in-the-loop.
Sends deletions / stimulus injections / manual valve overrides to the Orin,
receives telemetry, writes trial logs to Testing_Data\\walker_trials\\.

Usage (once the protocol is live):
    python lab_console.py --trial G5_del_ia_01
    then at the prompt:  del ia_in=0      |  stim HC-RG-E_R 2.0 0.2
                         valve VASTI_R 0.4 |  stop
"""
from __future__ import annotations
import argparse, os, time
from protocol import UdpLink, Telemetry, Command, CMD_PORT, TELEM_PORT, ORIN_IP


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", required=True, help="trial id, e.g. G5_del_ia_01")
    ap.add_argument("--outdir", default=r"Testing_Data\walker_trials")
    args = ap.parse_args()

    outdir = os.path.join(args.outdir, time.strftime("%Y%m%d") + "_" + args.trial)
    os.makedirs(outdir, exist_ok=True)
    log = open(os.path.join(outdir, "telemetry.jsonl"), "ab")
    cmdlog = open(os.path.join(outdir, "commands.jsonl"), "ab")

    cmds = UdpLink(CMD_PORT, remote=(ORIN_IP, CMD_PORT))
    telem = UdpLink(TELEM_PORT)

    print(f"trial {args.trial} -> {outdir}")
    print("commands: del GAIN=0 [...] | stim NEURON I DUR | valve NAME DUTY | stop")
    running = True
    while running:
        raw = input("> ").strip()
        if not raw:
            continue
        toks = raw.split()
        c: Command | None = None
        if toks[0] == "del":
            g = dict(t.split("=") for t in toks[1:])
            c = Command.override_gains({k: float(v) for k, v in g.items()})
        elif toks[0] == "stim":
            c = Command.stim(toks[1], float(toks[2]), float(toks[3]))
        elif toks[0] == "valve":
            c = Command.valve_override(toks[1], float(toks[2]))
        elif toks[0] == "stop":
            c = Command.stop()
            running = False
        if c:
            cmds.send(c.to_json())
            cmdlog.write(c.to_json() + b"\n")
            cmdlog.flush()
        # drain telemetry
        while (m := telem.recv(0.0)) is not None:
            if m.get("k") == "telem":
                log.write(m and (__import__("json").dumps(m) + "\n").encode())
        log.flush()

    log.close()
    cmdlog.close()
    print("trial closed:", outdir)


if __name__ == "__main__":
    main()
