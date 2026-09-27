# walker_io — Orin Nano walker I/O + computer-in-the-loop skeletons

SKELETON code for the G5 physical test (see
`Documentation\Program_Workflow\BPA_Walker_Program.md` for the master plan).

| File | Runs on | Purpose |
|---|---|---|
| `protocol.py` | both | message schemas + UDP endpoints (the wire contract) |
| `valves.py` | Orin | valve driver abstraction: bang-bang + taped-proportional, watchdog |
| `sensors.py` | Orin | pressure / encoder / Liquid Wire / IMU / insole readers (hardware TODOs inside) |
| `afferents.py` | both | sensor → SNS afferent encoding (real formulas; constants from `bpa_actuators_27.json`) |
| `orin_daemon.py` | Orin | 1 kHz loop: sensors → afferents → SNS → valves; telemetry; safety watchdogs |
| `lab_console.py` | lab PC | operator console: deletions, stimulus injection, live plots |

Conventions:
- The SNS itself is NOT in these files — the Orin loads the same spinal network config
  (params.G json) the MuJoCo/Simscape hosts use. `sns_rt.py` hook in `orin_daemon.py` is where
  the compiled/exported SNS step function plugs in.
- All hardware access points are marked `# TODO(hw)` — bench-bring-up fills them one at a time.
- Safety logic lives in `orin_daemon.py` (`watchdogs`) and is deliberately dumb.
