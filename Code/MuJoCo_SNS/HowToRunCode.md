# HowToRunCode — MuJoCo/SNS spinal walker on EB475WS4 (Spyder + `myo` env)

Written 2026-09-20. Verified against this machine's installs:
Spyder **6.1.5** (base Anaconda at `C:\ProgramData\anaconda3`) and the
spinal env **`myo`** at `C:\Users\Ben Bolen\.conda\envs\myo`
(Python 3.10.21, NumPy 1.22.4, SciPy 1.9.3, MuJoCo 2.3.7, SNS-Toolbox
1.5.2). The `myo` env already carries `spyder-kernels` 3.1.6, which is
the version Spyder 6 requires — **nothing needs to be installed**.

All scripts live in `Code\MuJoCo_SNS\spinal\` and must run with that
folder as the working directory (they open their jsons and the optuna
db `optuna_walk.db` by relative path).

---

## 1. Point Spyder at the `myo` environment (one time)

1. In Spyder: **Tools → Preferences → Application → Python interpreter**
   (left list inside the Preferences dialog).
2. Select **"Use the following Python interpreter:"** and browse to:

   ```
   C:\Users\Ben Bolen\.conda\envs\myo\python.exe
   ```

3. Click **OK**. If Spyder offers to restart the kernel/console, accept;
   otherwise open a new console (**Consoles → New console**).
4. Verify in the console (do not skip this gate):

   ```python
   import sys; print(sys.executable)
   import numpy, mujoco, sns_toolbox
   print(numpy.__version__, mujoco.__version__)
   ```

   Must print the `...\envs\myo\python.exe` path and `1.22.4 2.3.7`.
   Spyder's GUI keeps running on the base install — only the kernels
   switch. That is correct.

5. Once per console session, also run:

   ```python
   import os
   os.environ["CONDA_PREFIX"] = r"C:\Users\Ben Bolen\.conda\envs\myo"
   ```

   (Only needed if a `mujoco` import ever complains. If you later run
   the diagram scripts such as `spinal\draw_circuit.py`, also do
   `os.environ["PATH"] += r";C:\Users\Ben Bolen\.conda\envs\myo\Library\bin"`
   so graphviz is found.)

**Do not** let Spyder create or "fix" any environment — no further
packages are needed in `myo`, and that env's pins are fragile (see the
conda-transaction warning below).

## 2. Set the run options (one time per file)

For each file you run: **Run → Configuration per file… (Ctrl+F6)** →

- **Working directory** → choose **"The directory of the file being
  executed"** (or set it explicitly to
  `D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal`).
- Put any arguments in the **Command line options** box (table below).

## 3. The files and what to type

All paths relative to `D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\`.

| File | Command line options | Runtime | What you get |
|---|---|---|---|
| `_diag_stage3.py` | *(none)* | ~3 min | Reproduces the **ground-walk winner** from `curriculum_stage3.json` (2026-09-20: trial 52, score −85.2): prints duty/knee/tilt metrics; updates `spinal_run.npz`. Eval mode does NOT write the png. |
| `_run_stage2_air.py` | *(none, or `16`)* | ~5 min | **Air-walk winner** run; writes `spinal_run.png` + npz and prints joint ranges / RG activity. |
| `_onset_smoke.py` | *(none)* | ~5 min | The two `contact_onset` regression gates (gain 0 bit-identity; gain 0.6 finite). |
| `_debug_rhythm_tuned.py` | *(none, or `2.0 20` = drive, seconds)* | ~1 min | Network-only constant-drive self-sustain probe of the tuned config. |
| `_curriculum.py` | **required**: `3 15` (stage, n_trials) | ~1.5 min **per trial** | Resumes optuna studies `curr_s2b_air_aff` / `curr_s3b_ground`. **Overwrites** `curriculum_stage2/3.json` with the new best (winners also live in the db). For long runs, prefer `run_curriculum_20260920.bat` in a terminal. |
| `runner.py` (3D view) | `--fitted --best10 --view --time 22` | live | Real-time 3D view of the **v10** winner (`--view` opens the MuJoCo window; `--realtime` = 1x playback; `--scope` = live neural traces). The curriculum ground winner has no `--view` loader yet — use `_diag_stage3.py` for it. |
| `check_rhythm.py` | *(none, or `2.0 12`)* | ~1 min | Rhythm-layer calibration with DEFAULT params (no fitted baseline) — diagnostic only, a tonic result there is not a regression. |

## 4. Terminal alternative (long runs)

```bat
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" _diag_stage3.py
```

The full curriculum chain (stage 2 then 3) is `run_curriculum_20260920.bat`
in the same folder — run it from a terminal, not Spyder, because it takes
~1.5 h.

## 5. Cautions

- **Never pip/conda anything else into `myo`** without checking
  `python.exe` still exists and the NumPy/SciPy/MuJoCo pins survive
  (a conda graphviz install once clobbered `python.exe`; see the
  `python` skill / AGENTS.md).
- Eval-mode runs (`_diag_stage3.py`, curriculum trials) refresh
  `spinal_run.npz` but not `spinal_run.png` — only non-eval runs
  (`_run_stage2_air.py`, `runner.py` without `--eval`) write the png.
- `spinal_run.npz` `q` channels are the 13 named KEY_JOINTS **already in
  degrees** (`qfull` is the raw qpos) — do not apply `np.degrees` again.
- `_curriculum.py` with no arguments defaults to stage 1, 30 trials —
  always pass the stage explicitly.
