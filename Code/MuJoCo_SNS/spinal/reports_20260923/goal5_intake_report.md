# Gait-library intake report (2026-09-23)

- OK Case_40_motion: T_r 1.113 s, duty 0.57/0.57, knee_min -59.2 deg, hip range 51.6 deg

References saved: 1 -> D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\gait_refs
Use: load a .npz and rebuild the dict (arrays r_/l_ + scalar keys) as kine_ref.REF_CACHE for scoring/training against it (see gait_lib_score.py / gait_lib_pilot2.py).

muscfib fiber-length files (no GRF columns) will list as UNPAIRED - they integrate through the F-L-V validation route, not the reference route.
