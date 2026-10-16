"""Clean up the spinal folder: move one-off diagnostic/probe/audit scripts to _archive/."""
import shutil
from pathlib import Path

SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
ARCHIVE = SPINAL / "_archive"
ARCHIVE.mkdir(exist_ok=True)

# Patterns for one-off scripts to archive
ARCHIVE_PATTERNS = [
    # Diagnostics and probes (one-off debugging)
    "_diag_*.py", "_probe_*.py", "_audit_*.py", "_check_*.py",
    # Literature fetching (one-off)
    "_zot_*.py", "_dirusso_pdf*.py", "_shev_*.py", "_lit_fetch*.py",
    "_fetch_*.py", "_mccrea_*.py",
    # Superseded extraction iterations
    "_w2l_extract*.py", "_onset_smoke*.py",
    # One-off figure generation (superseded by _figs_20261003.py)
    "_fig_*.py", "_figure_*.py",
    # One-off stage runs and repros (superseded by _replay_best/_resume)
    "_repro_*.py", "_run_t*.py", "_stage2_*.py", "_stage_remaining.py",
    "_update_stage3*.py",
    # Old curriculum resets
    "_curriculum_reset*.py",
    # Old experiment scripts
    "_pose_*.py", "_nap_test.py", "_nap_opt.py", "_airwalk_*.py",
    "_air_capture_*.py", "_bpa_win_*.py", "_crop_*.py",
    "_gw_test.py", "_jpf_*.py", "_kine_v2_test.py",
    "_s3d_probe.py", "_s3f_probe.py", "_s3g_probe.py", "_s3j_probe.py",
    "_s3k_probe.py", "_seed_s3i.py",
    "_smoke_nap.py", "_sto_diag.py",
    "_test_max_net_numpy.py", "_test_max_realistic.py",
    # Misc one-offs
    "_activation_provenance_audit.py", "_aarl_enum.py",
    "_deaff_matrix.py", "_emit_group_table.py", "_eval_once.py",
    "_fsa_rank_robustness.py", "_gen_joint_pf_weights.py",
    "_jointlayer_waveforms.py", "_kine_check_final.py",
    "_live_plant_audit.py", "_model_geom_check.py",
    "_normal_so_audit.py", "_npz_layout_check.py",
    "_overlay_refresh.py", "_panels_check.py", "_patch_fig.py",
    "_pdf_text.py", "_peek_progress.py", "_pf_basis_audit.py",
    "_pf_hypo_probe.py", "_phase_query.py", "_purge_s3g.py",
    "_reg_fingerprint.py", "_render_diagram.py", "_render_gif.py",
    "_render_panels.py", "_rg_collapse_check.py", "_rom_audit.py",
    "_run_nocross.py", "_synergy_count_audit.py", "_synergy_nmf.py",
    "_synergy_pf_study.py", "_tau_h_check.py",
    "_validate_subject_model.py", "_w2l_to_connectome.py",
]

moved = 0
skipped = 0
for pat in ARCHIVE_PATTERNS:
    for f in SPINAL.glob(pat):
        if f.is_file() and f.suffix == ".py":
            dest = ARCHIVE / f.name
            if dest.exists():
                skipped += 1
            else:
                shutil.move(str(f), str(dest))
                moved += 1

print(f"Moved {moved} scripts to _archive/; {skipped} already existed there.")
remaining = [f.name for f in SPINAL.glob("_*.py")]
print(f"\nRemaining _-prefixed scripts ({len(remaining)}):")
for f in sorted(remaining):
    print(f"  {f}")
