from pathlib import Path
import zipfile,hashlib
repo=Path(r'D:\GitHub\Bipedal_Robot');dis=repo/'Documentation/Reports and Papers/Dissertation';p=dis/'Notes/sns_rebuild_session1_20261030/structure_revision_20261031';arch=dis/'Archive/session1_structure_20261031';arch.mkdir(exist_ok=True)
base=dis/'Archive/session1_submission_20261031/Overleaf_source.zip'
replacement={'chapters/60-sns.tex':p/'60-sns.tex','chapters/95-AppendixD.tex':p/'95-AppendixD.tex','figs/60-sns/sns_library_clean.pdf':p/'sns_library_clean.pdf','figs/60-sns/knee_transient_current_layout.pdf':dis/'Figures/60-sns/knee_transient_current_layout.pdf'}
with zipfile.ZipFile(base) as src,zipfile.ZipFile(arch/'Overleaf_source_updated.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
 names=set(src.namelist())
 for n in src.namelist(): z.writestr(n,replacement[n].read_bytes() if n in replacement else src.read(n))
 for n,q in replacement.items():
  if n not in names:z.writestr(n,q.read_bytes())
report='''# Session 1: chapter opening, figures, and final correction

Machine: easteregg2.

Live project: https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb

Completed: merged the duplicate opening into 6.1 Motivation: The AnimatLab Story (about 220 words including the caption); 6.2 Knowledge Gap and Framing is about 150 words; 6.3 is Gait2392 and Controller Transfer. The older controller remains a standalone sideways figure. The short AnimatLab neural/mechanical plot moved to Appendix D; its brief measured result remains in the motivation. The newer treadmill example moved with the later W2L tests. All existing labels are preserved across the chapter and appendix.

Removed all bold Verification pseudo-headings and the MATLAB paragraph heading, blended procedures and results into connected prose, and added explicit references to the existing Rybak, Shevtsova, Deng, Di Russo, and Shinohara diagrams in Chapter 2. Replaced pilot with simulation. Added 6 pt spacing between curriculum and model-table rows.

Replaced the dated SNS library image with sns_library_clean.pdf: twelve plain, legible functional blocks, no decorative pennation/mesh artwork, and a readable vector neuron signal-flow diagram checked against Code/Matlab/SNS_Simscape/sns_build_library.m. No model equations changed.

Final urgent correction: knee_transient_current_layout.pdf panel C now has FLEXOR pathways on the LEFT and EXTENSOR pathways on the RIGHT, preserving the connections, signs, and plotted results. Visually checked the standalone corrected PDF. Uploaded this replacement to figs/60-sns and compiled. Final live compiler: Errors 0 / Warnings 0. The preceding structural build is 311 pages. Repeated final browser PDF downloads timed out; use the live Overleaf download for the final orientation-corrected PDF rather than assuming an older local PDF includes that last correction.

Speed table now includes measured 100.17x native-SNS neural acceleration (24.533797 to 0.244924 s for 500 updates, 1027 neurons), quiet repeat 56.55x (14.306414 to 0.253001 s), and Python-to-Numba Di Russo equation-kernel acceleration 584.33x (1977.854 to 3.38485 microseconds per update), alongside 9.21x/17.75x whole coupled-loop results. No exact 200x/240x record was recovered. Neural-only measurements exclude physics and per-step conductance updates. Sources: spinal/hunt_training/runs/fourcell_hunt_20261029_v1/acceleration_check.json; fourcell_dynamic_20261029_v1/acceleration_check.json; spinal/reports_20261028_publication_flow/dirusso_equations/dirusso_neural_benchmark.json.

Local final sources and new vector figure are in this folder. The active corrected knee figure is at Dissertation/Figures/60-sns/knee_transient_current_layout.pdf. Full source package reconstructed from the authoritative downloaded live ZIP plus exactly these four uploaded replacements: Dissertation/Archive/session1_structure_20261031/Overleaf_source_updated.zip. The package is a local reconstruction, not a newly downloaded final ZIP. No published muscle or torque chapter changes. No new simulations or controller training.
'''
(p/'SESSION1_STRUCTURE_HANDOFF.md').write_text(report,encoding='utf-8')
for q in [repo/'CHATGPT_REPORT.md',repo/'SESSION1_HANDOFF_20261031.md']:
 with q.open('a',encoding='utf-8') as f:f.write('\n\n## Session 1 — urgent chapter opening and final knee orientation correction\n\nLive Ch. 6 now has 6.1 compact AnimatLab motivation, 6.2 compact gap/framing, and 6.3 Gait2392. Bold Verification labels removed; Background circuit figures cross-referenced; clean vector SNS library diagram replaces distorted icons/tiny wiring; curriculum rows spaced. Neural-only 100.17x/56.55x and Python-to-Numba 584.33x timings added with explicit scopes. Final knee panel C: flexor LEFT, extensor RIGHT. Live compile Errors 0 / Warnings 0. Use live download for final PDF; last browser download attempt timed out. Full handoff: '+str(p/'SESSION1_STRUCTURE_HANDOFF.md')+'.\n')
print(p/'SESSION1_STRUCTURE_HANDOFF.md')
