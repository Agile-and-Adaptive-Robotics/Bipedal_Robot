@echo off
rem Post-curriculum deliverables (laminated, 2026-09-16): stage-3 winner
rem full 22 s run + figures + diagrams + Dissertation copies.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
set LOG=post_curriculum_20260916.log
set DISS="D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\CPG_airstepping_figs"

echo === FINAL RUN START %date% %time% === >> %LOG%
%PY% _final_run_lam.py >> %LOG% 2>&1 || goto :fail
%PY% _final_metrics.py ground_curriculum_final.npz >> %LOG% 2>&1
echo === RENDER GIF === >> %LOG%
%PY% render_video.py ground_v11_curriculum.gif >> %LOG% 2>&1 || goto :fail
echo === PLOT RUN === >> %LOG%
%PY% plot_run.py >> %LOG% 2>&1
echo === HINDLIMB === >> %LOG%
%PY% _figure_hindlimb_style.py ground_curriculum_final.npz v11_curriculum >> %LOG% 2>&1 || goto :fail
echo === OVERLAY === >> %LOG%
%PY% _overlay_refresh.py ground_curriculum_final.npz >> %LOG% 2>&1 || goto :fail
echo === DENG FIGURE (winner gains, pdf svg png) === >> %LOG%
%PY% draw_circuit.py --which deng --vclasses --fmt pdf,svg,png >> %LOG% 2>&1 || goto :fail
echo === PANELS DIAGRAM === >> %LOG%
%PY% _render_panels.py >> %LOG% 2>&1 || goto :fail
echo === COPY TO DISSERTATION === >> %LOG%
copy /Y figures\circuit_dengstyle.pdf %DISS% >> %LOG%
copy /Y figures\circuit_dengstyle.svg %DISS% >> %LOG%
copy /Y figures\circuit_dengstyle.png %DISS% >> %LOG%
copy /Y figures\sns_diagram_panels.png %DISS% >> %LOG%
copy /Y ground_v11_curriculum.gif %DISS% >> %LOG%
copy /Y figures\hindlimb_style_v11_curriculum.png %DISS% >> %LOG%
copy /Y figures\opensim_overlay_gait_cycles.png %DISS% >> %LOG%
copy /Y opensim_overlay_gait_cycles.png %DISS% >> %LOG%
echo === ALL DELIVERABLES DONE %date% %time% === >> %LOG%
exit /b 0

:fail
echo === DELIVERABLES FAILED %date% %time% === >> %LOG%
exit /b 1
