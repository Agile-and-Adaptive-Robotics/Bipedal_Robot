%run_ExtCompare_record_20261026.m -- Ben's pipeline ruling (2026-10-26,
% evening): the MEASURED extensor 20 mm torques of record are whatever the
% compare section of Knee_Extensor_20mm.m computes, run exactly as it
% stands (script READ-ONLY; this dated runner only runs it and saves the
% resulting workspace variables). The xlsx rows 16/17 are sanity checks
% only. run() executes the script in this runner's workspace; the script's
% leading clear/restore is expected and the variables it creates (TorqueZ,
% predAP, ...) survive after run() returns, when they are saved to a dated
% mat and printed.
scriptPath = 'D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results/Knee_Extensor_20mm.m';
outMat = 'D:/GitHub/Bipedal_Robot/tmp/paper_bridge/figscripts/ext20_compare_record_20261026.mat';
run(scriptPath);
fprintf('\n=== RECORD (compare section output, verbatim variables) ===\n');
fprintf('K_deg  : '); fprintf('%9.1f', K_deg); fprintf('\n');
fprintf('Load N : '); fprintf('%9.1f', Load); fprintf('\n');
fprintf('LC tab : '); fprintf('%9.1f', -LC_deg); fprintf('\n');
fprintf('Presh  : '); fprintf('%9.1f', Presh); fprintf('\n');
fprintf('TorqueZ: '); fprintf('%9.4f', TorqueZ); fprintf('\n');
fprintf('predAP : '); fprintf('%9.4f', predAP); fprintf('\n');
fprintf('human@K: '); fprintf('%9.4f', interp1(ctx.humanAngleD(:), ctx.humanTorque(:), K_deg, 'pchip')); fprintf('\n');
fprintf('pred620@K: '); fprintf('%9.4f', interp1(phiD, Torque2_adj, K_deg, 'pchip')); fprintf('\n');
save(outMat, 'K_deg', 'Load', 'LC_deg', 'Presh', 'TorqueZ', 'predAP', ...
    'Torque2_adj', 'phiD', 'ctx', 'Name', 'LocationAdj', 'bendMeasureAdj', ...
    'tendonAdj', 'CrossPoint', 'Dia', 'T_Pam', 'rest', 'kmax', 'fitting', ...
    'Xi0', 'Xi1', 'Xi2', 'Xi3', 'wraps', 'BPAcount', '-v7.3');
fprintf('SAVED %s\n', outMat);
disp('DONE_run_ExtCompare_record');
