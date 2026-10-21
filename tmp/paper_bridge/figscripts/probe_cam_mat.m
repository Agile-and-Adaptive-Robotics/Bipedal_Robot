% probe_cam_mat.m — dump the structure of the cam result mat of record
% (read-only; scratch per spec: figscripts only).
matFile = 'D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result_cam_20261014_0416.mat';
S = load(matFile);
fprintf('TOP-LEVEL vars:\n');
disp(fieldnames(S).');

ctx = S.ctxCamBest;
fprintf('ctxCamBest fields:\n');
disp(fieldnames(ctx).');
fprintf('phiD: %d angles, [%g .. %g] deg\n', numel(ctx.phiD), min(ctx.phiD), max(ctx.phiD));
fprintf('N = %d\n', ctx.N);

fprintf('geo fields:\n');
g = ctx.geo;
fn = fieldnames(g);
for k = 1:numel(fn)
    v = g.(fn{k});
    if isnumeric(v) && numel(v) <= 4
        fprintf('  %-28s %s  %s\n', fn{k}, mat2str(size(v)), strjoin(string(v), ','));
    else
        fprintf('  %-28s %s\n', fn{k}, mat2str(size(v)));
    end
end

ri = S.predBest.routeInfo;
fprintf('routeInfo fields:\n');
disp(fieldnames(ri).');
fprintf('raw %s  active %s\n', mat2str(size(ri.raw)), mat2str(size(ri.active)));
fprintf('eliminatedAngleD (rows 1-9):\n');
disp(ri.eliminatedAngleD(:).');
fprintf('active at full flexion (pose 1):\n');
disp(ri.active(:,1).');
fprintf('active at full extension (pose end):\n');
disp(ri.active(:,end).');
fprintf('active per pose (row = p1..p9, col = pose):\n');
disp(ri.active);

fprintf('predBest fields:\n');
disp(fieldnames(S.predBest).');
fprintf('seed20 (native frames):\n');
disp(ri.seed20);
fprintf('tendonMax range: %.4f .. %.4f m\n', min(ri.tendonMax), max(ri.tendonMax));
fprintf('wrap fields:\n');
disp(fieldnames(ri.wrap).');
fprintf('wrap.radius size %s, wrap angle deg row2 (tibia) range %.1f..%.1f\n', ...
    mat2str(size(ri.wrap.radius)), min(ri.wrap.angleDeg(2,:)), max(ri.wrap.angleDeg(2,:)));
disp('DONE_probe');
