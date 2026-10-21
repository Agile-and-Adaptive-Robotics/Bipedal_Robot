% Read-only: cam-route pose-schedule + wrap-radius facts for the A4 figure caption.
S = load('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result_cam_20261014_0416.mat');
ctx = S.ctxCamBest; ri = S.predBest.routeInfo; phiD = ctx.phiD(:); N = numel(phiD);
fprintf('phiD: %d poses, [%.1f .. %.1f] deg\n', N, min(phiD), max(phiD));
transitionIdx = find(any(ri.active(:,1:end-1) & ~ri.active(:,2:end), 1));
transAngles = phiD(transitionIdx(:)).';
transRows = arrayfun(@(t) find(ri.active(:,t) & ~ri.active(:,t+1)), transitionIdx);
for k = 1:numel(transitionIdx)
    fprintf('transition %d: p%d releases at %.2f deg\n', k, transRows(k), transAngles(k));
end
cand = unique([1, transitionIdx, N], 'stable');
fprintf('cand poses = %d -> schedule: %s\n', numel(cand), mat2str(phiD(cand).'));
fprintf('cam-drum wrap radius ri.wrap.radius(2,1) = %.4f m = %.1f mm\n', ri.wrap.radius(2,1), 1000*ri.wrap.radius(2,1));
disp('DONE');
