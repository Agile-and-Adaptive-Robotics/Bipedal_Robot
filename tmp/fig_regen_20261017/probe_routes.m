function probe_routes()
% Probe for the figure-regen task: confirm the dated cam result mat holds
% the needed variables, and print the route contact-elimination transitions
% for the run-8 design, the run-8b final design, and the stock baseline,
% so the regenerated route figures can pick informative poses.
root = 'D:/Github/Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

moDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization');
matFile = fullfile(moDir, 'Results', 'Vas_Pam_20mm_Result_cam_20261014_0416.mat');

fprintf('--- whos -file %s ---\n', matFile);
w = whos('-file', matFile);
for k = 1:numel(w)
    fprintf('  %-22s %-10s %s\n', w(k).name, w(k).class, mat2str(w(k).size));
end

S = load(matFile);
ctx = S.ctx;
fprintf('\nphiD: %+.1f .. %+.1f deg (%d poses)\n', ...
    min(ctx.phiD), max(ctx.phiD), numel(ctx.phiD));
fprintf('KMAX = %g, requiredTorqueMargin = %g\n', ctx.KMAX, ctx.requiredTorqueMargin);
fprintf('targetName = %s\n', ctx.targetName);

reportTrans('run8b FINAL (from mat)', S.predBest.routeInfo, ctx.phiD);

% run-8 winner geometry verbatim from tmp/optcam_run8_20261013.log
% lines 2236-2239 (p6/p7/p8/tendon).
xCamRun8 = [0.116323, 0.015040, ...   % p6
            0.131323, -0.018974, ...  % p7
            0.066359, -0.018981, ...  % p8
            0.090440];                % tendon
[xFull8, ctxCam8] = camToExtDesign20mm(xCamRun8, ctx, S.xFixed);
pred8 = predictKneeExt20mm(xFull8, ctxCam8);
reportTrans('run8 (recomputed)', pred8.routeInfo, ctx.phiD);

reportTrans('stock baseline', S.routeInfoBase, ctx.phiD);

% candidateRaw sanity: rows, and whether row 8 differs from row 7 (the
% "p8 repeats p7" note applies to the routed points, not CAD candidates).
ri = S.predBest.routeInfo;
fprintf('\nrun8b candidateRaw size %s, raw size %s\n', ...
    mat2str(size(ri.candidateRaw)), mat2str(size(ri.raw)));
fprintf('candidate p7 vs p8 at pose 1: p7=[%.4f %.4f] p8=[%.4f %.4f]\n', ...
    ri.candidateRaw(7,1,1), ri.candidateRaw(7,2,1), ...
    ri.candidateRaw(8,1,1), ri.candidateRaw(8,2,1));
end

function reportTrans(name, ri, phiD)
trans = find(any(ri.active(:,1:end-1) & ~ri.active(:,2:end), 1));
fprintf('\n%s: %d transition poses\n', name, numel(trans));
for ii = trans(:).'
    removed = find(ri.active(:,ii) & ~ri.active(:,ii+1));
    fprintf('  pose %4d: theta=%+7.1f deg, removes: p%s\n', ii, phiD(ii), ...
        strjoin(compose('p%d', removed), ' '));
end
a1 = find(ri.active(:,1)); aE = find(ri.active(:,end));
fprintf('  active @pose1: p%s ; active @end: p%s\n', ...
    strjoin(compose('p%d', a1), ' '), strjoin(compose('p%d', aE), ' '));
end
