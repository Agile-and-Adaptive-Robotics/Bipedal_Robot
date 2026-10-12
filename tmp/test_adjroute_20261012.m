%% Smoke harness: clip-part adj-route override plumbing (2026-10-12)
% Exercises the geo.seedP6XY/seedP7XY/seedP8XY overrides in
% buildDistalRingLocation20mm with NON-default values (p6 rotated down
% in theta about the t1 origin, longer tendon) and reports what moves.

root = 'D:/GitHub/Bipedal_Robot';
scriptDir = fullfile(root,'Code','Matlab','Mesh_Optimization','Results');
addpath(genpath(fullfile(root,'Code','Matlab')));
addpath(fullfile(root,'Code','Matlab','Mesh_Optimization'));
addpath(fullfile(root,'Testing_Data','2022_02_Festo'),'-end');

S = load(fullfile(scriptDir,'Vas_Pam_20mm_Result.mat'),'xBest','XiUsed','ctx');
if isfield(S,'XiUsed')
    xiRun = S.XiUsed;
elseif isfield(S,'ctx') && isfield(S.ctx,'Xi0')
    xiRun = [S.ctx.Xi0,S.ctx.Xi1,S.ctx.Xi2,S.ctx.Xi3];
else
    error('no Xi record in result mat')
end

ctx = buildKneeExtContext20mm();
xiBuilt = [ctx.Xi0,ctx.Xi1,ctx.Xi2,ctx.Xi3];
assert(isequal(xiBuilt,xiRun),'Xi mismatch');

xBest = reshape(S.xBest,1,[]);
p1 = xBest(1:3);
pEnd = xBest(4:6);
tendon = xBest(8);

[LocationBase,~,routeInfoBase] = buildDistalRingLocation20mm(p1,pEnd,tendon,ctx);

% adj: p6 rotated DOWN in theta (theta' < theta) at constant radius
% about the t1 frame origin; p7/p8 unchanged; tendon lengthened.
p6Base = routeInfoBase.seed20(6,1:2);
th0 = atan2(p6Base(2),p6Base(1));
r0 = hypot(p6Base(1),p6Base(2));
th1 = th0 - deg2rad(8);
p6Adj = r0*[cos(th1),sin(th1)];
tendonAdj = tendon + 0.004;

ctxAdj = ctx;
ctxAdj.geo.seedP6XY = p6Adj;
ctxAdj.geo.seedP7XY = routeInfoBase.seed20(7,1:2);
ctxAdj.geo.seedP8XY = routeInfoBase.seed20(8,1:2);

[LocationAdj,~,routeInfoAdj] = buildDistalRingLocation20mm(p1,pEnd,tendonAdj,ctxAdj);

fprintf('theta(p6): %.3f -> %.3f deg\n', ...
    atan2d(p6Base(2),p6Base(1)), atan2d(p6Adj(2),p6Adj(1)));
fprintf('adj p6 override honored: %d\n', ...
    isequal(routeInfoAdj.seed20(6,1:2),p6Adj));
fprintf('p7/p8 override honored: %d\n', ...
    isequal(routeInfoAdj.seed20(7,1:2),routeInfoBase.seed20(7,1:2)) && ...
    isequal(routeInfoAdj.seed20(8,1:2),routeInfoBase.seed20(8,1:2)));
fprintf('max |dLocation| vs baseline: %.6g\n', ...
    max(abs(LocationAdj(:)-LocationBase(:))));
if ~isequal(routeInfoAdj.seed20(3:5,1:2),routeInfoBase.seed20(3:5,1:2))
    fprintf('p3:p5 seed rows MOVED (femur-chain repair fired):\n');
    for r = 3:5
        fprintf('  p%d: [%.5f %.5f] -> [%.5f %.5f] m\n', r, ...
            routeInfoBase.seed20(r,1:2), routeInfoAdj.seed20(r,1:2));
    end
else
    fprintf('p3:p5 seed rows unchanged (repair did not fire)\n');
end
fprintf('elimination angles, baseline -> adj (p2:p8):\n');
for r = 2:8
    fprintf('  p%d: %+7.2f -> %+7.2f deg\n', r, ...
        routeInfoBase.eliminatedAngleD(r), routeInfoAdj.eliminatedAngleD(r));
end

% Baseline bit-identity guard: overrides ABSENT must reproduce the
% pre-edit builder output exactly (LocationBase was just rebuilt with
% the edited builder; compare its release schedule to the hardcoded
% record from Knee_Extensor_20mm.m).
hardcodedD = [-89.797979628815; -55.656565451176; -47.777777564029; ...
              -20.202019959013; 6.060606331479; 7.373737646003];
rowsHard = [5; 4; 6; 7; 3; 8];
for k = 1:numel(rowsHard)
    row = rowsHard(k);
    trans = find(routeInfoBase.active(row,1:end-1) ...
        & ~routeInfoBase.active(row,2:end), 1);
    if isempty(trans)
        derivedD = NaN;
    else
        derivedD = ctx.phiD(trans+1);
    end
    fprintf('p%d release: hardcoded %+.4f, derived %+.4f\n', ...
        row, hardcodedD(k), derivedD);
    assert(abs(derivedD-hardcodedD(k)) < 1e-9, 'baseline drift');
end

disp('SMOKE OK');
