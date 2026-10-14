%% Numeric probe of the adj-route overlay in the display route figure:
%% per tile, where do the converted p6/p7/p8 markers actually land,
%% and do repeats stack on top of each other?
root = 'D:/GitHub/Bipedal_Robot';
moDir = fullfile(root,'Code','Matlab','Mesh_Optimization');
addpath(genpath(fullfile(root,'Code','Matlab')));
addpath(moDir);
addpath(fullfile(root,'Testing_Data','2022_02_Festo'),'-end');

ctx = buildKneeExtContext20mm();
S = load(fullfile(moDir,'Results','Vas_Pam_20mm_Result.mat'),'xBest');
xBest = reshape(S.xBest,1,[]);
p1 = xBest(1:3); pEnd = xBest(4:6); tendon = xBest(8);

[~,~,routeInfo] = buildDistalRingLocation20mm(p1,pEnd,tendon,ctx);

% adj block (run-8 winner) verbatim
tendonAdj = 0.090440;
p6Adj = [0.116323, 0.015040];
p7Adj = [0.131323, -0.018974];
p8Adj = [0.066359, -0.018981];
ctxAdj = ctx;
ctxAdj.geo.seedP6XY = p6Adj;
ctxAdj.geo.seedP7XY = p7Adj;
ctxAdj.geo.seedP8XY = p8Adj;
ctxAdj.geo.tibiaWrapCircumcircle = true;
[~,~,routeInfoAdj] = buildDistalRingLocation20mm(p1,pEnd,tendonAdj,ctxAdj);

transitionIdx = find(any(routeInfo.active(:,1:end-1) & ~routeInfo.active(:,2:end), 1));
plotIdx = unique([1, transitionIdx, numel(ctx.phiD)], 'stable');

fprintf('tiles at deg: %s\n', mat2str(round(ctx.phiD(plotIdx),1)));
for t = 1:numel(plotIdx)
    ii = plotIdx(t);
    PAdj = routeInfoAdj.raw(:,:,ii);
    for j = 6:9
        q = RowVecTrans(ctx.T_ICR_t1(:,:,ii), PAdj(j,:));
        PAdj(j,:) = RowVecTrans(ctx.T_Pam(:,:,ii), q);
    end
    act6 = routeInfoAdj.active(6,ii);
    act8 = routeInfoAdj.active(8,ii);
    fprintf('tile %4.1f deg: adj p6 [%7.1f %7.1f] active=%d | p7 [%7.1f %7.1f] | p8 [%7.1f %7.1f] active=%d | p8==p7? %d | p6==p7? %d\n', ...
        ctx.phiD(ii), 1000*PAdj(6,1:2), act6, 1000*PAdj(7,1:2), ...
        1000*PAdj(8,1:2), act8, isequal(PAdj(8,:),PAdj(7,:)), isequal(PAdj(6,:),PAdj(7,:)));
end
disp('PROBE DONE');
