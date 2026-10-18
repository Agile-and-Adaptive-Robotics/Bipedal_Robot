function probe_window()
% Print femur-frame coordinates of route points at the three zoom poses
% (1, 50, 100) for run8b + candidates, so the distal zoom window can be
% sized to the contact cluster instead of the full leg span.
root = 'D:/Github/Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
moDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization');
S = load(fullfile(moDir, 'Results', 'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), ...
    'ctx', 'predBest', 'routeInfoBase');
ctx = S.ctx; phiD = ctx.phiD(:);
ri = S.predBest.routeInfo;

for ii = [1, 50, 100]
    fprintf('\n=== pose %d, theta %+.1f deg ===\n', ii, phiD(ii));
    for src = 1:3
        if src == 1, tag = 'cam  '; Q0 = ri.raw(:,:,ii);
        elseif src == 2, tag = 'cand'; Q0 = ri.candidateRaw(:,:,ii);
        else, tag = 'base '; Q0 = S.routeInfoBase.raw(:,:,ii);
        end
        Q = Q0;
        for j = 6:9
            qICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), Q(j,:));
            Q(j,:) = RowVecTrans(ctx.T_Pam(:,:,ii), qICR);
        end
        if src == 1, act = ri.active(:,ii); else, act = nan(9,1); end
        fprintf('%s rows 3:9 femur-frame (active 0/1 cam only):\n', tag);
        for j = 3:9
            fprintf('  p%d [%+7.4f %+7.4f]  act=%g\n', j, Q(j,1), Q(j,2), act(j));
        end
    end
end

% Bounding boxes for candidate windows
fprintf('\n--- candidate windows (over poses 1/50/100) ---\n');
allP = []; allP569 = []; allP69 = [];
for ii = [1, 50, 100]
    for src = [ri.raw(:,:,ii), ri.candidateRaw(:,:,ii)]
    end
end
for ii = [1, 50, 100]
    for Q0 = {ri.raw(:,:,ii), ri.candidateRaw(:,:,ii), S.routeInfoBase.raw(:,:,ii)}
        Q = Q0{1};
        for j = 6:9
            qICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), Q(j,:));
            Q(j,:) = RowVecTrans(ctx.T_Pam(:,:,ii), qICR);
        end
        allP   = [allP;   Q(3:9,1:2)]; %#ok<AGROW>
        allP569= [allP569;Q(5:9,1:2)]; %#ok<AGROW>
        allP69 = [allP69; Q(6:9,1:2)]; %#ok<AGROW>
    end
end
boxes = {allP, allP569, allP69};
names = {'rows 3:9', 'rows 5:9', 'rows 6:9'};
for k = 1:3
    P = boxes{k};
    fprintf('%s: x [%+.4f %+.4f] (%.3f m) y [%+.4f %+.4f] (%.3f m)\n', ...
        names{k}, min(P(:,1)), max(P(:,1)), max(P(:,1))-min(P(:,1)), ...
        min(P(:,2)), max(P(:,2)), max(P(:,2))-min(P(:,2)));
end
% tibia clearance circle centers/radii (femur frame, pose 1) for context
fprintf('\ntibiaLowerCenter t1 [%+.4f %+.4f] R %.4f\n', ...
    ctx.geo.tibiaLowerCenter(1), ctx.geo.tibiaLowerCenter(2), ...
    ctx.geo.tibiaLowerClearRadius);
fprintf('tibiaUpperCenter t1 [%+.4f %+.4f] R %.4f\n', ...
    ctx.geo.tibiaUpperCenter(1), ctx.geo.tibiaUpperCenter(2), ...
    ctx.geo.tibiaUpperClearRadius);
fprintf('femurCylCenter [%+.4f %+.4f] R %.4f\n', ...
    ctx.geo.femurCylCenter(1), ctx.geo.femurCylCenter(2), ...
    ctx.geo.femurCylClearRadius);
end
