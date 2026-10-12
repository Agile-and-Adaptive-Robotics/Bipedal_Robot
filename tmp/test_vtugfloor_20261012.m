%% VTUG floor harness (2026-10-12): does each design need >= 200 kPa
% everywhere it carries a target? Prints the new cVtug constraint and a
% coarse required-pressure scan for the baseline construction, the
% 25.4 mm push, and Ben's 1000+1500-budget cam winner.

root = 'D:/GitHub/Bipedal_Robot';
moDir = fullfile(root,'Code','Matlab','Mesh_Optimization');
addpath(genpath(fullfile(root,'Code','Matlab')));
addpath(moDir);
addpath(fullfile(root,'Testing_Data','2022_02_Festo'),'-end');

ctx = buildKneeExtContext20mm();

S = load(fullfile(moDir,'Results','Vas_Pam_20mm_Result.mat'),'xBest');
xRef8 = reshape(S.xBest,1,[]);
xFixed = xRef8(1:7);
tendonRef = xRef8(8);

[~,~,routeInfoBase] = buildDistalRingLocation20mm( ...
    xFixed(1:3), xFixed(4:6), tendonRef, ctx);
seedBase = routeInfoBase.seed20;
thetaP6BaseD = atan2d(seedBase(6,2), seedBase(6,1));

adjPush = 0.0254;
refPts = seedBase(6:8,1:2);
refPts(1,:) = refPts(1,:) + adjPush*refPts(1,:)/norm(refPts(1,:));
refPts(2,:) = refPts(2,:) + adjPush*refPts(2,:)/norm(refPts(2,:));

xCamBase   = [reshape(seedBase(6:8,1:2),1,6), tendonRef];
xCamPush   = [reshape(refPts,1,6), tendonRef];
xCamWinner = [0.114964, 0.014015, 0.131610, -0.030990, ...
              0.058361, -0.031022, 0.090453];   % Ben's 1000+1500 run

vtugFloorKPa = 200;

designs = { ...
    'baseline xBest ', xCamBase; ...
    'push first cut', xCamPush; ...
    'cam winner    ', xCamWinner};

for iD = 1:size(designs,1)
    xCamD = designs{iD,2};

    [cD, ~] = nonlconExtCam20mm(xCamD, ctx, xFixed, thetaP6BaseD, vtugFloorKPa);
    [pminD, pmaxD] = reqPressureScan(xCamD, ctx, xFixed);

    fprintf('\n%s\n', designs{iD,1});
    fprintf('  cVtug@%dkPa = %+.6f (<=0 feasible)\n', vtugFloorKPa, cD(end));
    fprintf('  required pressure: min %.0f kPa, max %.0f kPa\n', pminD, pmaxD);
    fprintf('  (max constraint overall %+.6f)\n', max(cD));
end

disp('VTUG HARNESS OK');

function [pmin, pmax] = reqPressureScan(xCamD, ctx, xFixed)
Pgrid = (0:20:620)';
humanTgt = interp1(ctx.humanAngleD, ctx.humanTorqueAbs, ctx.phiD, ...
    'pchip', 'extrap');
hasTgt = isfinite(humanTgt) & humanTgt > 0;

[xFull, ctxCam] = camToExtDesign20mm(xCamD, ctx, xFixed);
restV = xFull(7);
kmaxV = restV*(1 - ctxCam.KMAX);
[LocationV, bendV] = buildDistalRingLocation20mm( ...
    xFull(1:3), xFull(4:6), xFull(8), ctxCam);

Tq = nan(numel(Pgrid), numel(ctx.phiD));
for k = 1:numel(Pgrid)
    b = MonoPamDataExplicit_balanceX3(ctxCam.Name, LocationV, ...
        ctxCam.CrossPoint, ctxCam.Dia, ctxCam.T_Pam, restV, kmaxV, ...
        xFull(8), ctxCam.fitting, Pgrid(k), ctxCam.Xi0, ctxCam.Xi1, ...
        ctxCam.Xi2, ctxCam.Xi3, ctxCam.wraps, ctxCam.phiD, ...
        ctxCam.BPAcount, bendV);
    Tq(k,:) = abs(b.Torque_p(:,3)).';
end

Preq = nan(numel(ctx.phiD),1);
for i = find(hasTgt)
    col = Tq(:,i);
    fin = isfinite(col);
    if ~any(fin) || max(col(fin)) < humanTgt(i)
        continue
    end
    k = find(col >= humanTgt(i), 1);
    if k <= 1
        Preq(i) = 0;
    else
        pq = Pgrid([k-1 k]);
        [xs, isrt] = sort(col([k-1 k]));
        Preq(i) = interp1(xs, pq(isrt), humanTgt(i));
    end
end

pmin = min(Preq(find(hasTgt)),[],'omitnan');
pmax = max(Preq(find(hasTgt)),[],'omitnan');
end
