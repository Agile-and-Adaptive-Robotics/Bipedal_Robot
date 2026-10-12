%% Reproduce the requiredPressureScanCam line-338 || error with the
%% real winner values, printing shapes along the way.
root = 'D:/GitHub/Bipedal_Robot';
moDir = fullfile(root,'Code','Matlab','Mesh_Optimization');
addpath(genpath(fullfile(root,'Code','Matlab')));
addpath(moDir);
addpath(fullfile(root,'Testing_Data','2022_02_Festo'),'-end');

ctx = buildKneeExtContext20mm();
S = load(fullfile(moDir,'Results','Vas_Pam_20mm_Result.mat'),'xBest');
xRef8 = reshape(S.xBest,1,[]);
xFixed = xRef8(1:7);

xCamW = [0.100710, 0.021345, 0.131610, -0.014208, 0.068434, -0.016017, 0.088381];
[xFull,ctxCam] = camToExtDesign20mm(xCamW, ctx, xFixed);
ctxEval = ctxCam;

Pgrid = (0:20:620)';
humanTgt = interp1(ctxEval.humanAngleD, ctxEval.humanTorqueAbs, ...
    ctxEval.phiD, 'pchip', 'extrap');
hasTgt = isfinite(humanTgt) & humanTgt > 0;

fprintf('shapes: humanAngleD %s, humanTorqueAbs %s, phiD %s, humanTgt %s, hasTgt %s\n', ...
    mat2str(size(ctxEval.humanAngleD)), mat2str(size(ctxEval.humanTorqueAbs)), ...
    mat2str(size(ctxEval.phiD)), mat2str(size(humanTgt)), mat2str(size(hasTgt)));
fprintf('any hasTgt: %d, class: %s\n', any(hasTgt), class(hasTgt));

restV = xFull(7);
kmaxV = restV*(1 - ctxEval.KMAX);
[LocationV, bendV] = buildDistalRingLocation20mm( ...
    xFull(1:3), xFull(4:6), xFull(8), ctxEval);

Tq = nan(numel(Pgrid), numel(ctxEval.phiD));
for k = 1:numel(Pgrid)
    b = MonoPamDataExplicit_balanceX3(ctxEval.Name, LocationV, ...
        ctxEval.CrossPoint, ctxEval.Dia, ctxEval.T_Pam, restV, kmaxV, ...
        xFull(8), ctxEval.fitting, Pgrid(k), ctxEval.Xi0, ctxEval.Xi1, ...
        ctxEval.Xi2, ctxEval.Xi3, ctxEval.wraps, ctxEval.phiD, ...
        ctxEval.BPAcount, bendV);
    Tq(k,:) = abs(b.Torque_p(:,3)).';
end
fprintf('Tq shape: %s, NaN count: %d\n', mat2str(size(Tq)), sum(isnan(Tq(:))));

% The exact failing loop, with shape probes on the first iteration.
iFirst = find(hasTgt, 1);
col = Tq(:,iFirst);
fin = isfinite(col);
fprintf('first iter: i=%d, col %s, fin %s, any(fin) %s, maxCol %s, tgt(i) %s class %s\n', ...
    iFirst, mat2str(size(col)), mat2str(size(fin)), mat2str(any(fin)), ...
    mat2str(size(max(col(fin)))), mat2str(size(humanTgt(iFirst))), class(humanTgt(iFirst)));

try
    if ~any(fin) || max(col(fin)) < humanTgt(iFirst)
        fprintf('combined || : true branch\n');
    else
        fprintf('combined || : false branch\n');
    end
catch ME
    fprintf('combined || FAILED: %s\n', ME.message);
end

% Full scan with the fixed element-wise loop.
Preq = nan(numel(ctxEval.phiD),1);
for i = find(hasTgt).'
    col = Tq(:,i);
    fin = isfinite(col);
    if ~any(fin)
        continue
    end
    if max(col(fin)) < humanTgt(i)
        continue
    end
    kk = find(col >= humanTgt(i), 1);
    if kk <= 1
        Preq(i) = 0;
    else
        pq = Pgrid([kk-1 kk]);
        [xs, isrt] = sort(col([kk-1 kk]));
        Preq(i) = interp1(xs, pq(isrt), humanTgt(i));
    end
end
valid = find(hasTgt);
fprintf('separate-if scan: min %.1f kPa, max %.1f kPa over %d target angles\n', ...
    min(Preq(valid),[],'omitnan'), max(Preq(valid),[],'omitnan'), numel(valid));

disp('DEBUG DONE');
