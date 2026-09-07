%% Debug_RouteElim_Ext.m
% READ-ONLY diagnostic for the 20 mm extensor route-elimination state machine
% in buildDistalRingLocation20mm. Does not modify any model file and does not
% run the optimizer. Safe to delete.
%
% What it does:
%   1. Builds the context and evaluates the initial design x0.
%   2. Prints the elimination table, active-route transitions, p7 seed
%      colinearity angle, and the max eliminations in one knee-angle step.
%   3. Reconstructs, for every sweep angle and every relevant optional row,
%      the builder's native-frame rotated-vector margin (femur +90 deg,
%      tibia -90 deg, principal atan2 values; positive margin = gate opens)
%      and cross-checks it against the machine's own records.
%   4. Audits each elimination event physically: was the bypass chord
%      actually collision-free against the clearance geometry when the row
%      was removed?  Also reports rows whose bypass chord is physically
%      clear over angle ranges where the angle gate still keeps them.

clc

%% Repo root and path setup (derived from this file's location)
dbgDir = fileparts(mfilename('fullpath'));
root = dbgDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root)
        error('Could not locate the Bipedal_Robot repo root from %s', dbgDir)
    end
    root = parent;
end

addpath(genpath(fullfile(root, 'Code', 'Matlab')));
% Mesh_Optimization must win any shadowing contest against data subfolders.
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));

td = fullfile(root, 'Testing_Data');
addpath(td);
sub = dir(td);
sub = sub([sub.isdir] & ~ismember({sub.name}, {'.', '..'}));
for k = 1:numel(sub)
    p = fullfile(td, sub(k).name);
    if isfolder(p)
        addpath(p);
    end
end

fprintf('Builder : %s\n', which('buildDistalRingLocation20mm'))
fprintf('Context : %s\n', which('buildKneeExtContext20mm'))
fprintf('BPA cls : %s\n', which('MonoPamDataExplicit_balanceX3'))
fprintf('Colors  : %s\n', which('Colors'))

ctx = buildKneeExtContext20mm();
pred0 = predictKneeExt20mm(ctx.x0, ctx);
ri = pred0.routeInfo;

if ~pred0.ok
    error('Initial prediction failed: %s', pred0.failReason)
end

%% 1. Elimination table and active transitions
fprintf('\n========== ELIMINATION TABLE (initial design x0) ==========\n')
fprintf('%-5s | %-11s | %-14s\n', 'Point', 'Added, deg', 'Eliminated, deg')
for j = 1:9
    if isnan(ri.eliminatedAngleD(j))
        elimText = 'active';
    else
        elimText = sprintf('% .2f', ri.eliminatedAngleD(j));
    end
    fprintf('p%-4d | %11.2f | %-14s\n', j, ri.addedAngleD(j), elimText)
end

fprintf('Elimination order: ')
fprintf('p%d ', ri.eliminationOrder)
fprintf('\n')

if isfield(ri, 'seedColinearAngleD')
    fprintf('\np7 seed colinearity with the pEnd->p8 ray: %.3f deg\n', ...
        ri.seedColinearAngleD)
    fprintf('Guard tolerance geo.seedColinearTolD = %.3f deg\n', ...
        ctx.geo.seedColinearTolD)
    if ri.initiallyEliminated(7)
        fprintf('p7 STARTS ELIMINATED (colinear within tolerance).\n')
    else
        fprintf('p7 starts active.\n')
    end
end

fprintf('\n========== ACTIVE ROUTE TRANSITIONS ==========\n')
prevActive = false(9,1);
for ii = 1:ctx.N
    a = ri.active(:,ii);
    if ii == 1 || any(a ~= prevActive)
        fprintf('%8.2f deg |', ctx.phiD(ii))
        fprintf(' p%d', find(a))
        fprintf('\n')
    end
    prevActive = a;
end

nPerStep = sum(ri.active(:,1:end-1) & ~ri.active(:,2:end), 1);
fprintf('Max route points eliminated in one knee-angle step: %d\n', ...
    max(nPerStep))

%% 2. Reconstructed margins for every relevant optional row at every angle
[mRaw, mWrap, aBigAll, aSmallAll, actMask] = routeMargins(ri, ctx);

% Cross-check the reconstruction against the machine's own record.
fprintf('\n========== RECONSTRUCTION CROSS-CHECK (s>1 indices) ==========\n')
nMismatch = 0;
for ii = 2:ctx.N
    tested = ri.angleCull.tested(:,ii) & isfinite(mRaw(:,ii));
    d = abs(mRaw(tested,ii) - ri.angleCull.marginD(tested,ii));
    if any(d > 1e-6)
        nMismatch = nMismatch + 1;
        if nMismatch <= 5
            fprintf('ii %3d (%6.2f deg): rows %s differ by up to %.2f deg\n', ...
                ii, ctx.phiD(ii), mat2str(find(tested)'), max(d))
        end
    end
end
if nMismatch == 0
    fprintf('Reconstruction matches the machine at every tested index.\n')
end

%% 3. Per-row gate summary with clearance audit
fprintf('\n========== PER-ROW GATE SUMMARY (x0) ==========\n')
fprintf('%-5s | %-12s | %14s | %14s | %14s\n', ...
    'Point', 'Status', 'gate opens at', 'clear at elim?', ...
    'first clear+kept')
fprintf('%s\n', repmat('-',1,74))

for j = 2:8
    if isnan(ri.eliminatedAngleD(j))
        statusText = 'active';
    else
        statusText = sprintf('elim@%7.2f', ri.eliminatedAngleD(j));
    end

    [iOpen, ~, ~] = gateOpenIndices(mRaw(j,:), mWrap(j,:), actMask(j,:));

    % Clearance audit for this row over the sweep.
    [clearMask, evalMask] = rowBypassClearance(ri, ctx, j);

    % Was the bypass chord physically clear at the elimination index?
    elimIdx = ri.eliminatedSweepIndex(j);
    if isfinite(elimIdx) && elimIdx >= 1
        clearAtElim = clearMask(elimIdx);
        elimClearText = 'yes';
        if ~clearAtElim
            elimClearText = 'NO - blocked';
        end
    else
        elimClearText = 'n/a';
    end

    % First angle where the chord is clear but the angle gate still keeps
    % the row ("not eliminated when it could be").
    kept = evalMask & clearMask & ~gateOpens(mRaw(j,:));
    iKept = find(kept, 1, 'first');

    % First angle where the row was removed although the chord was blocked.
    removedBlocked = evalMask & ~clearMask & gateOpens(mRaw(j,:));
    iRB = find(removedBlocked, 1, 'first');

    fprintf('p%-4d | %-12s | %14s | %14s | %14s\n', ...
        j, statusText, angleText(ctx.phiD, iOpen), ...
        elimClearText, angleText(ctx.phiD, iKept))

    if ~isempty(iRB)
        fprintf('       removed while chord BLOCKED first at %s\n', ...
            angleText(ctx.phiD, iRB))
    end
end

% Full-length audit: how much of the sweep has a physically clear bypass
% while the angle gate keeps the point?
fprintf('\n========== CLEAR-BUT-KEPT FRACTION PER ROW (x0) ==========\n')
fprintf('%-5s | %10s | %12s | %12s\n', 'Point', 'eval N', ...
    'clear count', 'clear+kept N')
for j = 2:8
    [clearMask, evalMask] = rowBypassClearance(ri, ctx, j);
    kept = evalMask & clearMask & ~gateOpens(mRaw(j,:));
    fprintf('p%-4d | %10d | %12d | %12d\n', j, ...
        sum(evalMask), sum(evalMask & clearMask), sum(kept))
end

%% 4. Bridge diagnostics: is the straight p1->p9 chord ever physically clear?
fprintf('\n========== P1->P9 BRIDGE CLEARANCE (x0) ==========\n')
bridgeFemur = false(ctx.N,1);
bridgeT1 = false(ctx.N,1);

tol = 1e-8;
for ii = 1:ctx.N
    raw = ri.raw(:,:,ii);

    p9Femur = RowVecTrans(ctx.T_Pam(:,:,ii)*ctx.T_ICR_t1(:,:,ii), raw(9,:));
    bridgeFemur(ii) = femurChordClear(raw(1,1:2), p9Femur(1:2), ctx.geo, tol, 2);

    p1T1 = RowVecTrans(ctx.T_t1_ICR(:,:,ii), ...
        RowVecTrans(ctx.T_Pam_inv(:,:,ii), raw(1,:)));
    bridgeT1(ii) = t1ChordClear(p1T1(1:2), raw(9,1:2), ctx.geo, tol, 8);
end

iBF = find(bridgeFemur, 1, 'first');
iBT = find(bridgeT1, 1, 'first');
fprintf('Femur chord p1->p9 clear from %s onward (%d of %d indices)\n', ...
    angleText(ctx.phiD, iBF), sum(bridgeFemur), ctx.N)
fprintf('T1 chord p1->p9 clear from %s onward (%d of %d indices)\n', ...
    angleText(ctx.phiD, iBT), sum(bridgeT1), ctx.N)

% Where both the bridge is clear and p3:p7 are gone (bridge test condition).
bothGone = all(~ri.active(3:7,:), 1).';
iBridgeEligible = find(bridgeFemur & bridgeT1 & bothGone, 1, 'first');
fprintf('Bridge-eligible (p3:p7 inactive AND both chords clear) from %s\n', ...
    angleText(ctx.phiD, iBridgeEligible))

%% 5. Sample-angle margin tables
sampleAnglesD = [-120, -100, -80, -69.8, -60, -44, -30, -9.19, 0, 10];
sampleIdx = zeros(size(sampleAnglesD));
for kk = 1:numel(sampleAnglesD)
    [~, sampleIdx(kk)] = min(abs(ctx.phiD - sampleAnglesD(kk)));
end

fprintf('\n========== MARGINS AT SAMPLE ANGLES (x0) ==========\n')
fprintf('margin = builder rule, principal atan2 values, deg.\n')
fprintf('Positive margin means the angle gate allows removal.\n')
fprintf('%-5s |', 'Point')
for kk = 1:numel(sampleIdx)
    fprintf(' %8.2f |', ctx.phiD(sampleIdx(kk)))
end
fprintf('\n')
for j = 2:8
    fprintf('p%-4d |', j)
    for kk = 1:numel(sampleIdx)
        if isfinite(mRaw(j,sampleIdx(kk)))
            fprintf(' %8.2f |', mRaw(j,sampleIdx(kk)))
        else
            fprintf(' %8s |', '--')
        end
    end
    fprintf('\n')
end

fprintf('\n========== BYPASS CHORD CLEAR AT SAMPLE ANGLES (x0) ==========\n')
fprintf('1 = physically clear bypass chord, 0 = blocked, - = not evaluated\n')
fprintf('%-5s |', 'Point')
for kk = 1:numel(sampleIdx)
    fprintf(' %8.2f |', ctx.phiD(sampleIdx(kk)))
end
fprintf('\n')
for j = 2:8
    fprintf('p%-4d |', j)
    [clearMask, evalMask] = rowBypassClearance(ri, ctx, j);
    for kk = 1:numel(sampleIdx)
        ii = sampleIdx(kk);
        if evalMask(ii)
            fprintf(' %8d |', clearMask(ii))
        else
            fprintf(' %8s |', '--')
        end
    end
    fprintf('\n')
end

%% 6. Scan p1 and pEnd across optimizer bounds
fprintf('\n========== P1 / PEND BOUND SCAN ==========\n')
fprintf('One coordinate at a time; all others held at x0.\n')
fprintf('elim = rows eliminated with angles; step = max same-step removals;\n')
fprintf('badElim = eliminations while the bypass chord was blocked.\n\n')

scanSpec = { ...
    'p1.x',   1, ctx.lb(1), ctx.ub(1); ...
    'p1.y',   2, ctx.lb(2), ctx.ub(2); ...
    'pEnd.x', 4, ctx.lb(4), ctx.ub(4); ...
    'pEnd.y', 5, ctx.lb(5), ctx.ub(5)};

nScan = 7;
x0 = ctx.x0;

for cSpec = 1:size(scanSpec,1)

    label = scanSpec{cSpec,1};
    slot  = scanSpec{cSpec,2};
    lo    = scanSpec{cSpec,3};
    hi    = scanSpec{cSpec,4};

    fprintf('--- %s sweep, %d designs ---\n', label, nScan)

    for cVal = 1:nScan

        val = lo + (hi-lo)*(cVal-1)/(nScan-1);

        xTry = x0;
        xTry(slot) = val;

        try
            [~, ~, infoS] = buildDistalRingLocation20mm( ...
                xTry(1:3), xTry(4:6), xTry(8), ctx);
        catch ME
            fprintf('%8.4f | builder error: %s\n', val, ME.message)
            continue
        end

        nStep = max(sum(infoS.active(:,1:end-1) & ~infoS.active(:,2:end), 1));

        nBad = 0;
        badText = '';
        for j = 2:8
            [clearMask, evalMask] = rowBypassClearance(infoS, ctx, j);
            bad = evalMask & ~clearMask & gateOpens(mR2row(infoS, ctx, j));
            if any(bad)
                nBad = nBad + 1;
                iBad = find(bad, 1, 'first');
                badText = [badText, sprintf(' p%d@%+.1f', ...
                    j, ctx.phiD(iBad))]; %#ok<AGROW>
            end
        end

        elimText = '';
        for j = 2:8
            if isfinite(infoS.eliminatedAngleD(j))
                elimText = [elimText, sprintf(' p%d@%+.1f', ...
                    j, infoS.eliminatedAngleD(j))]; %#ok<AGROW>
            end
        end
        if isempty(elimText)
            elimText = ' none';
        end

        fprintf('%8.4f |%s | step %d | blocked-removals:%s\n', ...
            val, elimText, nStep, badText)
    end
    fprintf('\n')
end

fprintf('DEBUG_RUN_COMPLETE\n')

%% ======================================================================
%% Local functions
%% ======================================================================

function margins = mR2row(infoS, ctx, j)
% Convenience: reconstructed margins for one row over the whole sweep.

m = routeMargins(infoS, ctx);
margins = m(j,:);

end


function [mRaw, mWrap, aBigAll, aSmallAll, actMask] = routeMargins(infoS, ctx)
% Reconstruct the builder's native-frame rotated-vector rule for every
% sweep index and relevant optional row:
%   femur rows (2:5): vectors from the previous active row, rotated +90 deg,
%     margin = aBig - aSmall (principal atan2 values, degrees).
%   tibia rows (6:8): vectors from the next active row, rotated -90 deg,
%     margin = aSmall - aBig.
% Positive margin means the angle gate allows removal. mRaw and mWrap are
% identical because the rule uses principal atan2 values directly.
%
% A row is evaluated at index ii if it is active there OR if that index is
% its recorded elimination index.

N = size(infoS.active, 2);
rows = 9;

mRaw      = nan(rows, N);
mWrap     = nan(rows, N);
aBigAll   = nan(rows, N);
aSmallAll = nan(rows, N);

actMask = infoS.active;

elimIdx = infoS.eliminatedSweepIndex;

for ii = 1:N

    act = infoS.active(:,ii);
    raw = infoS.raw(:,:,ii);

    for j = 2:8

        evalHere = act(j) || (isfinite(elimIdx(j)) && elimIdx(j) == ii);
        if ~evalHere
            continue
        end

        iPrev = find(act(1:j-1), 1, 'last');
        iNext = j + find(act(j+1:9), 1, 'first');
        if isempty(iPrev) || isempty(iNext)
            continue
        end

        A = raw(iPrev,:);
        B = raw(iNext,:);
        C = raw(j,:);

        if j <= 5
            if iNext >= 6
                B = RowVecTrans(ctx.T_Pam(:,:,ii)*ctx.T_ICR_t1(:,:,ii), B);
            end
            vBig   = B(1:2) - A(1:2);
            vSmall = C(1:2) - A(1:2);
            aBig   = atan2d(-vBig(2),   vBig(1));    % +90 rotation
            aSmall = atan2d(-vSmall(2), vSmall(1));
            marginD = aBig - aSmall;
        else
            if iPrev <= 5
                A = RowVecTrans(ctx.T_t1_ICR(:,:,ii), ...
                    RowVecTrans(ctx.T_Pam_inv(:,:,ii), A));
            end
            vBig   = A(1:2) - B(1:2);
            vSmall = C(1:2) - B(1:2);
            aBig   = atan2d(vBig(2),   -vBig(1));    % -90 rotation
            aSmall = atan2d(vSmall(2), -vSmall(1));
            marginD = aSmall - aBig;
        end

        mRaw(j,ii)  = marginD;
        mWrap(j,ii) = marginD;
        aBigAll(j,ii)   = aBig;
        aSmallAll(j,ii) = aSmall;
    end
end

end


function [clearMask, evalMask] = rowBypassClearance(infoS, ctx, j)
% For one optional row j, evaluate at every sweep index whether the bypass
% chord that would replace it is physically collision-free, using the same
% nearest-active-row anchors as the builder's elimination test.
%
% The row is evaluated where it is active or at its elimination index.

N = size(infoS.active, 2);
clearMask = false(1,N);
evalMask  = false(1,N);

tol = 1e-8;
geo = ctx.geo;

for ii = 1:N

    act = infoS.active(:,ii);
    raw = infoS.raw(:,:,ii);

    isElimHere = isfinite(infoS.eliminatedSweepIndex(j)) && ...
        infoS.eliminatedSweepIndex(j) == ii;

    if ~(act(j) || isElimHere)
        continue
    end

    iPrev = find(act(1:j-1), 1, 'last');
    iNext = j + find(act(j+1:9), 1, 'first');

    if isempty(iPrev) || isempty(iNext)
        continue
    end

    if j <= 5
        % ---- femur-side candidate ----
        A = raw(iPrev,1:2);
        B = raw(iNext,:);
        if iNext >= 6
            B = RowVecTrans(ctx.T_Pam(:,:,ii)*ctx.T_ICR_t1(:,:,ii), B);
        end

        clearMask(ii) = femurChordClear(A, B(1:2), geo, tol, j);
        evalMask(ii) = true;

    else
        % ---- t1-side candidate ----
        A = raw(iPrev,:);
        if iPrev <= 5
            A = RowVecTrans(ctx.T_t1_ICR(:,:,ii), ...
                RowVecTrans(ctx.T_Pam_inv(:,:,ii), A));
        end
        B = raw(iNext,1:2);

        clearMask(ii) = t1ChordClear(A(1:2), B, geo, tol, j);
        evalMask(ii) = true;
    end
end

end


function clear = femurChordClear(A, B, geo, tol, j)
% Femur-side bypass chord clearance, replicating the builder's gates with
% the bypass tolerance and endpoint trimming.

gt = geo.bypassTol;
[A, B] = trimChordEnds(A, B);

switch j
    case {2, 3}
        clear = ...
            ~segmentPenetratesCircle(A, B, ...
                geo.femurCylCenter, geo.femurCylClearRadius + gt, tol) && ...
            ~segmentPenetratesFemurOffset(A, B, geo, tol) && ...
            ~segmentIntersectsVerticalSpan(A, B, ...
                geo.femurLineX - gt, geo.femurLineY, tol);
    case {4, 5}
        clear = ~segmentPenetratesFemurOffset(A, B, geo, tol);
    otherwise
        clear = true;
end

end


function clear = t1ChordClear(A, B, geo, tol, j)
% Tibia-side bypass chord clearance, replicating the builder's gates with
% the bypass tolerance and endpoint trimming.

gt = geo.bypassTol;
[A, B] = trimChordEnds(A, B);

switch j
    case 6
        clear = ~segmentPenetratesCircle(A, B, ...
            geo.tibiaUpperCenter, geo.tibiaUpperClearRadius + gt, tol);
    case 7
        clear = ~segmentIntersectsVerticalSpan(A, B, ...
            geo.tibiaWallX - gt, geo.tibiaWallY, tol);
    case 8
        clear = ~segmentPenetratesCircle(A, B, ...
                geo.tibiaLowerCenter, geo.tibiaLowerClearRadius + gt, tol) && ...
            ~segmentIntersectsVerticalSpan(A, B, ...
                geo.tibiaWallX - gt, geo.tibiaWallY, tol);
    otherwise
        clear = true;
end

end


function hit = segmentPenetratesFemurOffset(A, B, geo, tol)
% Faithful replication of the builder's clipped-offset-polygon test.

poly = geo.femurOffsetBoundaryGate;

[inA, onA] = inpolygon(A(1), A(2), poly(:,1), poly(:,2));
[inB, onB] = inpolygon(B(1), B(2), poly(:,1), poly(:,2));

if (inA && ~onA) || (inB && ~onB)
    hit = true;
    return
end

tCuts = segmentPolygonIntersectionParams(A, B, poly, tol);

if isempty(tCuts)
    hit = false;
    return
end

tCuts = sort(unique([0; tCuts(:); 1]));

hit = false;

for k = 1:numel(tCuts)-1

    ta = tCuts(k);
    tb = tCuts(k+1);

    if tb-ta < 1e-10
        continue
    end

    tm = 0.5*(ta+tb);
    Pm = A + tm*(B-A);

    [in, on] = inpolygon(Pm(1), Pm(2), poly(:,1), poly(:,2));

    if in && ~on
        hit = true;
        return
    end
end

end


function tVals = segmentPolygonIntersectionParams(A, B, poly, tol)

tVals = zeros(0,1);
n = size(poly,1);

for e = 1:n
    C = poly(e,:);
    D = poly(mod(e,n)+1,:);
    t = segmentSegmentIntersectionParam(A, B, C, D, tol);
    if isfinite(t)
        tVals(end+1,1) = t; %#ok<AGROW>
    end
end

end


function t = segmentSegmentIntersectionParam(A, B, C, D, tol)

r = B-A;
s = D-C;
rxs = r(1)*s(2) - r(2)*s(1);
qmp = C-A;

t = NaN;

if abs(rxs) <= tol
    return
end

tt = (qmp(1)*s(2) - qmp(2)*s(1))/rxs;
u  = (qmp(1)*r(2) - qmp(2)*r(1))/rxs;

if tt >= -tol && tt <= 1+tol && u >= -tol && u <= 1+tol
    t = min(max(tt,0),1);
end

end


function hit = segmentPenetratesCircle(A, B, C, R, tol)

d = B-A;
d2 = dot(d,d);

if d2 < 1e-18
    hit = norm(A-C) < R-tol;
    return
end

t = dot(C-A,d)/d2;
t = min(max(t,0),1);

q = A + t*d;

hit = norm(q-C) < R-tol;

end


function hit = segmentIntersectsVerticalSpan(A, B, xValue, yRange, tol)

dx = B(1)-A(1);

if abs(dx) < 1e-14
    hit = false;
    return
end

t = (xValue-A(1))/dx;

if t <= tol || t >= 1-tol
    hit = false;
    return
end

y = A(2)+t*(B(2)-A(2));
hit = y >= min(yRange)-tol && y <= max(yRange)+tol;

end


function [iRaw, iWrap, nFlip] = gateOpenIndices(mR, mW, act)
% First sweep index at which the gate opens for one row under both recorded
% variants (identical by construction) and the seam-flip count (always 0).

both = isfinite(mR) & isfinite(mW);

openRaw  = both & gateOpens(mR);
openWrap = both & gateOpens(mW);

iRaw  = find(openRaw,  1, 'first');
iWrap = find(openWrap, 1, 'first');

nFlip = 0;

end


function s = angleText(phiD, ii)

if isempty(ii) || ~isfinite(ii)
    s = 'never';
else
    s = sprintf('%7.2f deg', phiD(ii));
end

end


function [A, B] = trimChordEnds(A, B)

d = B - A;
L = norm(d);

if L < 4e-3
    return
end

e = min(0.45, 2e-3/L);
A = A + e*d;
B = B - e*d;

end


function opens = gateOpens(marginD)
% Removal-rule replica: the angle gate opens when the principal-value
% margin is positive (big more counterclockwise on the femur side, big
% more clockwise on the tibia side, per the builder convention).

opens = marginD > 0;

end
