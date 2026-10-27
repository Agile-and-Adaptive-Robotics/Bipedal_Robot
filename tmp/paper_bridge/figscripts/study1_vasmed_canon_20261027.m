% study1_vasmed_canon_20261027.m -- STUDY 1 (2026-10-27, easteregg2)
% CANONICAL FUNCTION-FILE RESTRUCTURE of study1_vasmed_minBPA_20261027.m
% (kept untouched as the failure record; see its study1_run*_20261027.log).
%
% ROOT CAUSE fixed here: the 10-27 repair inserted an extra end at the old
% line 294. evalConfig already closes at its own end (old line 279), so the
% inserted end terminated the PRIMARY function early; the entire search/
% report/save tail then sat at file scope outside any function, which
% forces MATLAB to classify the file as a SCRIPT (statements outside
% functions are illegal in a function file). Invoked as a script, MATLAB
% only DEFINES the primary function (never calls it) and then executes
% the tail in an empty workspace -> 'Lchar' undefined at line 297; called
% by the renamed primary name -> not found (filename mismatch); run() ->
% script-classification errors. Same numbers, one clean layout: primary
% function with ONLY body statements, ALL helpers as LOCAL functions after
% its end, explicit argument passing (a P struct) replacing nested
% workspace sharing.
%
% ONE numeric-neutral repair beyond the end fix: cableFirstLeg's interior
% used ROW vectors, making the dot product rel*u and the rotation
% rot*base dimensionally invalid (1x2 * 1x2 / 2x2 * 1x2). The function
% never executed in ANY earlier run (the 01:55 run died at its paren-dot
% indexing bug first). Points are now forced to 2x1 columns inside; the
% math (dot-product projection, external tangent, same iteration) is
% unchanged.
%
% TWO further numeric-neutral repairs in never-executed tail code (first
% canonical run, 2026-10-27): (1) info.tendonMax is an N-pose VECTOR of
% one pose-independent scalar -- the hard reject now uses min(), exactly
% nonlconExt20mm's cTendonMax convention; (2) ctx.phiD(1) =
% -119.999994 deg (kneeMin truncated at 7 decimals) leaves the -120 deg
% surface node 6e-6 deg outside the grid -- every node interp now uses
% linear 'extrap' (the endpoint value to ~1e-7 relative).
%
% Study ask unchanged (BiPulley studies directive, 2026-10-27):
%   Minimum 20 mm BPAs (parallel, up to 3 tackle line parts each) to
%   meet/exceed the human vastus medialis SINGLE-muscle torque surface
%   over the graded knee RoM, per moment-arm method:
%     (a) OpenSim ds/dtheta arm (vas_med_surface_5node_20261027.csv
%         arm_m column, regenerated 2026-10-27 from pinned opensim 4.4.1
%         stock gait2392_simbody routing, computeMomentArm x Fmax 1294 N);
%     (b) OUR cable arm: pEnd -> old proximal ring p6 (ctxCamBest
%         routeSeed row 6) as a cable contacting the tibia clearance
%         circles 0.00075 m off the surfaces, moment about the robot
%         knee ICR from the measured t1-ICR table (fcn13/fcn14 of
%         buildKneeExtContext20mm.m).
%   Search: nBPA {1,2,3} x tackle G {1,2,3}; sizes [rest; tendon]
%   re-solved per configuration; STIFFNESS-UNAWARE first, then
%   STIFFNESS-AWARE (adopted extensor pick 32: Xi0 -0.00637802 m,
%   Xi1 39978.7 N/m axial, Xi2 14733.8 N/m normal; chi3 NOT carried --
%   the pulley closure family has no Xi3 term, stated). Tendon bounds
%   [0, 0.3] m with the routeInfo.tendonMax hard reject; rest within
%   [0.30, 0.90]*Lchar. Pressure 620 kPa; KMAX = 0.255.
%
% EQUILIBRIUM GATE: the STIFFNESS-AWARE mode evaluates THROUGH
% MonoPam_pulley itself (machinery of record -- see the evalConfig
% note; three 1-D port attempts failed the draft's 1e-6 closure-port
% gate at the 25 collapsed-segment poses, max residual 0.99%). The gate
% asserts the study's class-path plumbing (cached route, arg order,
% pulley config) against a FRESH MonoPam_pulley at (nB=1, G=1) with
% the adopted chi, to 1e-6 relative, before any search. The UNAWARE
% mode keeps the chi-zeroed 1-D closure port (no class analog: the
% constructor offers no kSpr-off switch).
%
% Outputs (unchanged names):
%   Results/BiPulley_surveys/VasMed_minBPA_study1_20261027.mat
%   tmp/paper_bridge/capability_table_20261022/VasMed_minBPA_study1_summary_20261027.csv
%
% Invocation (canonical): -batch "cd(figscripts); study1_vasmed_canon_20261027"

function study1_vasmed_canon_20261027()

R = 'D:\GitHub\Bipedal_Robot';
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));
capDir = fullfile(R, 'tmp', 'paper_bridge', 'capability_table_20261022');
resDir = fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'BiPulley_surveys');

%% ---------------- Fixed data of record --------------------------------
C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), 'xFullBest', 'xBestCam', ...
    'ctxCamBest', 'XiUsed');
xRec = C.xFullBest(:).';           % [p1(1:3) pEnd(1:3) rest tendon]
p1Rec = xRec(1:3);
pEndRec = xRec(4:6);
restRec = xRec(7);
tendonRec = xRec(8);
ctx = C.ctxCamBest;
Xi = C.XiUsed(:).';                % [Xi0 Xi1 Xi2 Xi3] adopted pick 32
assert(numel(Xi) == 4, 'XiUsed layout');
% ctxCamBest.routeSeed is the context's ORIGINAL seed table (its row 1 is
% the pre-pin p1 [0.040 0.035 0]); the design-of-record p1 comes from
% xFullBest ([0.040 0.137 0], pinned in the cam ctx lb/ub). The old
% proximal ring (row 6) is read from the context seed table.
ringOld = ctx.routeSeed(6, :);     % OLD PROXIMAL RING (stock pre-cam p6)
fprintf('Design of record: p1 [%s] pEnd [%s] rest %.6f tendon %.6f\n', ...
    num2str(p1Rec, 5), num2str(pEndRec, 5), restRec, tendonRec);
fprintf('Adopted chi (pick 32): Xi0 %.6g Xi1 %.6g Xi2 %.6g Xi3 %.6g\n', Xi);
fprintf('Old proximal ring (routeSeed row 6, t1): [%s]\n', num2str(ringOld, 5));

Dia = 20; Pres = 620; KMAX = 0.255; fitting = 0.021; wraps = 3;

%% ---------------- Human target + OpenSim arm (a) ----------------------
T5 = readtable(fullfile(capDir, 'vas_med_surface_5node_20261027.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
[~, ord] = sort(T5.knee_angle_deg(:));
nodeDeg = T5.knee_angle_deg(ord);           % [-120 -87.5 -55 -22.5 10]
tauH = T5.tau_h_Nm(ord);                    % +30.2 .. +58.7 N*m (ext +)
armOsim = T5.arm_m(ord);                    % OpenSim ds/dtheta arm (a)
nNode = numel(nodeDeg);
fprintf('vas_med nodes [%s] deg\n', num2str(nodeDeg.'));
fprintf('  tau_h  [%s] N*m\n', num2str(tauH.', 5));
fprintf('  arm(a) [%s] m\n', num2str(armOsim.', 5));

%% ---------------- OUR cable arm (b) ------------------------------------
% Path: pEnd (t1 frame, the given muscle insertion of the design of
% record) -> old proximal ring p6 (t1 frame), a cable contacting the
% tibia clearance circles 0.00075 m off the physical surfaces (upper cyl
% r 7.5 mm, lower cyl r 5 mm). Taut tangent construction from pEnd; the
% moment arm = signed z moment per unit tension, M = (pEnd - ICR) x u,
% with u = unit direction from pEnd toward the first contact (or the
% ring if the chord clears) and ICR from the measured t1-ICR table.
cableOff = 0.00075;
cylUp = struct('c', ctx.geo.tibiaUpperCenter, 'r', ...
    ctx.geo.tibiaUpperRadius + cableOff);
cylLo = struct('c', ctx.geo.tibiaLowerCenter, 'r', ...
    ctx.geo.tibiaLowerRadius + cableOff);
phiD = ctx.phiD(:);                          % 100-pose grid, -120..10 deg
armCable = zeros(size(phiD));
contactLog = zeros(numel(phiD), 2);
for i = 1:numel(phiD)
    th = phiD(i) * pi / 180;
    ICR_t1 = [t1icrVal('x', th), t1icrVal('y', th), 0];
    [uFirst, nC] = cableFirstLeg(pEndRec(1:2), ringOld(1:2), {cylUp, cylLo});
    contactLog(i, :) = nC;
    r = pEndRec - ICR_t1;
    m = cross(r, [uFirst(1), uFirst(2), 0]);
    armCable(i) = m(3);                      % signed per-unit-torque, m
end
% 'extrap' repair (canonical restructure): ctx.phiD(1) = -119.999994 deg
% (kneeMin -2.0943951 rad is truncated at 7 decimals), so the -120 deg
% surface node sits 6e-6 deg below the grid; without extrapolation every
% node interp returns NaN there. Linear extrapolation over that gap is
% the endpoint value to ~1e-7 relative. Same repair on every node interp
% below (tension, arm, strain, slack).
armCableNode = interp1(phiD, armCable, nodeDeg, 'linear', 'extrap');
fprintf('arm(b) cable (signed, nodes) [%s] m\n', num2str(armCableNode.', 5));
fprintf('  contact-active poses: upper %d/%d, lower %d/%d\n', ...
    nnz(contactLog(:, 1)), numel(phiD), nnz(contactLog(:, 2)), numel(phiD));

%% ---------------- Route cache + evaluator ------------------------------
routeCache = containers.Map('KeyType', 'double', 'ValueType', 'any');
[Location0, ~, info0] = getRoute(tendonRec, p1Rec, pEndRec, ctx, routeCache); %#ok<ASGLU>
mL0 = routeLengths(Location0, ctx);
Lchar = max(mL0);
fprintf('Lchar (max rigid route length over grid) = %.6f m; rest bounds [%.4f %.4f]\n', ...
    Lchar, 0.30 * Lchar, 0.90 * Lchar);

FmaxOf = @(rest) maxBPAforce(rest, '20', Pres);

% Shared evaluator context (explicit argument passing replaces the nested
% workspace sharing of the broken draft; values identical).
P = struct('ctx', ctx, 'Xi', Xi, 'Dia', Dia, 'Pres', Pres, 'KMAX', KMAX, ...
    'fitting', fitting, 'wraps', wraps, 'phiD', phiD, 'nodeDeg', nodeDeg, ...
    'tauH', tauH, 'armOsim', armOsim, 'armCable', armCable, ...
    'p1Rec', p1Rec, 'pEndRec', pEndRec, 'routeCache', routeCache, ...
    'FmaxOf', FmaxOf);

%% ---------------- Equilibrium gate vs MonoPam_pulley -------------------
gate = evalConfig([restRec, tendonRec], 1, 1, 'aware', P);
assert(gate.ok, 'gate eval failed');
pulCfg = struct('nPulleyBPA', 1, 'tackleLineParts', 1, 'gain', 1, ...
    'routingMode', 'moving_via');
obj = MonoPam_pulley('gate', getRoute(tendonRec, p1Rec, pEndRec, ctx, ...
    routeCache), ctx.CrossPoint, Dia, ctx.T_Pam, restRec, ...
    restRec * (1 - KMAX), tendonRec, fitting, Pres, Xi(1), Xi(2), Xi(3), ...
    wraps, pulCfg);
ftClass = double(obj.Ftendon(:));
ftMine = gate.tensionGrid;
relErr = max(abs(ftMine - ftClass) ./ max(1, abs(ftClass)));
fprintf('EQUILIBRIUM GATE (1,1 aware): max rel err vs MonoPam_pulley = %.3g\n', relErr);
assert(relErr < 1e-6, 'class-path gate failed: evalConfig vs fresh MonoPam_pulley');

%% ---------------- Search: 9 configs x 2 modes x 2 arm methods ----------
lb = [0.30 * Lchar; 0];
ub = [0.90 * Lchar; 0.3];
seeds = [0.60 * Lchar, 0.04; 0.30 * Lchar, 0.02; 0.90 * Lchar, 0.02; ...
    restRec, tendonRec];

opts = optimoptions('patternsearch', 'Display', 'none', ...
    'MeshTolerance', 1e-4, 'StepTolerance', 1e-4, ...
    'MaxFunctionEvaluations', 100);
rows = {};
tAll = tic;
for md = {'unaware', 'aware'}
    for mh = {'a', 'b'}
        for nB = 1:3
            for G = 1:3
                fSeed = zeros(size(seeds, 1), 1);
                for sIdx = 1:size(seeds, 1)
                    fSeed(sIdx) = objFun(seeds(sIdx, :), nB, G, md{1}, mh{1}, P);
                end
                [~, ixSeed] = min(fSeed);
                [xB, fB] = patternsearch(@(x) objFun(x, nB, G, md{1}, mh{1}, P), ...
                    seeds(ixSeed, :), [], [], [], [], lb, ub, [], opts);
                E = evalConfig(xB, nB, G, md{1}, P);
                if strcmp(mh{1}, 'a'), tq = E.tau_a; else, tq = E.tau_b; end
                meet = isfinite(tq) & sign(tq) == sign(tauH) & ...
                    abs(tq) >= abs(tauH) - 1e-12;
                margin = tq - tauH;
                minPct = 100 * (min(abs(tq) ./ abs(tauH)) - 1);
                rows(end + 1, :) = {md{1}, mh{1}, nB, G, xB(1), xB(2), ...
                    sum(meet), nNode, min(margin), minPct, ...
                    max(E.infeas), max(E.slack), fB, E.tension, tq}; %#ok<AGROW>
                fprintf(['%-8s arm(%s) nB%d G%d: meet %d/%d, min margin ' ...
                    '%+.2f N*m (%+.1f%%), rest %.4f tendon %.4f\n'], ...
                    md{1}, mh{1}, nB, G, sum(meet), nNode, min(margin), ...
                    minPct, xB(1), xB(2));
            end
        end
    end
end
fprintf('Search wall time: %.1f min\n', toc(tAll) / 60);

%% ---------------- Minimum-config report + save -------------------------
hdr = {'mode', 'arm_method', 'nBPA', 'G', 'rest_m', 'tendon_m', ...
    'meetCells', 'nCells', 'minMargin_Nm', 'minMarginPct', ...
    'anyInfeasible', 'anySlack', 'score'};
pickRow = zeros(4, 1); pickTag = cell(4, 1); pk = 0;
for md = {'unaware', 'aware'}
    for mh = {'a', 'b'}
        pk = pk + 1;
        sel = find(strcmp(rows(:, 1), md{1}) & strcmp(rows(:, 2), mh{1}));
        met = sel([rows{sel, 7}] == nNode);
        if ~isempty(met)
            [~, ix] = min([rows{met, 3}] + 0.01 * [rows{met, 4}]);
            pickRow(pk) = met(ix); pickTag{pk} = 'MEETS';
        else
            m7 = [rows{sel, 7}];
            cand = sel(m7 == max(m7));
            [~, ix] = max([rows{cand, 10}]);
            pickRow(pk) = cand(ix); pickTag{pk} = 'SHORTFALL';
        end
        pr = pickRow(pk);
        fprintf(['STUDY1 RESULT %-8s arm(%s): %s -> nBPA %d, G %d ' ...
            '(meet %d/%d, min margin %+.2f N*m / %+.1f%%)\n'], md{1}, ...
            mh{1}, pickTag{pk}, rows{pr, 3}, rows{pr, 4}, rows{pr, 7}, ...
            nNode, rows{pr, 9}, rows{pr, 10});
    end
end

marginsPick = struct();
pk = 0;
for md = {'unaware', 'aware'}
    for mh = {'a', 'b'}
        pk = pk + 1;
        pr = pickRow(pk);
        tq = rows{pr, 15};
        marginsPick.(sprintf('%s_%s', md{1}, mh{1})) = struct( ...
            'tau_r', tq, 'margin', tq - tauH, 'marginPct', ...
            100 * (abs(tq) ./ abs(tauH) - 1), 'meet', isfinite(tq) & ...
            sign(tq) == sign(tauH) & abs(tq) >= abs(tauH) - 1e-12, ...
            'tension', rows{pr, 14}, 'nBPA', rows{pr, 3}, 'G', rows{pr, 4}, ...
            'rest', rows{pr, 5}, 'tendon', rows{pr, 6});
        fprintf('  pick %-8s arm(%s) per-node margins N*m: [%s]\n', ...
            md{1}, mh{1}, num2str((tq - tauH).', '%.2f '));
    end
end

save(fullfile(resDir, 'VasMed_minBPA_study1_20261027.mat'), 'rows', 'hdr', ...
    'pickRow', 'pickTag', 'marginsPick', 'tauH', 'nodeDeg', 'armOsim', ...
    'armCableNode', 'armCable', 'ringOld', 'pEndRec', 'p1Rec', ...
    'restRec', 'tendonRec', 'Xi', 'Lchar', 'cableOff', 'contactLog', ...
    'relErr', 'tendonRec');

fid = fopen(fullfile(capDir, 'VasMed_minBPA_study1_summary_20261027.csv'), 'w');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for k = 1:size(rows, 1)
    fprintf(fid, '%s,%s,%d,%d,%.6g,%.6g,%d,%d,%.6g,%.6f,%d,%d,%.6g\n', ...
        rows{k, 1}, rows{k, 2}, rows{k, 3}, rows{k, 4}, rows{k, 5}, ...
        rows{k, 6}, rows{k, 7}, rows{k, 8}, rows{k, 9}, rows{k, 10}, ...
        rows{k, 11}, rows{k, 12}, rows{k, 13});
end
fprintf(fid, '# picks of record (minimum config per stiffness mode x arm method)\n');
for k = 1:4
    pr = pickRow(k);
    fprintf(fid, '# PICK,%s,%s,%s,nBPA %d,G %d,meet %d/%d,minMarginNm %.4f,minMarginPct %.4f\n', ...
        rows{pr, 1}, rows{pr, 2}, pickTag{k}, rows{pr, 3}, rows{pr, 4}, ...
        rows{pr, 7}, rows{pr, 8}, rows{pr, 9}, rows{pr, 10});
end
fclose(fid);
fprintf('STUDY1 DONE\n');
end

%% ================= LOCAL FUNCTIONS (all helpers, explicit args) =======
%% --------------- measured t1-ICR table (fcn13/fcn14) -------------------
function v = t1icrVal(xy, th)
% The measured t1-ICR table (fcn13/fcn14) pchip-interpolated at knee
% angle th (rad). Tables verbatim from buildKneeExtContext20mm.m:30-40.
ka = [0.17, 0.09, 0.03, 0.00, -0.09, -0.17, -0.26, -0.52, -0.79, ...
    -1.05, -1.31, -1.57, -1.83, -2.09, -2.36, -2.62];
if strcmp(xy, 'x')
    tab = [29.66, 28.54, 27.86, 27.40, 26.23, 25.03, 23.81, 20.03, ...
        16.17, 12.34, 8.67, 5.24, 2.04, -1.01, -4.10, -7.58] / 1000;
else
    tab = [25.97, 25.74, 25.61, 25.53, 25.35, 25.19, 25.03, 24.57, ...
        24.04, 23.39, 22.66, 21.93, 21.32, 20.99, 21.20, 22.33] / 1000;
end
v = interp1(ka, tab, th, 'pchip');
end

%% --------------- taut-cable first leg (column vectors) -----------------
function [u, nC] = cableFirstLeg(pFrom, pTo, circs)
% Unit direction of the FIRST leg of the taut cable pFrom -> pTo that
% contacts the clearance circles it cannot pass through. Iterative: take
% the first circle hit along the current leg, replace the leg with the
% tangent from the current point to that circle (tangent point chosen on
% the side the straight chord passes), repeat. nC = contacts used.
% Points forced to 2x1 columns (canonical restructure): the draft's row
% vectors made proj = rel*u and u = rot*base dimensionally invalid.
cur = pFrom(:);
dest = pTo(:);
used = false(1, numel(circs));
nC = zeros(1, numel(circs));
u = (dest - cur) / norm(dest - cur);
for it = 1:numel(circs)
    d = dest - cur;
    L = norm(d); u = d / L;
    hit = 0; hitT = inf;
    for k = 1:numel(circs)
        if used(k), continue; end
        rel = circs{k}.c(:) - cur;
        proj = rel' * u;
        if proj <= 0 || proj >= L, continue; end
        perp = norm(rel - proj * u);
        if perp < circs{k}.r && proj < hitT
            hitT = proj; hit = k;
        end
    end
    if hit == 0, break; end
    c = circs{hit}.c(:); r = circs{hit}.r;
    dd = norm(c - cur);
    if dd <= r, break; end                 % degenerate; keep last leg
    % External tangent from cur to circle (c, r), on the side of the
    % chord (sign of the cross product chord x (c - cur)).
    base = (c - cur) / dd;
    side = sign((dest(1) - cur(1)) * (c(2) - cur(2)) - ...
        (dest(2) - cur(2)) * (c(1) - cur(1)));
    if side == 0, side = 1; end
    alpha = acos(min(max(r / dd, -1), 1));
    rot = [cos(side * alpha), -sin(side * alpha); ...
        sin(side * alpha), cos(side * alpha)];
    u = rot * base;                        % tangent direction from cur
    cur = cur + u * sqrt(max(dd^2 - r^2, 0));   % tangent point
    used(hit) = true; nC(hit) = 1;
end
end

%% --------------- cached route builder ---------------------------------
function [Location, bendM, info] = getRoute(tendon, p1Rec, pEndRec, ctx, ...
    routeCache)
key = round(tendon * 1e9);
if isKey(routeCache, key)
    v = routeCache(key);
    Location = v{1}; bendM = v{2}; info = v{3};
    return
end
[Location, bendM, info] = buildDistalRingLocation20mm(p1Rec, ...
    pEndRec, tendon, ctx);
routeCache(key) = {Location, bendM, info};
end

%% --------------- muscle lengths (muscleLengthNormal verbatim) ---------
function mL = routeLengths(Location, ctx)
% muscleLengthNormal verbatim (buildKneeExtContext20mm local):
% segment (Cross-1)->Cross crosses frames; the crossing row is
% mapped into the parent frame by T_Pam.
N = size(Location, 3);
mL = zeros(N, 1);
for ii = 1:N
    tot = 0;
    for j = 1:size(Location, 1) - 1
        pA = Location(j, :, ii);
        pB = Location(j + 1, :, ii);
        if j + 1 == ctx.CrossPoint
            pB = RowVecTrans(ctx.T_Pam(:, :, ii), pB);
        end
        tot = tot + norm(pA - pB);
    end
    mL(ii) = tot;
end
end

%% --------------- crossed-segment span (MonoPam_pulley verbatim) -------
function span = spanSeg(Location, ctx)
% MonoPam_pulley.updatePulleyGeometry verbatim: the rigid tendon span is
% the SINGLE crossed segment -- exit row (PulleyExitIndex = Cross-1)
% mapped into the child frame, to the Cross row. (The draft's
% multi-segment exit->insertion span was wrong; gate-rejected.)
N = size(Location, 3);
span = zeros(N, 1);
for ii = 1:N
    exitT1 = RowVecTrans(ctx.T_Pam(:, :, ii) \ eye(4), ...
        Location(ctx.CrossPoint - 1, :, ii));
    span(ii) = norm(exitT1 - Location(ctx.CrossPoint, :, ii));
end
end


%% --------------- per-configuration equilibrium evaluator --------------
function E = evalConfig(x, nB, G, mode, P)
rest = x(1); tendon = x(2);
E = struct('ok', false);
if tendon < 0 || tendon > 0.3
    return
end
try
    [Location, ~, info] = getRoute(tendon, P.p1Rec, P.pEndRec, P.ctx, ...
        P.routeCache);
catch
    return
end
% info.tendonMax is the N-pose vector; tendonLimit20mm(pEnd, geo) is
% pose-independent, so every element is the same scalar. min() over the
% vector reproduces nonlconExt20mm's cTendonMax = tendon - tendonMaxGeom
% <= 0 hard reject exactly (the draft's scalar compare never executed).
if isfield(info, 'tendonMax') && tendon > min(info.tendonMax) + 1e-9
    return                       % hard reject (nonlcon convention)
end
mL = routeLengths(Location, P.ctx);
N = numel(mL);
if strcmp(mode, 'aware')
    % AWARE = THROUGH THE CLASS OF RECORD. Three 1-D port attempts of
    % fortz/Lok (multi-segment span 71.6 -> class single-segment span +
    % scalar s0ref 0.0709 -> verbatim computeForceVector direction
    % 0.00995 max rel err) each failed the 1e-6 gate at the 25
    % collapsed-segment poses (76-100, knee -21.5..+10 deg, route rows
    % 1-5 coincident): the class's stored Funit there is row-unit of the
    % row-6 direction, which no reading of its normalize(F_vec) source
    % line reproduces (MATLAB normalize defaults to column z-score,
    % which would NaN the zero z-column). Rather than port drift, the
    % study evaluates the machinery itself; the port residual
    % (75/100 poses < 1e-6; max 0.99% at the 25 collapsed poses,
    % probe_gate2_20261027.m) is FLAGGED to Ben, not silently adopted.
    pulCfg = struct('nPulleyBPA', nB, 'tackleLineParts', G, 'gain', G, ...
        'routingMode', 'moving_via');
    obj = MonoPam_pulley('study1', Location, P.ctx.CrossPoint, P.Dia, ...
        P.ctx.T_Pam, rest, rest * (1 - P.KMAX), tendon, P.fitting, ...
        P.Pres, P.Xi(1), P.Xi(2), P.Xi(3), P.wraps, pulCfg);
    tension = double(obj.Ftendon(:));       % class convention: 0 at
    slack = logical(obj.PulleySlack(:));    % slack, FINITE at infeasible
    infeas = logical(obj.PulleyInfeasible(:));
    s = double(obj.sContraction(:));
    E.tensionGrid = tension;
    E.slackGrid = slack; E.infeasGrid = infeas;
    E.strainGrid = s / rest;
    E.ok = true;
else
    % UNAWARE = the same closure chi-zeroed (xi0 = 0, cb = 0, tendon
    % spring off), ported 1-D; MonoPam_pulley offers no kSpr-off switch,
    % so this mode has no class analog by construction. SPAN per the
    % class: single crossed segment, SCALAR s0ref = s0 at the
    % longest-span pose.
    xi0 = 0; cSpr = 0;
    span = spanSeg(Location, P.ctx);
    iRef = find(span == max(span), 1);
    DeltaL = span(iRef) - span;
    s0 = rest - mL + xi0 + tendon + 2 * P.fitting;
    s0ref = s0(iRef);
    Fm = P.FmaxOf(rest);
    Ffun = @(q) festo4(P.Dia, (q / rest) / P.KMAX, P.Pres) * Fm;
    sHi = P.KMAX * rest;
    cb = zeros(N, 1);
    s = zeros(N, 1); F = zeros(N, 1);
    slack = false(N, 1); infeas = false(N, 1);
    for i = 1:N
        psi = @(qq) G * qq - G * cb(i) * Ffun(qq) ...
            - nB * Ffun(qq) * cSpr / G - DeltaL(i) - G * s0ref;
        if psi(0) > 0
            s(i) = 0; F(i) = Ffun(0); slack(i) = true;
        elseif psi(sHi) < 0
            s(i) = sHi; F(i) = Ffun(sHi); infeas(i) = true;
        else
            try
                s(i) = max(fzero(psi, [0, sHi]), 0);
            catch
                s(i) = 0; slack(i) = true;  % base r=0 escape
            end
            F(i) = Ffun(s(i));
        end
    end
    tension = nB * F / G;
    tension(slack) = 0;
    tension(infeas) = NaN;
    E.tensionGrid = tension;
    E.slackGrid = slack; E.infeasGrid = infeas;
    E.strainGrid = s / rest;
    E.ok = true;
end
% Node quantities: the class leaves Ftendon FINITE at infeasible frames;
% the study's meet convention NaNs them (the torque-NaN contract).
infN = interp1(P.phiD, double(infeas), P.nodeDeg, 'linear', 'extrap') > 0.5;
E.tension = interp1(P.phiD, tension, P.nodeDeg, 'linear', 'extrap');
E.tension(infN) = NaN;
E.tau_a = E.tension .* P.armOsim;
E.tau_b = E.tension .* interp1(P.phiD, P.armCable, P.nodeDeg, 'linear', ...
    'extrap');
E.strain = interp1(P.phiD, E.strainGrid, P.nodeDeg, 'linear', 'extrap');
E.infeas = ~isfinite(E.tension);
E.slack = interp1(P.phiD, double(slack), P.nodeDeg, 'linear', 'extrap') > 0.5;
E.x = x; E.nB = nB; E.G = G; E.mode = mode;
end

%% --------------- search objective -------------------------------------
function f = objFun(x, nB, G, mode, meth, P)
E = evalConfig(x, nB, G, mode, P);
if ~E.ok
    f = 1e6; return
end
if strcmp(meth, 'a'), tq = E.tau_a; else, tq = E.tau_b; end
meet = isfinite(tq) & sign(tq) == sign(P.tauH) & ...
    abs(tq) >= abs(P.tauH) - 1e-12;
rat = tq ./ P.tauH;
rat(~isfinite(rat)) = -1;
f = -(1000 * mean(meet) + mean(min(rat, 1)));
end
