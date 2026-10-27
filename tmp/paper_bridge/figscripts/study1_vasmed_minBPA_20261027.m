% study1_vasmed_minBPA_20261027.m -- STUDY 1 (2026-10-27, easteregg2)
% Minimum 20 mm BPAs (parallel, up to 3 tackle line parts each) to
% meet/exceed the human vastus medialis SINGLE-muscle torque surface over
% the graded knee RoM, per moment-arm method.
%
% Ask (BiPulley studies directive, 2026-10-27):
%   - Human target: vas_med_r 5-node knee surface (1-D graded grid, the
%     act-11 p-node convention [-120 -87.5 -55 -22.5 10] deg), REGENERATED
%     2026-10-27 by gen_vas_med_surface_20261027.py (pinned opensim 4.4.1,
%     sha-pinned stock gait2392_simbody routing, computeMomentArm x Fmax
%     1294 N) -- re-run and verified reproduced before this script; read
%     from the CSV (nothing hardcoded).
%   - Moment arms BOTH ways, reported both:
%       (a) OpenSim ds/dtheta = arm_m column of that surface CSV;
%       (b) OUR cable arm: the path from the given muscle insertion point
%           (the cam design-of-record tibial insertion pEnd, xFullBest
%           rows 4:6) to the patellofibular insertion point (the OLD
%           PROXIMAL RING location from the cam context = ctxCamBest
%           routeSeed row 6, the stock pre-cam p6), treated as a cable
%           contacting the tibia surfaces 0.00075 m off them (task value;
%           the .m files' own bypass tolerance is 0.5 mm -- stated in the
%           report), moment about the robot knee ICR from the measured
%           t1-ICR table (fcn13/fcn14 of buildKneeExtContext20mm.m).
%   - Search: nBPA {1,2,3} x tackle G {1,2,3}; sizes [rest; tendon]
%     re-solved per configuration; STIFFNESS-UNAWARE first, then
%     STIFFNESS-AWARE. Tendon lower bound 0, upper the listed family
%     ceiling 0.3 m (also geometry-capped by routeInfo.tendonMax, hard
%     reject, the nonlconExt20mm convention). Rest upper-limited per the
%     standing rule (rest <= 0.90*Lchar; lower 0.30*Lchar, survey2 S2-5).
%   - Stiffness-aware chi = the ADOPTED EXTENSOR PICK of record
%     (minimizeExt10mmX3_results_20260920_flxr77.mat pick 32, carried by
%     the cam context as XiUsed): Xi0 = -0.0063780 m, Xi1 = 39978.7 N/m
%     (axial, along the BPA/femur line), Xi2 = 14733.8 N/m (normal,
%     bending). Xi1/Xi2 enter exactly as the pipeline's diagonal bracket
%     K = [X1 X2 X2] in the LokBi bracket frame (axis 1 rotated onto the
%     BPA line from the bracket centroid: chi1 along the line, chi2
%     normal, per the ruling). Chi3 (wrap-loss, 0.158) is NOT carried:
%     the pulley closure family (MonoPam_pulley / BiPam_pulley2_20261022)
%     has no Xi3 term -- stated, never silently dropped. Tendon spring
%     kSpr = wraps(3) * 1.51e-6 * 193e9 / tendon (extensor ctx of record
%     wraps = 3, fitting = 0.021 m). Pressure 620 kPa; KMAX = 0.255.
%
% Machinery (read-only, called, never edited): the cam result of record
%   Results/Vas_Pam_20mm_Result_cam_20261014_0416.mat (xFullBest = [p1,
%   pEnd, rest 0.8277, tendon 0.09249], xBestCam, ctxCamBest),
%   buildDistalRingLocation20mm, festo4, maxBPAforce, MonoPam_pulley
%   (regression gate only).
%
% Outputs (dated 20261027):
%   Results/BiPulley_surveys/VasMed_minBPA_study1_20261027.mat
%   tmp/paper_bridge/capability_table_20261022/VasMed_minBPA_study1_summary_20261027.csv
%
% EQUILIBRIUM GATE: the 1-D tackle closure below is a verbatim port of
% fortzBi2's per-cell solve (BiPam_pulley2_20261022.m). At (nB=1, G=1)
% with the adopted chi it must reproduce MonoPam_pulley's own Ftendon to
% 1e-6 relative (asserted before any search).

function study1_vasmed_main_20261027()

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
armCableNode = interp1(phiD, armCable, nodeDeg, 'linear');
fprintf('arm(b) cable (signed, nodes) [%s] m\n', num2str(armCableNode.', 5));
fprintf('  contact-active poses: upper %d/%d, lower %d/%d\n', ...
    nnz(contactLog(:, 1)), numel(phiD), nnz(contactLog(:, 2)), numel(phiD));

%% ---------------- Route cache + evaluator ------------------------------
routeCache = containers.Map('KeyType', 'double', 'ValueType', 'any');
[Location0, ~, info0] = getRoute(tendonRec); %#ok<ASGLU>
mL0 = routeLengths(Location0);
Lchar = max(mL0);
fprintf('Lchar (max rigid route length over grid) = %.6f m; rest bounds [%.4f %.4f]\n', ...
    Lchar, 0.30 * Lchar, 0.90 * Lchar);

FmaxOf = @(rest) maxBPAforce(rest, '20', Pres);

    function [Location, bendM, info] = getRoute(tendon)
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

    function mL = routeLengths(Location)
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

    function span = tendonSpan(Location)
        % Class convention: exit row = Cross-1 (row 5) to the insertion,
        % measured in the child (ICR) frame; row 5 mapped by T_Pam^-1.
        N = size(Location, 3);
        span = zeros(N, 1);
        for ii = 1:N
            T = ctx.T_Pam(:, :, ii);
            pPrev = RowVecTrans(inv(T), Location(ctx.CrossPoint - 1, :, ii));
            tot = 0;
            for j = ctx.CrossPoint:size(Location, 1) - 1
                tot = tot + norm(pPrev - Location(j, :, ii));
                pPrev = Location(j, :, ii);
            end
            span(ii) = tot;
        end
    end

    function cb = bracketCb(Location, xi1, xi2)
        % LokBi2's bracket projection, 1-D verbatim (Pbr 20 mm, RhbrZ*Ry
        % frame; cb = u' diag(1/[x1 x2 x2]) u per pose, u = BPA force
        % direction from the origin row).
        N = size(Location, 3);
        cb = zeros(N, 1);
        Pbr = [9.48, -36.21, 30.86] / 1000;
        phbrA = p1Rec - Pbr;
        th1 = atan2(phbrA(2), phbrA(1));
        RhbrZ = [cos(th1), -sin(th1), 0; sin(th1), cos(th1), 0; 0, 0, 1];
        pbrhA = RhbrZ' * phbrA.';
        th2 = atan2(pbrhA(3), pbrhA(1));
        Ry = [cos(th2), 0, sin(th2); 0, 1, 0; -sin(th2), 0, cos(th2)];
        Thbr = RpToTrans(RhbrZ * Ry', Pbr');
        Cb = diag(1 / [xi1, xi2, xi2]);
        for ii = 1:N
            u = Location(2, :, ii) - Location(1, :, ii);
            u = u / norm(u);
            ub = RowVecTrans(Thbr \ eye(4), u);
            cb(ii) = (ub * Cb) * ub.';
        end
    end

    function E = evalConfig(x, nB, G, mode)
        rest = x(1); tendon = x(2);
        E = struct('ok', false);
        if tendon < 0 || tendon > 0.3
            return
        end
        try
            [Location, ~, info] = getRoute(tendon);
        catch
            return
        end
        if isfield(info, 'tendonMax') && tendon > info.tendonMax + 1e-9
            return                       % hard reject (nonlcon convention)
        end
        mL = routeLengths(Location);
        span = tendonSpan(Location);
        N = numel(mL);
        DeltaL = max(span) - span;
        switch mode
            case 'unaware'
                xi0 = 0; cb = zeros(N, 1); cSpr = 0;
            case 'aware'
                xi0 = Xi(1); cb = bracketCb(Location, Xi(2), Xi(3));
                if tendon > 1e-9
                    cSpr = 1 / (wraps * 1.51e-6 * 193e9 / tendon);
                else
                    cSpr = 0;
                end
            otherwise
                error('mode')
        end
        s0 = rest - mL + xi0 + tendon + 2 * fitting;
        Fm = FmaxOf(rest);
        Ffun = @(s) festo4(Dia, (s / rest) / KMAX, Pres) * Fm;
        sHi = KMAX * rest;
        s = zeros(N, 1); F = zeros(N, 1);
        slack = false(N, 1); infeas = false(N, 1);
        for i = 1:N
            s0ref = s0(i) - DeltaL(i);
            psi = @(qq) G * qq - G * cb(i) * Ffun(qq) ...
                - nB * Ffun(qq) * cSpr / G - DeltaL(i) - G * s0ref;
            if psi(0) > 0
                s(i) = 0; F(i) = Ffun(0); slack(i) = true;
            elseif psi(sHi) < 0
                s(i) = sHi; F(i) = Ffun(sHi); infeas(i) = true;
            else
                s(i) = max(fzero(psi, [0, sHi]), 0);
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
        E.tension = interp1(phiD, tension, nodeDeg, 'linear');
        E.tau_a = E.tension .* armOsim;
        E.tau_b = E.tension .* interp1(phiD, armCable, nodeDeg, 'linear');
        E.strain = interp1(phiD, E.strainGrid, nodeDeg, 'linear');
        E.infeas = ~isfinite(E.tension);
        E.slack = interp1(phiD, double(slack), nodeDeg, 'linear') > 0.5;
        E.x = x; E.nB = nB; E.G = G; E.mode = mode;
    end

%% ---------------- Equilibrium gate vs MonoPam_pulley -------------------
gate = evalConfig([restRec, tendonRec], 1, 1, 'aware');
assert(gate.ok, 'gate eval failed');
pulCfg = struct('nPulleyBPA', 1, 'tackleLineParts', 1, 'gain', 1, ...
    'routingMode', 'moving_via');
obj = MonoPam_pulley('gate', getRoute(tendonRec), ctx.CrossPoint, Dia, ...
    ctx.T_Pam, restRec, restRec * (1 - KMAX), tendonRec, fitting, Pres, ...
    Xi(1), Xi(2), Xi(3), wraps, pulCfg);
ftClass = double(obj.Ftendon(:));
ftMine = gate.tensionGrid;
relErr = max(abs(ftMine - ftClass) ./ max(1, abs(ftClass)));
fprintf('EQUILIBRIUM GATE (1,1 aware): max rel err vs MonoPam_pulley = %.3g\n', relErr);
assert(relErr < 1e-6, 'closure port failed the MonoPam_pulley gate');
    end

%% ---------------- Search: 9 configs x 2 modes x 2 arm methods ----------
lb = [0.30 * Lchar; 0];
ub = [0.90 * Lchar; 0.3];
seeds = [0.60 * Lchar, 0.04; 0.30 * Lchar, 0.02; 0.90 * Lchar, 0.02; ...
    restRec, tendonRec];

    function f = objFun(x, nB, G, mode, meth)
        E = evalConfig(x, nB, G, mode);
        if ~E.ok
            f = 1e6; return
        end
        if strcmp(meth, 'a'), tq = E.tau_a; else, tq = E.tau_b; end
        meet = isfinite(tq) & sign(tq) == sign(tauH) & ...
            abs(tq) >= abs(tauH) - 1e-12;
        rat = tq ./ tauH;
        rat(~isfinite(rat)) = -1;
        f = -(1000 * mean(meet) + mean(min(rat, 1)));
    end

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
                    fSeed(sIdx) = objFun(seeds(sIdx, :), nB, G, md{1}, mh{1});
                end
                [~, ixSeed] = min(fSeed);
                [xB, fB] = patternsearch(@(x) objFun(x, nB, G, md{1}, mh{1}), ...
                    seeds(ixSeed, :), [], [], [], [], lb, ub, [], opts);
                E = evalConfig(xB, nB, G, md{1});
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

%% ---------------- Local helpers ----------------------------------------
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

function [u, nC] = cableFirstLeg(pFrom, pTo, circs)
% Unit direction of the FIRST leg of the taut cable pFrom -> pTo that
% contacts the clearance circles it cannot pass through. Iterative: take
% the first circle hit along the current leg, replace the leg with the
% tangent from the current point to that circle (tangent point chosen on
% the side the straight chord passes), repeat. nC = contacts used.
cur = pFrom;
dest = pTo;
used = false(1, numel(circs));
nC = zeros(1, numel(circs));
u = (dest - cur) / norm(dest - cur);
for it = 1:numel(circs)
    d = dest - cur;
    L = norm(d); u = d / L;
    hit = 0; hitT = inf;
    for k = 1:numel(circs)
        if used(k), continue; end
        rel = circs{k}.c - cur;
        proj = rel * u;
        if proj <= 0 || proj >= L, continue; end
        perp = norm(rel - proj * u);
        if perp < circs{k}.r && proj < hitT
            hitT = proj; hit = k;
        end
    end
    if hit == 0, break; end
    c = circs{hit}.c; r = circs{hit}.r;
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
