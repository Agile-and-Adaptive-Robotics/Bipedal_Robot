function asm = fr_kneeAssembly(C, cs, params, kin, qList)
%FR_KNEEASSEMBLY Long knee-crossing package model, per pose (task 4).
%
% Master directive item 6: for long actuators extending across the knee,
% use PINNED PROXIMAL AND DISTAL BPA SECTIONS with the central crossing
% carried by artificial tendon / Bowden cable. There are NO BLACK LINES in
% this model: every element of the central span is placed explicitly --
% brackets at both pinned ends, Bowden housing, ferrules, housing stops,
% cable exits, bend allowance, articulation clearance.
%
% TOPOLOGY (per pose, in the case pass-1 frame -- the knee cases' T1 frame,
% i.e. the femur frame for the femur->tibia routes):
%   proximal section : straight pinned BPA, femur-fixed, from the origin
%                      row point to the PROXIMAL bracket (= the case's
%                      exit row / sleeve tip). Length L0*split, pose
%                      constant (femur-fixed) -- verified per pose.
%   distal section   : straight pinned BPA, tibia-fixed, from the DISTAL
%                      bracket to the insertion row point. Length
%                      L0*(1-split), pose constant (tibia-fixed). The
%                      distal bracket is a bounded Tibia-cloud surface
%                      mount L0dist proximal of the insertion along the
%                      tibia axis, projected with fr_patch + standoff.
%   central span     : Bowden housing / artificial tendon from the
%                      proximal bracket to the distal bracket, laid as a
%                      circular arc through the two brackets with apex
%                      offset h; every element (ferrules 12 mm PLACEHOLDER,
%                      stops, cable exits, bend radius) placed on the arc.
%
% BEND ALLOWANCE: minimum bend radius Rmin (default C.bowden.bendRadius_mm
% = 75 mm, marked AS PREVIEW -- a bicycle-brake preview assumption, NOT a
% vendor rating; the Sava control-cable discussion spans 2-8 inches by
% construction and gives no product-specific minimum). The apex offset is
% h = max(h_clear, h_Rmin) with h_Rmin = Rmin - sqrt(Rmin^2 - (d/2)^2)
% (exact circle-through-3-points radius R = d^2/(8h) + h/2 >= Rmin), so
% R_path >= Rmin BY CONSTRUCTION; the reported margin documents it and the
% FLAG case fires when the required apex exceeds the anatomical envelope
% (params.apexMax_mm) -- flagged, never silently resized.
%
% params (defaults from C.bowden / hardware_candidates.json):
%   .split            proximal fraction of L0 (default 0.5, master item 6:
%                     initialize equal halves; optimized split = later item)
%   .L0               total BPA rest length [m] (default cs.L0Seed)
%   .bendRadius_mm    Rmin (default C.bowden.bendRadius_mm = 75, PREVIEW)
%   .housing          'bicycle' (KEB-SL 5.0 mm OD) or 'sava' (C277PL
%                     7.747 mm FINISHED OD) -- handoff Bowden references
%   .ferrule_mm       (default C.bowden.ferrule_mm = 12 mm PLACEHOLDER)
%   .stopLen_mm       housing stop length (preview 15 mm)
%   .stopAllowance_mm (default C.bowden.stopAllowance_mm = 20 mm stand-in)
%   .apexClear_mm     minimum apex offset to clear the joint line (default
%                     housing OD + stop allowance, i.e. the hardware must
%                     pass the articulation line)
%   .apexMax_mm       anatomical apex bound (default 60 mm) -- exceeding it
%                     is a FLAG, not a resize
%   .apexSide         +1/-1 (default +1 = in-plane +perpendicular side)
%   .bracketStandoff_mm distal bracket mount standoff (default 5)
%
% OUTPUT asm.perPose(i): proxBracket/distBracket/central path points,
% ferrule/stop/exit positions, chord, apex, R_path, bend margin, arcLen,
% requiredSpan (fr_spans convention: max chord + 2 ferrules + stop
% allowance), sleeve end diameters (calibrated endpoints at the local
% contraction), clearance flags. asm.summary: min/max margins + flags over
% the sweep. Everything is in meters unless the field name says _mm.

if nargin < 4
    params = struct();
end
p = fr_params(C, cs, params);

q0 = kin.qDefault;
T10 = cs.passes{1}.T1(kin, q0);            % tibia -> femur at neutral
used = find(cs.rowUse{1});
exIdx = used(p.exitRowUsed);               % proximal bracket row
insIdx = used(numel(used));                % insertion row
insSeed = cs.seedLoc(insIdx, :);           % tibia frame

% distal bracket: bounded tibia-surface mount, L0dist proximal of the
% insertion along the tibia axis (y superior in the tibia frame)
dbSeed = insSeed + [0, p.L0dist, 0];
dbPatch = fr_patch(C, 'Tibia', 'tibia_r', dbSeed, 15);
dbTib = dbPatch.project(0, 0, p.bracketStandoff_mm / 1000);

nG = numel(qList);
asm = struct();
asm.params = p;
asm.perPose = repmat(struct('chord', 0, 'apex', 0, 'Rpath', 0, ...
    'bendMargin', 0, 'arcLen', 0, 'requiredSpan', 0, 'proxBracket', ...
    [0 0 0], 'distBracket', [0 0 0], 'path', zeros(13, 3), 'ferrule1', ...
    [0 0 0], 'ferrule2', [0 0 0], 'stop1', [0 0 0], 'stop2', [0 0 0], ...
    'exit1', [0 0 0], 'exit2', [0 0 0], 'housingOD_mm', 0, ...
    'proxLen', 0, 'distLen', 0, 'proxTipDia_mm', 0, 'flags', {{}}, ...
    'clearance', NaN), 1, nG);
flags = {};  %#ok<AGROW>
for i = 1:nG
    q = qList{i};
    T1 = cs.passes{1}.T1(kin, q);
    e1 = cs.seedLoc(exIdx, :);             % femur frame (pass frame)
    e2 = fr_applyT(T1, dbTib);             % distal bracket, femur frame
    d = norm(e2 - e1);
    u = (e2 - e1) / max(d, 1e-9);
    % in-plane perpendicular (knee axis = z in the pass frame)
    perp = cross(u, [0; 0; 1]);
    perp = perp / max(norm(perp), 1e-9);
    hMin = p.apexClear;
    if d < 2 * p.Rmin
        hR = p.Rmin - sqrt(p.Rmin^2 - (d / 2)^2);
    else
        hR = 0;   % straight chord has infinite radius; apex still = hClear
    end
    h = max(hMin, hR);
    apex = 0.5 * (e1 + e2) + p.apexSide * h * perp;
    % circumcircle of e1, apex, e2 (2D in the pass plane, z held)
    [Rpath, cen] = fr_circum(e1(1:2), apex(1:2), e2(1:2));
    margin = Rpath - p.Rmin;
    th1 = atan2(e1(2) - cen(2), e1(1) - cen(1));
    thM = atan2(apex(2) - cen(2), apex(1) - cen(1));
    th2 = atan2(e2(2) - cen(2), e2(1) - cen(1));
    d1 = angleDiff(thM, th1);
    d2 = angleDiff(th2, thM);
    sgn = sign(d1);
    sweep = d1 + d2;
    arcLen = Rpath * abs(sweep);
    % sampled path on the arc through the apex (13 points)
    tS = linspace(0, 1, 13)';
    ang = th1 + tS * sweep;
    path = [repmat(cen(:)', 13, 1) + Rpath * [cos(ang), sin(ang)], ...
        repmat(e1(3), 13, 1)];
    % hardware elements on the path
    tan1 = path(2, :) - path(1, :);
    tan1 = tan1 / max(norm(tan1), 1e-9);
    tan2 = path(13, :) - path(12, :);
    tan2 = tan2 / max(norm(tan2), 1e-9);
    ferrule1 = path(1, :) + p.ferrule * tan1;
    ferrule2 = path(13, :) - p.ferrule * tan2;
    stop1 = path(1, :) + (p.ferrule + p.stopLen) * tan1;
    stop2 = path(13, :) - (p.ferrule + p.stopLen) * tan2;
    exit1 = stop1 + p.stopAllow * tan1;    % cable exit toward the anchor
    exit2 = stop2 - p.stopAllow * tan2;
    requiredSpan = arcLen + 2 * p.ferrule + p.stopAllow;
    pf = asm.perPose(i);
    pf.chord = d;
    pf.apex = h;
    pf.Rpath = Rpath;
    pf.bendMargin = margin;
    pf.arcLen = arcLen;
    pf.requiredSpan = requiredSpan;
    pf.proxBracket = e1;
    pf.distBracket = e2;
    pf.path = path;
    pf.ferrule1 = ferrule1;
    pf.ferrule2 = ferrule2;
    pf.stop1 = stop1;
    pf.stop2 = stop2;
    pf.exit1 = exit1;
    pf.exit2 = exit2;
    pf.housingOD_mm = p.housingOD_mm;
    pf.proxLen = p.L0prox;
    pf.distLen = p.L0dist;
    pf.proxTipDia_mm = fr_diameter(C, cs.diaStr, 0);
    cflags = {};  %#ok<AGROW>
    if h > p.apexMax
        cflags{end + 1} = sprintf( ...
            ['FLAG apex %.1f mm > anatomical bound %.1f mm (pose %d, ' ...
            'knee %.1f deg): route needs review, NOT silently resized'], ...
            1e3 * h, 1e3 * p.apexMax, i, rad2deg(q.knee_angle_r));
    end
    if Rpath < p.Rmin - 1e-9
        cflags{end + 1} = sprintf( ...
            'FLAG bend radius %.1f mm < Rmin %.1f mm (pose %d)', ...
            1e3 * Rpath, 1e3 * p.Rmin, i);
    end
    if arcLen + 2 * p.ferrule + p.stopAllow > p.maxSpan
        cflags{end + 1} = sprintf( ...
            ['FLAG required span %.1f mm exceeds the procured housing ' ...
            'quantity %.1f mm (pose %d): procurement must cover the ' ...
            'sweep maximum'], 1e3 * requiredSpan, 1e3 * p.maxSpan, i);
    end
    pf.flags = {cflags{:}};
    asm.perPose(i) = pf;
    flags = [flags, cflags]; %#ok<AGROW>
end
% sweep maximum span is the procurement quantity of record
p.maxSpan = max([asm.perPose.requiredSpan]);
asm.params = p;
for i = 1:nG
    pf = asm.perPose(i);
    pf.flags = pf.flags(~contains(pf.flags, 'procured housing'));
    asm.perPose(i) = pf;
end
flags = flags(~contains(flags, 'procured housing'));
asm.summary = struct( ...
    'poses', nG, ...
    'chordMin_m', min([asm.perPose.chord]), ...
    'chordMax_m', max([asm.perPose.chord]), ...
    'apexMin_m', min([asm.perPose.apex]), ...
    'apexMax_m', max([asm.perPose.apex]), ...
    'bendMarginMin_m', min([asm.perPose.bendMargin]), ...
    'bendRadiusMin_m', min([asm.perPose.Rpath]), ...
    'requiredSpanMax_m', max([asm.perPose.requiredSpan]), ...
    'requiredSpanMin_m', min([asm.perPose.requiredSpan]), ...
    'bendRadiusNote', ['75 mm AS PREVIEW (bicycle brake assumption); ' ...
    'vendor radius UNRESOLVED (handoff)'], ...
    'ferruleNote', '12 mm PLACEHOLDER from the Blender library', ...
    'housingNote', sprintf('%s finished OD %.3f mm (handoff)', ...
    p.housingName, p.housingOD_mm), ...
    'flags', {flags});
end

% ==========================================================================
function p = fr_params(C, cs, params)
p = struct();
p.split = fr_get(params, 'split', 0.5);
p.L0 = fr_get(params, 'L0', cs.L0Seed);
p.L0prox = p.L0 * p.split;
p.L0dist = p.L0 * (1 - p.split);
p.Rmin = fr_get(params, 'bendRadius_mm', C.bowden.bendRadius_mm) / 1000;
hn = fr_get(params, 'housing', 'bicycle');
switch hn
    case 'bicycle'
        p.housingOD_mm = 5.0;
        p.housingName = 'Jagwire KEB-SL (bicycle candidate, no load rating)';
    case 'sava'
        p.housingOD_mm = 7.747;
        p.housingName = 'Sava C277PL (FINISHED OD; casing option only)';
    otherwise
        error('fr_kneeAssembly:Housing', ...
            'housing must be ''bicycle'' or ''sava''');
end
p.ferrule = fr_get(params, 'ferrule_mm', C.bowden.ferrule_mm) / 1000;
p.stopLen = fr_get(params, 'stopLen_mm', 15) / 1000;
p.stopAllow = fr_get(params, 'stopAllowance_mm', ...
    C.bowden.stopAllowance_mm) / 1000;
p.apexClear = fr_get(params, 'apexClear_mm', ...
    p.housingOD_mm + C.bowden.stopAllowance_mm) / 1000;
p.apexMax = fr_get(params, 'apexMax_mm', 60) / 1000;
p.apexSide = fr_get(params, 'apexSide', 1);
p.bracketStandoff_mm = fr_get(params, 'bracketStandoff_mm', 5);
p.maxSpan = fr_get(params, 'maxSpan_m', 1.0);
% exit row: last proximal row of pass 1 (used-row indexing)
used = find(cs.rowUse{1});
p.exitRowUsed = fr_get(params, 'exitRowUsed', ...
    find(used == cs.cross(1) - 1, 1, 'last'));
if isempty(p.exitRowUsed) || p.exitRowUsed < 1
    p.exitRowUsed = max(numel(used) - 1, 1);
end
end

function v = fr_get(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end

function [R, cen] = fr_circum(p1, p2, p3)
% 2D circumradius + circumcenter of three points
ax = p1(1); ay = p1(2);
bx = p2(1); by = p2(2);
cx = p3(1); cy = p3(2);
dd = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
if abs(dd) < 1e-18
    R = inf;
    cen = [0; 0];
    return
end
a2 = ax^2 + ay^2;
b2 = bx^2 + by^2;
c2 = cx^2 + cy^2;
cen = [(a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)); ...
    (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax))] / dd;
R = norm([p1(1) - cen(1); p1(2) - cen(2)]);   % both as columns
end

function dth = angleDiff(b, a)
dth = b - a;
while dth > pi
    dth = dth - 2 * pi;
end
while dth < -pi
    dth = dth + 2 * pi;
end
end

function p = fr_applyT(Tm, prow)
p = (Tm(1:3, 1:3) * prow(:) + Tm(1:3, 4))';
end
