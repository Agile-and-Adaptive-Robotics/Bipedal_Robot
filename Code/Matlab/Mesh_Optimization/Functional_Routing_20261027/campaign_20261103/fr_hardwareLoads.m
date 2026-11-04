function res = fr_hardwareLoads(varargin)
%FR_HARDWARELOADS Rated-load check against the S/M/L tackle families.
%
% Task 5 of the physical-packaging pass (directive item 7 / handoff
% "Load sizing that the optimizer must supply"). Two call forms:
%
%   res = fr_hardwareLoads(T_N, beta_rad, familyId)
%       single-cable tension T_N (N) turning through wrap angle beta_rad;
%       radial sheave reaction R = 2*T*sin(beta/2) (handoff formula; at
%       beta = 180 deg this is the R = 2*T reversing-sheave worst case).
%
%   res = fr_hardwareLoads(C, mech, cs)
%       walk a fr_mechBuild result: worst pose/crossing cable tension
%       (mech.Ftendon) and the geometric deviation angle at each active
%       crossing (same wrap-arc construction as fr_mechBuild's
%       fr_bendMeasure), plus the 180-deg-equivalent 2*T worst case that
%       fr_objective uses for its rated-load constraint.
%
% CAPACITIES (hardware_candidates.json / HARDWARE_AND_ZCODE_HANDOFF.md
% table, retrieved 2026-11-03 Pacific): S 400.34 N (9466T44), M 667.23 N
% (9466T41), L 3002.55 N (3434T78).
%
% CAVEAT attached to EVERY use (never stripped): the catalog capacities are
% published weight capacities whose interpretation for dynamic radial shaft
% loading is UNRESOLVED (handoff Open items; hardware_candidates.json
% capacity_basis). res.caveat carries that text and res.approved is always
% false: final hardware approval is NOT claimed by this function.
%
% Overload policy (task 5): demand > assigned-family capacity is FLAGGED in
% res.flagged (cell of strings) with the exact numbers, and res.ok = false.
% The demand is never resized, rerouted, or silently absorbed.
%
% Also returns res.packaging: the handoff's two-sheave stack heights
% (S 21.70 mm / M 24.875 mm / L 36.40 mm) and Bowden finished ODs
% (bicycle Jagwire KEB-SL 5.0 mm; Sava C277PL 7.747 mm finished OD -- NOT
% its bare casing OD) as packaging inputs, with the handoff citations.

res.caveat = ['UNRESOLVED convention: the catalog capacities are published ' ...
    'weight capacities; their interpretation for dynamic radial shaft ' ...
    'loading, cyclic duty, and cable compatibility is unconfirmed ' ...
    '(HARDWARE_AND_ZCODE_HANDOFF.md Open items). Catalog pulley diameter ' ...
    'is not assumed to be the wire pitch diameter.'];
res.approved = false;  %#ok<NASGU> (documented API field)
res.approvalNote = ['final hardware approval is NOT claimed: no assembly ' ...
    'WLL exists until cable, terminations, pins, bearings, plates, and ' ...
    'brackets are checked (handoff Open items).'];

switch nargin
    case 3
        if isstruct(varargin{1}) && isstruct(varargin{2})
            res = fr_walkMech(varargin{:});
            return
        end
        res = fr_single(varargin{1}, varargin{2}, varargin{3});
    otherwise
        error('fr_hardwareLoads:Args', ...
            'expected (T_N, beta_rad, familyId) or (C, mech, cs)');
end
end

% ==========================================================================
function res = fr_single(T_N, beta_rad, familyId)
T = double(T_N);
beta = double(beta_rad);
res = fr_shell(familyId);
res.tension_N = T;
res.wrapAngle_rad = beta;
res.reaction_N = 2 * T * abs(sin(beta / 2));
res.reversing2T_N = 2 * T;
% the 2T reversing-sheave worst case applies to THIS wrap only when the
% wrap itself is near 180 deg; the build-walk form keeps the unconditional
% fr_objective convention (worst case over the whole mechanism)
res = fr_judge(res, beta >= deg2rad(170));
end

function res = fr_shell(familyId)
valid = {'S', 'M', 'L'};
if ~any(strcmp(familyId, valid))
    error('fr_hardwareLoads:Family', ...
        'familyId must be S, M, or L (handoff families of record)');
end
cap = struct('S', 400.339945, 'M', 667.233242, 'L', 3002.549589);
part = struct('S', '9466T44', 'M', '9466T41', 'L', '3434T78');
od = struct('S', 31.75, 'M', 38.10, 'L', 63.50);        % mm
rope = struct('S', 1.5875, 'M', 2.38125, 'L', 3.175);   % mm
fam = familyId(1);
res.family = fam;
res.part = part.(fam);
res.capacity_N = cap.(fam);
res.sheaveOD_mm = od.(fam);
res.ropeDiameter_mm = rope.(fam);
res.flagged = {};
res.caveat = ['UNRESOLVED convention: the catalog capacities are ' ...
    'published weight capacities; their interpretation for dynamic ' ...
    'radial shaft loading, cyclic duty, and cable compatibility is ' ...
    'unconfirmed (HARDWARE_AND_ZCODE_HANDOFF.md Open items). Catalog ' ...
    'pulley diameter is not assumed to be the wire pitch diameter.'];
res.approved = false;
res.approvalNote = ['final hardware approval is NOT claimed: no ' ...
    'assembly WLL exists until cable, terminations, pins, bearings, ' ...
    'plates, and brackets are checked (handoff Open items).'];
res.packaging = struct( ...
    'source', ['HARDWARE_AND_ZCODE_HANDOFF.md (2026-11-03) + ' ...
    'hardware_candidates.json'], ...
    'stackTwoSheave_mm', struct('S', 21.70, 'M', 24.875, 'L', 36.40), ...
    'stackNote', ['two sheave widths + 1 mm inter-sheave gap + 1 mm ' ...
    'outer clearance + two plates (3/3/4 mm); axle head/nut and mounts ' ...
    'EXCLUDED; planar block allowance = sheave OD + 4 mm'], ...
    'bowdenFinishedOD_mm', struct('bicycle_KEB_SL', 5.0, ...
    'sava_C277PL', 7.747), ...
    'bowdenNote', ['Sava value is the FINISHED assembly OD (7.747 mm); ' ...
    'the bare 7.0358 mm casing OD must not be used (handoff). Bicycle ' ...
    'KEB-SL 5 mm OD is a dimensional candidate with no load rating.']);
res.ok = true;
end

function res = fr_judge(res, applyReversing)
if nargin < 2
    applyReversing = true;   % fr_objective build-walk convention
end
if res.reaction_N > res.capacity_N
    res.ok = false;
    res.flagged{end + 1} = sprintf( ...
        ['FLAGGED (never resized): family %s (%s) catalog capacity ' ...
        '%.1f N < radial reaction %.1f N (T = %.1f N, beta = %.1f deg, ' ...
        'R = 2*T*sin(beta/2)). Require a load-separated branch, a ' ...
        'different family, or a rerouting decision by Ben.'], ...
        res.family, res.part, res.capacity_N, res.reaction_N, ...
        res.tension_N, rad2deg(res.wrapAngle_rad));
end
if applyReversing && res.reversing2T_N > res.capacity_N
    res.ok = false;
    res.flagged{end + 1} = sprintf( ...
        ['FLAGGED (never resized): family %s 180-deg reversing-sheave ' ...
        'equivalent 2T = %.1f N exceeds catalog capacity %.1f N ' ...
        '(fr_objective rated-load convention).'], ...
        res.family, res.reversing2T_N, res.capacity_N);
end
end

% ==========================================================================
function res = fr_walkMech(~, mech, cs)
if ~cs.hasTackle
    res = fr_shell('M');   % base struct with caveat/packaging attached
    res.family = '(no tackle in this case)';
    res.capacity_N = NaN;
    res.tension_N = max(mech.Ftendon);
    res.wrapAngle_rad = NaN;
    res.reaction_N = NaN;
    res.reversing2T_N = 2 * res.tension_N;
    res.ok = true;
    res.flagged = {};
    res.note = 'case has no tackle; only the informational 2T record is kept';
    return
end
res = fr_shell(cs.txConfig.tackleFamily);
nG = mech.nG;
Tbest = 0;
betaBest = pi;   % worst-case default (180 deg) when no angle resolvable
for k = 1:numel(cs.crossFor{1})
    kk = min(k, numel(mech.cell.ActiveCrossing));
    cx = mech.caseDef.cross(kk);
    rIn = cx - 1;
    rOut = min(cx + 1, mech.caseDef.nRows);
    for i = 1:nG
        T1 = mech.T1list{i};
        T2 = mech.T2list{i};
        pIn = fr_apply(T1, mech.loc(rIn, :));
        pX = fr_apply(T2 * T1, mech.loc(cx, :));
        pOut = fr_apply(T2 * T1, mech.loc(rOut, :));
        u1 = pX - pIn;
        u2 = pOut - pX;
        if norm(u1) < 1e-9 || norm(u2) < 1e-9
            beta = pi;
        else
            beta = acos(min(1, max(-1, dot(u1 / norm(u1), ...
                u2 / norm(u2)))));
        end
        R = 2 * mech.Ftendon(i) * abs(sin(beta / 2));
        if R > 2 * Tbest * abs(sin(betaBest / 2)) - 1e-12 && ...
                mech.Ftendon(i) >= Tbest
            if R >= 2 * Tbest * abs(sin(betaBest / 2))
                Tbest = mech.Ftendon(i);
                betaBest = beta;
            end
        end
    end
end
res.tension_N = Tbest;
res.wrapAngle_rad = betaBest;
res.reaction_N = 2 * Tbest * abs(sin(betaBest / 2));
res.reversing2T_N = 2 * max(mech.Ftendon);
res.poseGrid = nG;
res = fr_judge(res);
res.note = sprintf( ...
    ['worst (T, beta) over %d poses x %d crossing(s) of the built mech; ' ...
    '2T uses the max pose tension (fr_objective convention)'], ...
    nG, numel(cs.crossFor{1}));
end

function p = fr_apply(Tm, prow)
p = (Tm(1:3, 1:3) * prow(:) + Tm(1:3, 4))';
end
