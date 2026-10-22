function [W, info] = hipAdductionEnvelope_20261022(hipFlexDeg, adductionDeg, opts)
%HIPADDUCTIONENVELOPE_20261022 Hip adduction/abduction limit vs hip flexion.
%   [W, INFO] = HIPADDUCTIONENVELOPE_20261022(HIPFLEXDEG, ADDUCTIONDEG, OPTS)
%   implements the hip adduction/abduction coupling envelope the pulley
%   sweep directive orders (PULLEY_SWEEP_DIRECTIVE_20261021.md §5,
%   "limits hip adduction/abduction based on hip flexion/extension").
%
%   Repo search verdict (this session, grep over *.m *.py *.md *.json in
%   Code/ tmp/ Documentation/): NO coupling code exists anywhere. The
%   nearest artifacts are the fixed task table entry 'hip_adduction_r',
%   [-25, 20] (Routing27_v2/buildRoutingContext.m:74-81), the pure
%   rectangular AddMag sweep (Add_Mag_Mesh_Opt.m:29-45), secondary
%   coordinates pinned at zero (biPulleySpecsFromOpenSim.m:25-26), and
%   MuJoCo capsule stiffness (spinal/runner.py:514 — not a kinematic
%   limit). So per the directive's "or write it", this file IS the
%   coupling code, with every control point FLAGGED for Ben's review.
%
%   Envelope definition (FLAGGED, PENDING-BEN):
%     adducLimDeg(f), abdLimDeg(f) = piecewise-linear interpolation over
%     the control points in hipAdductionEnvelopeTable_20261022.csv (same
%     folder). The shipped default table:
%       * box = the Gait2392 campaign graded ranges (actuator_map_27.json
%         rom_table_deg: hip_flexion [-15, 110] deg, hip_adduction
%         [-45, 20] deg) — the stock gait2392_simbody.osim itself carries
%         NO coordinate ranges (23 Coordinate elements, 0 range_min), so
%         the campaign tables are the repo's Gait2392 hip range data;
%       * full adduction range wherever Ben's own AddMag grid paired full
%         adduction with flexion (flexion -25..85 deg);
%       * quarter-ellipse taper sqrt(1-u^2), u = (f-85)/25, from the band
%         edge to the graded flexion limit, reaching 0 at the box corner
%         (110 deg) — the combined-extreme corner buildRoutingContext.m:26
%         documents as impractical.
%     Every row's source column carries PENDING-BEN: these are defensible
%     defaults, NOT physiology Ben has approved. Ben edits the CSV only;
%     this code never needs to change.
%
%   W (soft mask): 1 while abdLim(f) <= a <= adducLim(f); cosine taper to
%   OPTS.wMin across OPTS.taperBandDeg beyond the limit pair; wMin beyond
%   (humans CAN reach combined extremes — down-weight, never forbid;
%   kernel shared with softCornerWeights_20261022 per the plan).
%
%   INFO: struct with addLim/abdLim per pose, table, trusted flag, and
%   INFO.banner. While any table row is PENDING-BEN, INFO.trusted is
%   false and INFO.banner reads 'HIP ENVELOPE UNTRUSTED — coupling table
%   pending Ben'; every hip-site driver run MUST print it and treat
%   hip-site extreme corners with the Item 6 floor weight (plan Item 5).
%
%   Guardrails: read-only (no Ben file touched or required on the path);
%   weights only down-weight requirement terms — slack cells (F_ins
%   zeroed, MonoPam_pulley.m:595; BiPam_pulley.m:671) are excluded by
%   mask from torque tests and routed to the RoM/slack side, never
%   re-weighted (plan review 2026-10-22, Item 2 blocker).

if nargin < 3
    opts = struct();
end
band = getDef(opts, 'taperBandDeg', 10);   % PENDING-BEN (Item 6 corner band)
wMin = getDef(opts, 'wMin', 0.1);       % PENDING-BEN (Item 6 floor)
verbose = getDef(opts, 'verbose', false);
tblPath = getDef(opts, 'table', ...
    fullfile(fileparts(mfilename('fullpath')), ...
    'hipAdductionEnvelopeTable_20261022.csv'));  %#ok<NASGU>

T = readtable(tblPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'string');
tf = T.hipFlexDeg(:)';
ta = T.adducLimDeg(:)';
tb = T.abdLimDeg(:)';
if any(diff(tf) <= 0)
    error('hipAdductionEnvelope_20261022:table', ...
        'Table hipFlexDeg column must be strictly increasing.');
end
if any(ta < 0) || any(tb > 0) || any(ta < tb)
    error('hipAdductionEnvelope_20261022:table', ...
        'Need adducLimDeg >= 0 >= abdLimDeg with adducLimDeg >= abdLimDeg.');
end

F = double(hipFlexDeg(:));
A = double(adductionDeg(:));
if ~isequal(numel(F), numel(A))
    error('hipAdductionEnvelope_20261022:args', ...
        'hipFlexDeg and adductionDeg must have the same number of elements.');
end

addLim = interp1(tf, ta, F, 'linear');
abdLim = interp1(tf, tb, F, 'linear');
addLim(F < tf(1)) = ta(1);   % clamp outside the table span (no extrapolation)
addLim(F > tf(end)) = ta(end);
abdLim(F < tf(1)) = tb(1);
abdLim(F > tf(end)) = tb(end);

W = taperKernel_20261022(A, abdLim, addLim, band, wMin);
W = reshape(W, size(hipFlexDeg));

info = struct();
info.addLim = reshape(addLim, size(hipFlexDeg));
info.abdLim = reshape(abdLim, size(hipFlexDeg));
info.table = T;
info.taperBandDeg = band;
info.wMin = wMin;
info.pendingRows = find(contains(T.source, 'PENDING-BEN'));
info.trusted = isempty(info.pendingRows);
if info.trusted
    info.banner = '';
else
    info.banner = 'HIP ENVELOPE UNTRUSTED — coupling table pending Ben';
end
if verbose
    fprintf('%s\n', info.banner);
    fprintf('  control points (%d rows, %d PENDING-BEN):\n', ...
        height(T), numel(info.pendingRows));
    for i = 1:height(T)
        fprintf('    f=%7.2f deg -> adduction in [%8.4f, %8.4f] deg  (%s)\n', ...
            T.hipFlexDeg(i), T.abdLimDeg(i), T.adducLimDeg(i), T.source(i));
    end
end
end

function v = getDef(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
