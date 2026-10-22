%TEST_HIPADDUCTIONENVELOPE_20261022 Gate tests, plan Item 5 verification gate.
%   Deterministic, pure-function asserts; depends only on files in this
%   folder (gate item d). Run from tmp/paper_bridge/opt_improve_20261022:
%     matlab -batch "test_hipAdductionEnvelope_20261022"
%   Gate items (plan lines 128-129):
%     (a) at the control points the interpolated limit equals the table
%         exactly;
%     (b) the full-flexion + full-abduction corner is masked (consistent
%         with buildRoutingContext.m:26's documented exclusion);
%     (c) with the PENDING table the untrusted-banner rule fires;
%     (d) zero dependence on any Ben file.
%   Extra checks: AddMag-band poses stay fully trusted; taper monotone and
%   inside the graded box.

here = fileparts(mfilename('fullpath'));
tblPath = fullfile(here, 'hipAdductionEnvelopeTable_20261022.csv');
T = readtable(tblPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'string');

opts = struct('taperBandDeg', 10, 'wMin', 0.1, 'verbose', false);
fprintf('--- test_hipAdductionEnvelope_20261022 (%d control points) ---\n', ...
    height(T));

% (a) interpolation exact at the control points
for i = 1:height(T)
    f = T.hipFlexDeg(i);
    [~, info] = hipAdductionEnvelope_20261022(f, 0, opts);
    assert(info.addLim == T.adducLimDeg(i), ...
        'addLim mismatch at control point %d (f=%g).', i, f);
    assert(info.abdLim == T.abdLimDeg(i), ...
        'abdLim mismatch at control point %d (f=%g).', i, f);
end
fprintf('PASS (a) interpolated limits equal the table at all %d control points\n', ...
    height(T));

% (b) full-flexion + full-abduction / full-adduction corners masked
WcornAbd = hipAdductionEnvelope_20261022(110, -45, opts);   % +abduction side
WcornAdd = hipAdductionEnvelope_20261022(110, 20, opts);    % +adduction side
assert(WcornAbd == opts.wMin, 'full-flexion+full-abduction corner must sit at the floor.');
assert(WcornAdd == opts.wMin, 'full-flexion+full-adduction corner must sit at the floor.');
Wband = hipAdductionEnvelope_20261022(92.5, -45, opts);     % just past band-edge limit
assert(Wband < 1 && Wband > opts.wMin, ...
    'pose just beyond the taper-edge limit must be partially down-weighted.');
assert(hipAdductionEnvelope_20261022(110, 0, opts) == 1, ...
    'pure full flexion (adduction neutral) stays fully trusted.');
fprintf(['PASS (b) corner masked: W(110,-45)=%g, W(110,20)=%g (== wMin); ' ...
    'W(92.5,-45)=%.4f (partial); W(110,0)=1\n'], WcornAbd, WcornAdd, Wband);

% Extra: Ben's AddMag pairing poses stay fully trusted
for f = [-25 -15 0 30 60 85]
    Wlo = hipAdductionEnvelope_20261022(f, -45, opts);
    Whi = hipAdductionEnvelope_20261022(f, 20, opts);
    assert(Wlo == 1 && Whi == 1, ...
        'AddMag-band pose (f=%g) with full adduction range must have W==1.', f);
end
fprintf('PASS (extra) full adduction range W==1 at f in [-25 -15 0 30 60 85] (AddMag band)\n');

% Extra: taper monotone toward the graded corner and inside the graded box
fgrid = 85:0.5:110;
[~, gi] = hipAdductionEnvelope_20261022(fgrid, zeros(size(fgrid)), opts);
assert(all(diff(gi.addLim) <= 0), 'addLim must be non-increasing toward full flexion.');
assert(all(diff(gi.abdLim) >= 0), 'abdLim must be non-decreasing toward full flexion.');
assert(all(gi.addLim >= 0 & gi.addLim <= 20 + eps), 'addLim outside graded box.');
assert(all(gi.abdLim <= 0 & gi.abdLim >= -45 - eps), 'abdLim outside graded box.');
fprintf('PASS (extra) taper monotone 85->110 deg and contained in the graded box [-45,20]\n');

% (c) PENDING table => untrusted banner rule fires
[~, info] = hipAdductionEnvelope_20261022(0, 0, opts);
assert(~info.trusted, 'PENDING-BEN rows present => trusted must be false.');
assert(~isempty(info.pendingRows), 'pendingRows must be non-empty.');
assert(startsWith(info.banner, 'HIP ENVELOPE UNTRUSTED'), ...
    'banner must start with the untrusted marker.');
assert(all(contains(T.source, 'PENDING-BEN')), ...
    'every shipped control point must carry the PENDING-BEN flag.');
fprintf('PASS (c) banner rule fires (%d/%d rows PENDING-BEN)\n', ...
    numel(info.pendingRows), height(T));
fprintf('drivers must print: %s\n', info.banner);

% (d) zero dependence on any Ben file
assert(~isempty(strfind(which('hipAdductionEnvelope_20261022'), ...
    'opt_improve_20261022')), 'helper must resolve inside the campaign folder.');
benFiles = {'Opt_run_pulley.m', 'Opt_run_BiPulley.m', 'MonoPam_pulley.m', ...
    'BiPam_pulley.m', 'biPulleySpecsFromOpenSim.m', 'Add_Mag_Mesh_Opt.m', ...
    'buildRoutingContext.m'};
onPath = benFiles(cellfun(@(s) ~isempty(which(s)), benFiles));
assert(isempty(onPath), 'Ben files must not be on the path for this test: %s.', ...
    strjoin(onPath, ', '));
fprintf('PASS (d) helper lives in the campaign folder; no Ben file on the path\n');

fprintf('TEST hipAdductionEnvelope_20261022: ALL PASS\n\n');
