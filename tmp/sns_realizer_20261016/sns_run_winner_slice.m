function sns_run_winner_slice()
% Headless verification of SNS_WinnerSlice_R8.slx (2026-10-16 campaign).
%
% Four checks against NUMPY-VERIFIED targets. References (all in
% reports_20261016/sns_simulink/): the FULL winner network (real
% SNS_NumpyFixedTau backend, dt = 2 ms, same inputs - probes 2/3) anchors
% A/B/D (the RG/PF core is identical in slice and full net - nothing else
% feeds it under these inputs); the SLICE-topology reference (probe5 ->
% slice_ref_targets.json, exported subgraph + exact toolbox stepper)
% anchors the MN levels in C - in the full net every Renshaw cell is
% inhibited by ~44 same-side pools, so full-net MNs sit higher (2.902)
% than the 2-pool slice (2.247) in BOTH numpy and Simulink.
%
%   A  base config simulates 22 s and entrains: RG_E_r ~5 bursts in
%      6..15 s, period 1.7285 s, duty 0.315, vmax ~5.75 mV (RG_F antiphase
%      duty 0.630).
%   B  GAMMA hand-drive knob modulates the rhythm: KNOB_GAMMA_R x4 ->
%      duty ~0.64 and near-full-swing bursting (vmin < 1 mV).
%   C  ALPHA hand-drive knob modulates the MN pools: KNOB_ALPHA_R x0.5/x2
%      -> MN_vas_lat_r window mean ~2.055 / ~2.591 mV (RG unchanged).
%   D  heel/toe afferent resets set the cadence: contact script at 2.1 s
%      cycle -> period 2.1 s (~4 bursts), 1:1 entrainment.
%
% Tolerances are statistics-level: the circuit is chaotic (pointwise
% agreement between integrators dies at ~0.4 s - README_SNS_Simscape).

here = fileparts(mfilename('fullpath'));
root = fileparts(here);                    % ...\Code\Matlab\SNS_Simscape
repo = fileparts(fileparts(root));         % ...\Code
repo = fileparts(repo);                    % ...\Bipedal_Robot (repo root)
repOut = fullfile(repo, 'Code', 'MuJoCo_SNS', 'spinal', 'reports_20261016', ...
                   'sns_simulink');
if ~exist(repOut, 'dir'), mkdir(repOut); end

mdl = 'SNS_WinnerSlice_R8';
load_system(fullfile(root, 'SNS_Library_R2025a.slx'));
load_system(fullfile(here, [mdl '.slx']));
set_param(mdl, 'SignalLogging', 'on', 'SignalLoggingName', 'sigs');

logs = {'RG_E_r', 1, 'RG_E_r_V'; 'RG_F_r', 1, 'RG_F_r_V'; ...
        'RG_E_l', 1, 'RG_E_l_V'; 'RG_F_l', 1, 'RG_F_l_V'; ...
        'PF_E1_r', 1, 'PF_E1_r_V'; 'PF_F1_r', 1, 'PF_F1_r_V'; ...
        'MN_vas_lat_r', 1, 'MN_vas_V'; 'MN_bifemsh_r', 1, 'MN_bif_V'; ...
        'MN_vas_lat_r', 2, 'S_ext'; 'tq_sum', 1, 'tq_proxy'};
for k = 1:size(logs, 1)
    ph = get_param([mdl '/' logs{k, 1}], 'PortHandles');
    set_param(ph.Outport(logs{k, 2}), 'DataLogging', 'on', ...
        'DataLoggingNameMode', 'Custom', 'DataLoggingName', logs{k, 3});
end

KNOB0 = struct( ...
    'GAMMA_R', get_param([mdl '/KNOB_GAMMA_R'], 'Value'), ...
    'ALPHA_R', get_param([mdl '/KNOB_ALPHA_R'], 'Value'));
res = struct();
verdict = true;

% ---------------- CHECK A: base ------------------------------------------
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
get = @(nm) out.sigs.get(nm).Values;
R = grab(get, 'RG_E_r_V'); F = grab(get, 'RG_F_r_V');
sA = bstats(R.t, R.y, 6, 15);
fA = bstats(F.t, F.y, 6, 15);
fprintf(['CHECK A base: cycles=%d period=%s duty=%.3f v=[%.3f,%.3f] | ' ...
         'RG_F duty=%.3f v=[%.3f,%.3f]\n'], sA.cycles, ...
        fmtper(sA.period), sA.duty, sA.vmin, sA.vmax, fA.duty, ...
        fA.vmin, fA.vmax);
okA = sA.cycles >= 4 && abs(sA.period - 1.7285) < 0.09 && ...
      abs(sA.duty - 0.315) < 0.10 && abs(sA.vmax - 5.746) < 0.6 && ...
      abs(fA.duty - 0.630) < 0.12;
fprintf('CHECK A %s (targets: 5 cyc, 1.7285 s, duty 0.315, vmax 5.75)\n', ...
        passfail(okA));
verdict = verdict && okA;
res.base = struct('RG_E', sA, 'RG_F', fA);

% ---------------- CHECK B: gamma knob x4 ----------------------------------
set_param([mdl '/KNOB_GAMMA_R'], 'Value', ...
    num2str(4*str2double(KNOB0.GAMMA_R), 16));
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
get = @(nm) out.sigs.get(nm).Values;
R4 = grab(get, 'RG_E_r_V');
sB = bstats(R4.t, R4.y, 6, 15);
fprintf(['CHECK B gamma x4: cycles=%d period=%s duty=%.3f v=[%.3f,%.3f]' ...
         '\n'], sB.cycles, fmtper(sB.period), sB.duty, sB.vmin, sB.vmax);
okB = sB.cycles >= 4 && abs(sB.duty - 0.643) < 0.12 && sB.vmin < 1.0;
fprintf(['CHECK B %s (numpy target: duty 0.643, vmin -0.204 - the rhythm' ...
         ' is MODULATED by the GAMMA hand-drive knob)\n'], passfail(okB));
verdict = verdict && okB;
res.gamma_x4 = sB;
set_param([mdl '/KNOB_GAMMA_R'], 'Value', KNOB0.GAMMA_R);

% ---------------- CHECK C: alpha knob x0.5 / x2 ---------------------------
aVals = [0.5, 1, 2];
mnMean = zeros(3, 1);
for k = 1:3
    set_param([mdl '/KNOB_ALPHA_R'], 'Value', ...
        num2str(aVals(k)*str2double(KNOB0.ALPHA_R), 16));
    out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
    get = @(nm) out.sigs.get(nm).Values;
    M = grab(get, 'MN_vas_V');
    m = M.t >= 6 & M.t <= 15;
    mnMean(k) = mean(M.y(m));
end
set_param([mdl '/KNOB_ALPHA_R'], 'Value', KNOB0.ALPHA_R);
fprintf('CHECK C MN_vas mean: a0.5=%.3f a1=%.3f a2=%.3f mV\n', ...
        mnMean(1), mnMean(2), mnMean(3));
% slice-equivalent targets (probe5, slice_ref_targets.json): the FULL-net
% reference (probe3) reads 2.902/3.089/3.424 because in the full winner
% every Renshaw cell is inhibited by ~44 same-side pools (RC sits lower ->
% less recurrent inhibition -> higher MN); the 2-pool slice sits lower in
% BOTH numpy-slice and Simulink identically.
okC = abs(mnMean(1) - 2.055) < 0.15 && abs(mnMean(3) - 2.591) < 0.15 && ...
      mnMean(3) > mnMean(1);
fprintf(['CHECK C %s (slice-ref targets 2.055 / 2.591; full-net ref ' ...
         '2.902 / 3.424 for the 46-pool Renshaw network) - the ALPHA ' ...
         'knob sets MN-pool level without touching the rhythm\n'], ...
        passfail(okC));
verdict = verdict && okC;
res.mn_mean = mnMean; res.aVals = aVals;

% ---------------- CHECK D: contact cadence 2.1 s ---------------------------
set_cycle(mdl, 2.1);
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
get = @(nm) out.sigs.get(nm).Values;
R21 = grab(get, 'RG_E_r_V');
sD = bstats(R21.t, R21.y, 6, 15);
fprintf('CHECK D cycle 2.1: cycles=%d period=%s duty=%.3f\n', ...
        sD.cycles, fmtper(sD.period), sD.duty);
okD = sD.cycles >= 3 && abs(sD.period - 2.1) < 0.11;
fprintf(['CHECK D %s (heel/toe afferent resets entrain the RG 1:1 to the' ...
         ' contact cadence)\n'], passfail(okD));
verdict = verdict && okD;
res.cycle_21 = sD;
set_cycle(mdl, 1.7284);   % restore production cadence

% ---------------- figure + save -------------------------------------------
fig = figure('Visible', 'off', 'Position', [60 60 1100 820]);
tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(tl); hold on;
plot(R.t, R.y, 'LineWidth', 1.0); plot(F.t, F.y, 'LineWidth', 1.0);
grid on; ylabel('RG V (mV)'); legend('RG\_E\_r', 'RG\_F\_r');
title(sprintf(['Winner slice R8 base: entrained period %.4f s ' ...
    '(numpy ref 1.7285)'], sA.period));
nexttile(tl); hold on;
plot(R.t, R.y, 'LineWidth', 1.0); plot(R4.t, R4.y, 'LineWidth', 1.0);
grid on; ylabel('RG\_E\_r V (mV)'); legend('base', 'GAMMA\_R x4');
title('GAMMA hand-drive knob modulates the rhythm (duty/amplitude)');
nexttile(tl); hold on;
% base MN traces at restored knobs (final confirm run)
outB = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
g3 = @(nm) outB.sigs.get(nm).Values;
Mv = g3('MN_vas_V'); Mf = g3('MN_bif_V'); Tq = g3('tq_proxy');
mv = [Mv.Time, Mv.Data]; mf = [Mf.Time, Mf.Data];
tq = [Tq.Time, Tq.Data];
plot(mv(:,1), mv(:,2), 'LineWidth', 1.0);
plot(mf(:,1), mf(:,2), 'LineWidth', 1.0);
grid on; ylabel('MN V (mV)'); legend('vas\_lat (ext)', 'bifemsh (flx)');
title(sprintf(['MN levels (alpha x0.5/1/2 means %.3f/%.3f/%.3f mV vs ' ...
    'numpy 2.902/3.089/3.424)'], mnMean(1), mnMean(2), mnMean(3)));
nexttile(tl);
plot(tq(:,1), tq(:,2), 'LineWidth', 1.0);
grid on; ylabel('knee torque proxy (N m)'); xlabel('t (s)');
title('Demo load: 190*S_{ext} - 60*S_{flx} (nominal gains)');
exportgraphics(fig, fullfile(repOut, 'winner_slice_verify.png'), ...
    'Resolution', 130);
close(fig);

res.verdict = verdict;
res.targets = struct( ...
    'A', 'cycles 5, period 1.7285, duty 0.315, vmax 5.75 (RG_F 0.630)', ...
    'B', 'duty ~0.64 (slice ref 0.641), vmin ~ -0.23', ...
    'C', 'MN means 2.055 / 2.591 (slice ref; full-net 46-pool 2.902/3.424)', ...
    'D', 'period 2.1 at 2.1 s script');
res.timestamp = datetime('now');
save(fullfile(repOut, 'winner_slice_verify.mat'), '-struct', 'res');
close_system(mdl, 0);
close_system('SNS_Library_R2025a', 0);
if verdict
    fprintf('=== WINNER SLICE VERIFICATION: ALL CHECKS PASS ===\n');
else
    fprintf('=== WINNER SLICE VERIFICATION: FAILURES PRESENT ===\n');
end
end

% ---------------------------------------------------------------------------
function s = grab(get, nm)
v = get(nm);
s.t = v.Time(:); s.y = v.Data(:);
end

function s = bstats(t, y, t0, t1)
m = t >= t0 & t <= t1;
tt = t(m); yy = y(m);
thr = 0.5*(min(yy) + max(yy));
on = yy > thr;
rises = find(diff(double(on)) == 1);   % 0->1 transitions
if numel(rises) < 2
    s = struct('cycles', numel(rises), 'period', NaN, 'duty', mean(on), ...
               'vmin', min(yy), 'vmax', max(yy));
    return;
end
s = struct('cycles', numel(rises), ...
           'period', mean(diff(tt(rises))), 'duty', mean(on), ...
           'vmin', min(yy), 'vmax', max(yy));
end

function pf = passfail(ok)
if ok, pf = 'PASS'; else, pf = 'FAIL'; end
end

function p = fmtper(p)
if isnan(p), p = 'None'; else, p = sprintf('%.4f', p); end
end

function set_cycle(mdl, cycle)
% rewrite the contact-script Repeating Sequence vectors at a new cycle
[yh, yt, yl, tvec] = contact_vectors(cycle);
keys = {'heel_seq', 'toe_seq', 'load_seq'};
vals = {yh, yt, yl};
ymap = containers.Map(keys, vals);
tv = mat2str(round(tvec.'*1e9)./1e9, 12);
for sfx = {'contact_r', 'contact_l'}
    for k = 1:3
        b = [mdl '/' sfx{1} '/' keys{k}];
        set_param(b, 'rep_seq_t', tv, ...
                     'rep_seq_y', mat2str(round(vals{k}.'*1e9)./1e9, 12));
    end
end
end

function [yh, yt, yl, tvec] = contact_vectors(cycle)
% identical to sns_build_winner_slice.contact_vectors (kept in sync)
heel_rge = 1.4827569817170443;
toe_rge = 0.36000460790520544;
ib_rge = 0.31221919098930606;
onset = 0.8172727211325486;
phi = (0:0.005:1)';
fh = zeros(size(phi)); ft = zeros(size(phi));
fh(phi < 0.03) = 0.95*phi(phi < 0.03)/0.03;
m = phi >= 0.03 & phi < 0.30; fh(m) = 0.95 - 0.60*(phi(m)-0.03)/0.27;
m = phi >= 0.30 & phi < 0.55; fh(m) = 0.35 - 0.35*(phi(m)-0.30)/0.25;
m = phi >= 0.30 & phi < 0.55; ft(m) = 0.55*(phi(m)-0.30)/0.25;
m = phi >= 0.55 & phi < 0.62; ft(m) = 0.55;
m = phi >= 0.62 & phi < 0.70; ft(m) = 0.55 + 0.40*(phi(m)-0.62)/0.08;
m = phi >= 0.70 & phi < 0.72; ft(m) = 0.95 - 0.95*(phi(m)-0.70)/0.02;
fh = min(max(fh, 0), 1.5); ft = min(max(ft, 0), 1.5);
heel_sig = min(fh/0.35, 1.5);
toe_sig = min(ft/0.50, 1.5);
load_sig = min((fh + ft)/0.60, 1.5);
kick = exp(-(phi*cycle)/0.25);
kick(phi >= 0.72) = -exp(-((phi(phi >= 0.72)-0.72)*cycle)/0.25);
yh = heel_rge*heel_sig + onset*kick;
yt = toe_rge*toe_sig + onset*kick;
yl = ib_rge*load_sig;
tvec = phi*cycle;
end
