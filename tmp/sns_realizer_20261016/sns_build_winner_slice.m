function sns_build_winner_slice()
% SNS_WINNER_SLICE - build SNS_WinnerSlice_R8.slx (demos\), a reduced but
% faithful Simulink realization of the 2026-10-16 walk-campaign WINNER
% (round8_best3.json = optuna curr_s8b_wean04 trial 35, drive_scale
% 0.3807379786697019) in the REDESIGNED SNS_Library architecture
% (2026-09-22: one-input synapses + INTERNAL neuronal summation).
%
% Slice (right leg detailed, left RG pair for the L/R drive demo):
%   - RG half-center pairs RG_E/RG_F per side = persistent-Na conditional
%     bursters (NapNeuron, FIXED tau_h 0.842 s - the toolbox
%     SNS_NumpyFixedTau semantics), IN-laminated mutual inhibition,
%     weak C1 flexor commissurals (winner no_cross: c1 0.1x, V3 absent)
%   - descending cells DRIVE/POSTURE + the campaign's four per-side HAND
%     DRIVE knobs GAMMA_l/r (-> own side's RG pair, descend_to_rg_e/f)
%     and ALPHA_l/r (-> ipsilateral MNs, alpha_to_mn), fed by the
%     production drive schedule (runner.drive_posture, walk drive
%     2.0773 nA) so u[GAMMA_c_s]=drive(t)*G[gamma_s] exactly as
%     runner.py:1372-1378
%   - one PF layer (PF_E1/E2/F1/F2 + PF_IN_E/F, Shevtsova lamination)
%   - heel/toe afferent reset INs + stance-Ib LBIN (winner heel_rge/
%     toe_rge/ib_rge/ib_e_central), driven by a stylized bilateral
%     antiphase foot-contact script at the production cadence 1.7284 s
%     (round8_full3.npz) shaped like the runner feeds it (incl. the
%     contact_onset edge kicks, 0.25 s decay)
%   - KINH swing quad suppression + IBEXC stance load sharing + Renshaw
%   - one knee MN pair (vas_lat_r ext / bifemsh_r flex) with the fitted
%     PF->MN / POSTURE->MN weights, driving a DEMO LOAD: net knee torque
%     proxy 190*S_ext - 60*S_flx N*m (nominal placeholder gains - NOT
%     tuned; the Xi-corrected model is the future replacement).
%
% All values come from slice_export.json (extracted from the REAL built
% winner network via the campaign's own loader - see
% spinal\_slice_probe_20261016.py); they are baked into block masks.
% RG somata use a purpose-built masked NapNeuron subsystem (SNS_Library
% has no persistent-Na block yet; candidate for promotion, Ben's call).
%
% Solver: fixed-step ode1 @ 2 ms = the production Euler semantics
% (runner dt 0.002; chaos note in README_SNS_Simscape - compare burst
% STATISTICS, not pointwise).
%
% Regenerate any time:  matlab -batch "run('...\sns_build_winner_slice.m')"

here = fileparts(mfilename('fullpath'));   % ...\Code\Matlab\SNS_Simscape\demos
root = fileparts(here);                    % ...\Code\Matlab\SNS_Simscape
repo = fileparts(fileparts(root));         % ...\Code
repo = fileparts(repo);                    % ...\Bipedal_Robot (repo root)

% ---- library (R2025a-native copy of the redesigned SNS_Library) ----
libFile = fullfile(root, 'SNS_Library_R2025a.slx');
assert(exist(libFile, 'file') == 2, 'SNS_Library_R2025a.slx not found');
load_system(libFile);
ph = get_param('SNS_Library_R2025a/NonSpikingNeuron', 'PortHandles');
assert(numel(ph.Inport) == 7 && numel(ph.Outport) == 2, ...
    ['SNS_Library_R2025a/NonSpikingNeuron does not have the redesigned ' ...
     'port signature (7 in / 2 out) - stale library?']);
phS = get_param('SNS_Library_R2025a/NonSpikingSynapse', 'PortHandles');
assert(numel(phS.Inport) == 1 && numel(phS.Outport) == 1, ...
    'SNS_Library_R2025a/NonSpikingSynapse is not one-input/one-output');
phY = get_param('SNS_Library_R2025a/SynSum', 'PortHandles');
assert(numel(phY.Inport) >= 2, 'SNS_Library_R2025a/SynSum missing');
LIBN = 'SNS_Library_R2025a';

% ---- slice export ----
J = jsondecode(fileread(fullfile(repo, 'Code', 'MuJoCo_SNS', 'spinal', ...
    'reports_20261016', 'sns_simulink', 'slice_export.json')));
sched = J.meta.schedule;
kn = sched.knob_G_values;
h0 = J.meta.h0;

% merge duplicate parallel synapses by SUMMING g (exact: the numpy
% backend accumulates parallel conductances; the 4x stacked
% IBEXC->LBIN edges = audit finding N1, kept effective value 4x0.5=2.0)
keys = cell(numel(J.synapses), 1);
for k = 1:numel(J.synapses)
    keys{k} = sprintf('%s|%s|%.12g|%.12g', J.synapses(k).src, ...
                      J.synapses(k).dst, J.synapses(k).Esyn_mV, ...
                      J.synapses(k).ThrPre_mV);
end
[uk, ia, ic] = unique(keys, 'stable');
syn = struct('src', cell(numel(uk), 1), 'dst', cell(numel(uk), 1), ...
             'g_uS', zeros(numel(uk), 1), 'Esyn_mV', zeros(numel(uk), 1), ...
             'ThrPre_mV', zeros(numel(uk), 1), ...
             'SlopePre_mV', zeros(numel(uk), 1));
for k = 1:numel(uk)
    m = (ic == k);
    syn(k).src = J.synapses(ia(k)).src;
    syn(k).dst = J.synapses(ia(k)).dst;
    syn(k).g_uS = sum([J.synapses(m).g_uS]);
    syn(k).Esyn_mV = J.synapses(ia(k)).Esyn_mV;
    syn(k).ThrPre_mV = J.synapses(ia(k)).ThrPre_mV;
    syn(k).SlopePre_mV = J.synapses(ia(k)).SlopePre_mV;
end
fprintf('slice: %d neurons, %d synapses (%d after merging parallel)\n', ...
        numel(J.neurons), numel(J.synapses), numel(syn));

% ---- model ----
mdl = 'SNS_WinnerSlice_R8';
if bdIsLoaded(mdl), close_system(mdl, 0); end
dst = fullfile(here, [mdl '.slx']);
if exist(dst, 'file'), delete(dst); end
new_system(mdl);
load_system(mdl);
set_param(mdl, 'UnconnectedInputMsg', 'none');

% ================= layout =================================================
% column x positions + per-column y walker (room for synapse stacks)
COL = struct('desc', 430, 'in', 620, 'rg', 800, 'sens', 990, ...
             'pf', 990, 'mn', 1290, 'out', 1560);
yCtr = struct('desc', 40, 'in', 40, 'rg', 40, 'sens', 40, 'pf', 660, ...
              'mn', 40, 'out', 40);
dY = 150;
colOf = containers.Map('KeyType', 'char', 'ValueType', 'char');
descCells = {'DRIVE', 'POSTURE', 'GAMMA_r', 'GAMMA_l', 'ALPHA_r', 'ALPHA_l'};
inCells = {'InE_r', 'InF_r', 'CIN_F_r', 'InE_l', 'InF_l', 'CIN_F_l'};
rgCells = {'RG_E_r', 'RG_F_r', 'RG_E_l', 'RG_F_l'};
sensCells = {'HEEL_r', 'TOE_r', 'LBIN_r', 'KINH_r', 'IBEXC_knee_ext_r', ...
             'HEEL_l', 'TOE_l', 'LBIN_l'};
pfCells = {'PF_E1_r', 'PF_E2_r', 'PF_F1_r', 'PF_F2_r', ...
           'PF_IN_E_r', 'PF_IN_F_r'};
mnCells = {'MN_vas_lat_r', 'RC_vas_lat_r', 'MN_bifemsh_r', 'RC_bifemsh_r'};
groups = {descCells, 'desc'; inCells, 'in'; rgCells, 'rg'; ...
          sensCells, 'sens'; pfCells, 'pf'; mnCells, 'mn'};
for g = 1:size(groups, 1)
    for nm = groups{g, 1}
        colOf(nm{1}) = groups{g, 2};
    end
end

ny = containers.Map('KeyType', 'char', 'ValueType', 'double');
npos = struct();
isCellN = iscell(J.neurons);
for k = 1:numel(J.neurons)
    if isCellN, nrec = J.neurons{k}; else, nrec = J.neurons(k); end
    nm = nrec.name;
    c = colOf(nm);
    y = yCtr.(c); yCtr.(c) = yCtr.(c) + dY; ny(nm) = y;
    if strcmp(nrec.kind, 'nap')
        build_nap(mdl, nm, [COL.(c) y COL.(c)+100 y+100], nrec, h0);
    else
        add_block([LIBN '/NonSpikingNeuron'], [mdl '/' nm], ...
            'Vrest', num2str(nrec.Vrest_mV, 12), ...
            'Gm', num2str(nrec.Gm_uS, 12), ...
            'Cm', num2str(nrec.Cm_nF, 12), ...
            'Thr', '0', 'Slope', '5', ...
            'Position', [COL.(c) y COL.(c)+80 y+80]);
    end
end

% ================= drive schedule + knobs ================================
% drive(t)   = D * f_up * (1 - f_dn),  f_up = clip((t-3)/2,0,1),
%                                        f_dn = clip((t-15)/2,0,1)
% posture(t) = 1 - 0.4*f_up + 0.4*f_dn     (runner.drive_posture, exact)
sx = 40; sy = 40;
add_block('simulink/Sources/Clock', [mdl '/clock1'], ...
    'Position', [sx sy sx+25 sy+25]);
add_block('simulink/Sources/Constant', [mdl '/DRIVE_MAX'], ...
    'Value', num2str(sched.walk_drive_nA, 16), ...
    'Position', [sx sy+40 sx+90 sy+70]);
mk = @(nm, pos, val) add_block('simulink/Sources/Constant', ...
    [mdl '/' nm], 'Value', num2str(val, 16), 'Position', pos);
mk('KNOB_GAMMA_L', [sx sy+90 sx+90 sy+120], kn.gamma_l);
mk('KNOB_GAMMA_R', [sx sy+130 sx+90 sy+160], kn.gamma_r);
mk('KNOB_ALPHA_L', [sx sy+170 sx+90 sy+200], kn.alpha_l);
mk('KNOB_ALPHA_R', [sx sy+210 sx+90 sy+240], kn.alpha_r);
% f_up / f_dn chains
fu = sched_clip(mdl, 'f_up', [sx+130 sy sx+40 sy+50], 3.0);
fd = sched_clip(mdl, 'f_dn', [sx+130 sy+70 sx+40 sy+120], 15.0);
% drive = D * f_up * (1-f_dn)
add_block('simulink/Math Operations/Gain', [mdl '/f_dn_neg'], ...
    'Gain', '-1', 'Position', [sx+200 sy+70 sx+230 sy+120]);
add_block('simulink/Math Operations/Sum', [mdl '/one_minus_fdn'], ...
    'Inputs', '++', 'Position', [sx+250 sy+70 sx+275 sy+120]);
add_block('simulink/Sources/Constant', [mdl '/one_c'], 'Value', '1', ...
    'Position', [sx+190 sy+130 sx+220 sy+160]);
add_block('simulink/Math Operations/Product', [mdl '/drive_prod'], ...
    'Inputs', '3', 'Position', [sx+310 sy+40 sx+340 sy+130]);
add_line(mdl, 'clock1/1', 'f_up_tsub/1', 'autorouting', 'on');
add_line(mdl, 'clock1/1', 'f_dn_tsub/1', 'autorouting', 'on');
add_line(mdl, 'f_dn/1', 'f_dn_neg/1', 'autorouting', 'on');
add_line(mdl, 'f_dn_neg/1', 'one_minus_fdn/1', 'autorouting', 'on');
add_line(mdl, 'one_c/1', 'one_minus_fdn/2', 'autorouting', 'on');
add_line(mdl, 'f_up/1', 'drive_prod/1', 'autorouting', 'on');
add_line(mdl, 'one_minus_fdn/1', 'drive_prod/2', 'autorouting', 'on');
add_line(mdl, 'DRIVE_MAX/1', 'drive_prod/3', 'autorouting', 'on');
% posture = 1 - 0.4 f_up + 0.4 f_dn
add_block('simulink/Math Operations/Gain', [mdl '/fu_g'], 'Gain', '-0.4', ...
    'Position', [sx+200 sy+180 sx+230 sy+210]);
add_block('simulink/Math Operations/Gain', [mdl '/fd_g'], 'Gain', '0.4', ...
    'Position', [sx+200 sy+220 sx+230 sy+250]);
add_block('simulink/Math Operations/Sum', [mdl '/posture_sum'], ...
    'Inputs', '+++', 'Position', [sx+270 sy+190 sx+295 sy+250]);
add_line(mdl, 'f_up/1', 'fu_g/1', 'autorouting', 'on');
add_line(mdl, 'f_dn/1', 'fd_g/1', 'autorouting', 'on');
add_line(mdl, 'one_c/1', 'posture_sum/1', 'autorouting', 'on');
add_line(mdl, 'fu_g/1', 'posture_sum/2', 'autorouting', 'on');
add_line(mdl, 'fd_g/1', 'posture_sum/3', 'autorouting', 'on');
% knob products: GAMMA_c_s = drive*K (ALPHA same)
kp = 1;
for knn = ["GAMMA_L", "GAMMA_R", "ALPHA_L", "ALPHA_R"]
    add_block('simulink/Math Operations/Product', ...
        sprintf('%s/knob_%s', mdl, knn), 'Position', ...
        [sx+370 sy+30*(kp-1) sx+400 sy+60+30*(kp-1)]);
    add_line(mdl, 'drive_prod/1', sprintf('knob_%s/1', knn), ...
        'autorouting', 'on');
    add_line(mdl, sprintf('KNOB_%s/1', knn), sprintf('knob_%s/2', knn), ...
        'autorouting', 'on');
    kp = kp + 1;
end

% wire the descending currents to the cells' Iapp
add_line(mdl, 'drive_prod/1', 'DRIVE/1', 'autorouting', 'on');
add_line(mdl, 'posture_sum/1', 'POSTURE/1', 'autorouting', 'on');
add_line(mdl, 'knob_GAMMA_R/1', 'GAMMA_r/1', 'autorouting', 'on');
add_line(mdl, 'knob_GAMMA_L/1', 'GAMMA_l/1', 'autorouting', 'on');
add_line(mdl, 'knob_ALPHA_R/1', 'ALPHA_r/1', 'autorouting', 'on');
add_line(mdl, 'knob_ALPHA_L/1', 'ALPHA_l/1', 'autorouting', 'on');

% ================= foot-contact script subsystems ========================
% Stylized bilateral antiphase stance pattern at the production cadence
% (1.7284 s), shaped exactly like the runner feeds the ports (incl.
% contact_onset edge kicks, 0.25 s decay). Signals are ZERO until t=5 s
% (walk onset). Vectors computed by local function contact_vectors().
CYCLE = 1.7284;
[cvh, cvt, cvl, cvtvec] = contact_vectors(CYCLE);
build_contact(mdl, 'contact_r', [40 1180 260 1420], 0.0, ...
              cvh, cvt, cvl, cvtvec);
build_contact(mdl, 'contact_l', [40 1450 260 1690], 0.5, ...
              cvh, cvt, cvl, cvtvec);
add_line(mdl, 'contact_r/1', 'HEEL_r/1', 'autorouting', 'on');
add_line(mdl, 'contact_r/2', 'TOE_r/1', 'autorouting', 'on');
add_line(mdl, 'contact_r/3', 'LBIN_r/1', 'autorouting', 'on');
add_line(mdl, 'contact_l/1', 'HEEL_l/1', 'autorouting', 'on');
add_line(mdl, 'contact_l/2', 'TOE_l/1', 'autorouting', 'on');
add_line(mdl, 'contact_l/3', 'LBIN_l/1', 'autorouting', 'on');

% MN standing-bias ports (production feeds 6*act_stand*stand_frac from the
% MuJoCo standing solve; out of slice -> 0 default, documented)
add_block('simulink/Sources/Constant', [mdl '/POST_vas_lat_r'], ...
    'Value', '0', 'Position', [1100 40 1190 70]);
add_block('simulink/Sources/Constant', [mdl '/POST_bifemsh_r'], ...
    'Value', '0', 'Position', [1100 90 1190 120]);
add_line(mdl, 'POST_vas_lat_r/1', 'MN_vas_lat_r/1', 'autorouting', 'on');
add_line(mdl, 'POST_bifemsh_r/1', 'MN_bifemsh_r/1', 'autorouting', 'on');

% ================= synapses ==============================================
% one one-input NonSpikingSynapse per (merged) edge, placed left of the
% destination neuron; >6 incoming -> chained SynSum into syn6 (the
% sns_build_from_json pattern)
incoming = struct();
for k = 1:numel(J.neurons)
    if isCellN, nmr = J.neurons{k}.name; else, nmr = J.neurons(k).name; end
    incoming.(nmr) = [];
end
for k = 1:numel(syn)
    incoming.(syn(k).dst) = [incoming.(syn(k).dst), k]; %#ok<AGROW>
end
NSYNPORT = 6;
for k = 1:numel(syn)
    d = syn(k).dst;
    c = colOf(d);
    j = find([incoming.(d)] == k);   % stack position of this synapse
    nm = sprintf('syn_%d_%s_to_%s', k, syn(k).src, d);
    add_block([LIBN '/NonSpikingSynapse'], [mdl '/' nm], ...
        'gmax', num2str(syn(k).g_uS, 12), ...
        'Esyn', num2str(syn(k).Esyn_mV, 12), ...
        'ThrPre', num2str(syn(k).ThrPre_mV, 12), ...
        'SlopePre', num2str(syn(k).SlopePre_mV, 12), ...
        'Position', [COL.(c)-150 ny(d)+10+26*(j-1) ...
                     COL.(c)-106 ny(d)+42+26*(j-1)]);
    set_param([mdl '/' nm], 'ShowName', 'off');
    add_line(mdl, [syn(k).src '/1'], [nm '/1'], 'autorouting', 'on');
end
for k = 1:numel(J.neurons)
    if isCellN, nm = J.neurons{k}.name; else, nm = J.neurons(k).name; end
    ks = incoming.(nm);
    ns = numel(ks);
    if ns == 0, continue; end
    if ns <= NSYNPORT
        nd = ns;
    else
        nd = 5;
    end
    for j = 1:nd
        add_line(mdl, sprintf('syn_%d_%s_to_%s/1', ks(j), ...
            syn(ks(j)).src, nm), sprintf('%s/%d', nm, j+1), ...
            'autorouting', 'on');
    end
    rest = ks(nd+1:end);
    if ~isempty(rest)
        prevSrc = '';
        i = 1; nss = 0;
        while i <= numel(rest)
            nss = nss + 1;
            ss = sprintf('ss_%d_%s', nss, nm);
            add_block([LIBN '/SynSum'], [mdl '/' ss], ...
                'Position', [COL.(colOf(nm))-240 ny(nm)+20+26*nss ...
                             COL.(colOf(nm))-196 ny(nm)+52+26*nss]);
            set_param([mdl '/' ss], 'ShowName', 'off');
            pn = 0;
            if ~isempty(prevSrc)
                pn = pn + 1;
                add_line(mdl, [prevSrc '/1'], sprintf('%s/%d', ss, pn), ...
                    'autorouting', 'on');
            end
            while pn < 8 && i <= numel(rest)
                pn = pn + 1;
                add_line(mdl, sprintf('syn_%d_%s_to_%s/1', rest(i), ...
                    syn(rest(i)).src, nm), sprintf('%s/%d', ss, pn), ...
                    'autorouting', 'on');
                i = i + 1;
            end
            prevSrc = ss;
        end
        add_line(mdl, [prevSrc '/1'], sprintf('%s/%d', nm, NSYNPORT+1), ...
            'autorouting', 'on');
    end
end

% ================= outputs: demo load + scopes ============================
% demo load: net knee-torque proxy with NOMINAL gains (placeholders, not
% tuned; swap for the Xi-corrected vasti/hamstring predictions later)
add_block('simulink/Math Operations/Gain', [mdl '/tq_ext'], ...
    'Gain', '190', 'Position', [1430 40 1460 80]);
add_block('simulink/Math Operations/Gain', [mdl '/tq_flx'], ...
    'Gain', '-60', 'Position', [1430 110 1460 150]);
add_block('simulink/Math Operations/Sum', [mdl '/tq_sum'], ...
    'Inputs', '++', 'Position', [1490 70 1515 120]);
add_line(mdl, 'MN_vas_lat_r/2', 'tq_ext/1', 'autorouting', 'on');
add_line(mdl, 'MN_bifemsh_r/2', 'tq_flx/1', 'autorouting', 'on');
add_line(mdl, 'tq_ext/1', 'tq_sum/1', 'autorouting', 'on');
add_line(mdl, 'tq_flx/1', 'tq_sum/2', 'autorouting', 'on');
outs = {'RG_E_r_V', 'RG_E_r/1'; 'RG_F_r_V', 'RG_F_r/1'; ...
        'S_ext', 'MN_vas_lat_r/2'; 'S_flx', 'MN_bifemsh_r/2'; ...
        'tq_proxy', 'tq_sum/1'};
for k = 1:size(outs, 1)
    add_block('simulink/Sinks/Out1', [mdl '/' outs{k, 1}], ...
        'Port', num2str(k), 'Position', [1560 40+60*(k-1) 1590 70+60*(k-1)]);
    add_line(mdl, outs{k, 2}, [outs{k, 1} '/1'], 'autorouting', 'on');
end
% scopes (GUI use; batch verification uses port DataLogging)
rgSrcs = {'RG_E_r/1', 'RG_F_r/1', 'RG_E_l/1', 'RG_F_l/1'};
add_block('simulink/Signal Routing/Mux', [mdl '/rg_mux'], 'Inputs', '4', ...
    'Position', [1430 200 1433 320]);
add_block('simulink/Sinks/Scope', [mdl '/RG_rhythm_scope'], ...
    'Position', [1480 230 1520 290]);
for j = 1:4
    add_line(mdl, rgSrcs{j}, sprintf('rg_mux/%d', j), ...
        'autorouting', 'on');
end
add_line(mdl, 'rg_mux/1', 'RG_rhythm_scope/1', 'autorouting', 'on');
add_block('simulink/Sinks/Scope', [mdl '/Torque_scope'], ...
    'Position', [1480 380 1520 430]);
add_line(mdl, 'tq_sum/1', 'Torque_scope/1', 'autorouting', 'on');

% ================= annotation + config + save =============================
try
    a = Simulink.Annotation(mdl, sprintf(['Winner circuit slice R8 (2026-10-16 campaign): ' ...
        'RG NaP half-centers + L/R GAMMA/ALPHA hand drive + PF layer + heel/toe resets ' ...
        '+ knee MN pair.\nSource: round8_best3.json (curr_s8b_wean04 t35, drive_scale 0.3807) ' ...
        'via slice_export.json - all values from the built winner network.\n' ...
        'ode1 @ 2 ms = production Euler semantics. Torque proxy gains 190/60 N*m are ' ...
        'NOMINAL placeholders.\nKNOB_* constants = per-side hand-drive knobs (verify: ' ...
        'gamma x4 flips RG duty; alpha shifts MN level).']));
    a.Position = [40 1720; 40 1720];
catch
end
set_param(mdl, 'Description', ['Reduced winner slice: round8_best3.json ' ...
    'human circuit in the redesigned SNS_Library architecture. See ' ...
    'reports_20261016/sns_simulink/.']);
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
    'FixedStep', '0.002', 'StopTime', '22');
save_system(mdl, dst);
fprintf('saved %s\n', dst);

% smoke: compile + 0.2 s
set_param(mdl, 'SignalLogging', 'off');
sim(mdl, 'StopTime', '0.2');
fprintf('smoke run 0.2 s OK\n');
close_system(mdl, 0);
close_system('SNS_Library_R2025a', 0);
end

% ===================== local functions ====================================
function out = sched_clip(mdl, pfx, pos, tEdge)
% clip((t - tEdge)/2, 0, 1): Sum(-) + Gain(1/2) + Saturation
add_block('simulink/Math Operations/Sum', [mdl '/' pfx '_tsub'], ...
    'Inputs', '+-', 'Position', [pos(1) pos(2) pos(1)+30 pos(2)+30]);
add_block('simulink/Sources/Constant', [mdl '/' pfx '_edge'], ...
    'Value', num2str(tEdge, 12), ...
    'Position', [pos(1) pos(2)+40 pos(1)+40 pos(2)+70]);
add_block('simulink/Math Operations/Gain', [mdl '/' pfx '_half'], ...
    'Gain', '0.5', 'Position', [pos(1)+50 pos(2) pos(1)+80 pos(2)+30]);
add_block('simulink/Discontinuities/Saturation', [mdl '/' pfx], ...
    'UpperLimit', '1', 'LowerLimit', '0', ...
    'Position', [pos(1)+100 pos(2) pos(1)+140 pos(2)+30]);
add_line(mdl, [pfx '_edge/1'], [pfx '_tsub/2'], 'autorouting', 'on');
add_line(mdl, [pfx '_tsub/1'], [pfx '_half/1'], 'autorouting', 'on');
add_line(mdl, [pfx '_half/1'], [pfx '/1'], 'autorouting', 'on');
% clock line added by the caller (tsub/1)
out = [pfx '/1'];
end

function build_contact(mdl, name, pos, lag, yh, yt, yl, tvec)
% One side's contact-script subsystem: 3 Repeating Sequence blocks
% (final port currents, incl. contact_onset kicks), gated to zero before
% t=5 s (walk onset). lag is informational (the vectors already carry the
% phase; contact_l gets the antiphase vectors).
blk = [mdl '/' name];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', pos);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/Clock', [blk '/clk'], ...
    'Position', [20 20 45 45]);
add_block('simulink/Sources/Constant', [blk '/t0'], 'Value', '5', ...
    'Position', [20 60 50 90]);
add_block('simulink/Logic and Bit Operations/Relational Operator', ...
    [blk '/ge5'], 'Operator', '>=', 'Position', [80 35 110 65]);
add_line(blk, 'clk/1', 'ge5/1', 'autorouting', 'on');
add_line(blk, 't0/1', 'ge5/2', 'autorouting', 'on');
sigs = {'heel_seq', yh, 'HEEL_c'; 'toe_seq', yt, 'TOE_c'; ...
        'load_seq', yl, 'LOAD_c'};
for k = 1:3
    s = sigs{k, 1}; y = sigs{k, 2};
    add_block('simulink/Sources/Repeating Sequence', [blk '/' s], ...
        'rep_seq_t', mat2str(round(tvec.'*1e9)./1e9, 12), ...
        'rep_seq_y', mat2str(round(y.'*1e9)./1e9, 12), ...
        'Position', [20 120+90*(k-1) 90 160+90*(k-1)]);
    add_block('simulink/Signal Routing/Switch', [blk '/' s '_gate'], ...
        'Position', [140 120+90*(k-1) 170 160+90*(k-1)]);
    add_block('simulink/Sources/Constant', [blk '/' s '_zero'], ...
        'Value', '0', 'Position', [90 200+90*(k-1) 120 230+90*(k-1)]);
    add_block('simulink/Sinks/Out1', [blk '/' sigs{k, 3}], ...
        'Port', num2str(k), 'Position', [210 130+90*(k-1) 240 150+90*(k-1)]);
    add_line(blk, [s '/1'], [s '_gate/1'], 'autorouting', 'on');
    add_line(blk, 'ge5/1', [s '_gate/2'], 'autorouting', 'on');
    add_line(blk, [s '_zero/1'], [s '_gate/3'], 'autorouting', 'on');
    add_line(blk, [s '_gate/1'], [sigs{k, 3} '/1'], 'autorouting', 'on');
end
% lag only affects annotation; vectors are pre-shifted by the caller
set_param(blk, 'Description', sprintf(['Foot-contact script (lag %.2f ' ...
    'cycle): stylized stance pattern, production cadence; runner-shaped ' ...
    'port currents incl. contact_onset kicks.'], lag));
end

function [yh, yt, yl, tvec] = contact_vectors(cycle)
% One cycle of the final port currents (nA) for ONE side - exactly the
% probe/runner shaping: heel_sig=min(Fh/0.35,1.5), toe_sig=min(Ft/0.50,
% 1.5), load_sig=min((Fh+Ft)/0.60,1.5); u_HEEL=heel_rge*sig+onset*kick
% (kick: +1 at loading edge phi=0, -1 at toe-off phi=0.72, decaying with
% 0.25 s). Winner gains: heel_rge 1.482757, toe_rge 0.360005, ib_rge
% 0.312219, contact_onset 0.817273 (from slice_export.json EFFECTIVE).
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

function build_nap(mdl, name, pos, nrec, h0)
% Persistent-Na RG soma in the REDESIGNED library conventions: port 1 =
% Iapp [nA], ports 2..7 = syn1..syn6 ([g; g*Esyn] lines, unconnected
% ground to 0), outs V [mV] + S(V). Dynamics = sns_toolbox
% NonSpikingNeuronWithPersistentSodiumChannel under SNS_NumpyFixedTau
% (CONSTANT tau_h - the production stepper):
%   Cm dV/dt = Gm(Vrest-V) + Iapp + (sum(gE) - V*sum(g))
%              + g_ion * a_inf(V) * h * (e_ion - V)
%   a_inf(V) = 1/(1+k_m*exp(slope_m*(e_m-V)))   [instantaneous]
%   dh/dt    = (b_inf(V) - h)/tau_h,  b_inf = 1/(1+k_h*exp(slope_h*(e_h-V)))
np = nrec.nap;
blk = [mdl '/' name];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', pos);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/Iapp'], 'Port', '1', ...
    'Position', [25 48 55 62]);
add_block('simulink/Math Operations/Sum', [blk '/IappSum'], ...
    'Inputs', '+', 'Position', [95 43 125 67]);
for k = 1:6
    add_block('simulink/Sources/In1', [blk '/syn' num2str(k)], ...
        'Port', num2str(k+1), 'PortDimensions', '2', ...
        'Position', [25 78+30*(k-1) 55 92+30*(k-1)]);
    add_block('simulink/Signal Routing/Demux', [blk '/D' num2str(k)], ...
        'Outputs', '2', 'Position', [95 72+30*(k-1) 98 106+30*(k-1)]);
    add_line(blk, ['syn' num2str(k) '/1'], ['D' num2str(k) '/1'], ...
        'autorouting', 'on');
end
add_block('simulink/Math Operations/Sum', [blk '/sumG'], ...
    'Inputs', '++++++', 'Position', [170 90 200 90+26*6]);
add_block('simulink/Math Operations/Sum', [blk '/sumE'], ...
    'Inputs', '++++++', 'Position', [170 90+26*6+40 200 90+52*6+40]);
for k = 1:6
    add_line(blk, ['D' num2str(k) '/1'], ['sumG/' num2str(k)], ...
        'autorouting', 'on');
    add_line(blk, ['D' num2str(k) '/2'], ['sumE/' num2str(k)], ...
        'autorouting', 'on');
end
add_block('simulink/Math Operations/Product', [blk '/gTimesV'], ...
    'Inputs', '2', 'Position', [270 95 300 125]);
add_block('simulink/Math Operations/Gain', [blk '/negGv'], 'Gain', '-1', ...
    'Position', [330 100 360 130]);
add_block('simulink/Math Operations/Sum', [blk '/Isyn'], 'Inputs', '++', ...
    'Position', [400 110 430 140]);
% NaP current
add_block('simulink/User-Defined Functions/Fcn', [blk '/aInf'], ...
    'Expression', '1/(1+kM*exp(sM*(eM-u(1))))', ...
    'Position', [270 180 340 210]);
add_block('simulink/User-Defined Functions/Fcn', [blk '/bInf'], ...
    'Expression', '1/(1+kH*exp(sH*(eH-u(1))))', ...
    'Position', [270 240 340 270]);
add_block('simulink/Math Operations/Sum', [blk '/dh'], 'Inputs', '+-', ...
    'Position', [400 240 430 270]);
add_block('simulink/Math Operations/Gain', [blk '/invTauH'], ...
    'Gain', '1000/tauH', 'Position', [470 240 500 270]);
add_block('simulink/Continuous/Integrator', [blk '/hInt'], ...
    'InitialCondition', 'h0', 'Position', [530 236 560 274]);
add_block('simulink/Math Operations/Sum', [blk '/naDrive'], ...
    'Inputs', '+-', 'Position', [400 300 430 330]);
add_block('simulink/Sources/Constant', [blk '/eion_c'], 'Value', 'eIon', ...
    'Position', [320 300 350 330]);
add_block('simulink/Math Operations/Product', [blk '/naProd'], ...
    'Inputs', '3', 'Position', [600 180 630 240]);
add_block('simulink/Math Operations/Gain', [blk '/gion_g'], ...
    'Gain', 'gIon', 'Position', [650 190 680 220]);
% membrane
add_block('simulink/Sources/Constant', [blk '/Vrest_c'], 'Value', 'Vrest', ...
    'Position', [330 360 360 390]);
add_block('simulink/Math Operations/Sum', [blk '/dVm'], 'Inputs', '-+', ...
    'Position', [400 355 430 385]);
add_block('simulink/Math Operations/Gain', [blk '/Gm'], 'Gain', 'Gm', ...
    'Position', [460 360 490 390]);
add_block('simulink/Math Operations/Sum', [blk '/dV'], 'Inputs', '++++', ...
    'Position', [740 140 770 220]);
add_block('simulink/Math Operations/Gain', [blk '/membrane'], ...
    'Gain', '1000/Cm', 'Position', [800 150 830 190]);
add_block('simulink/Continuous/Integrator', [blk '/membraneRC'], ...
    'InitialCondition', 'Vrest', 'Position', [870 146 900 184]);
% S output
add_block('simulink/Math Operations/Sum', [blk '/thrSub'], ...
    'Inputs', '+-', 'Position', [870 240 900 270]);
add_block('simulink/Sources/Constant', [blk '/Thr_c'], 'Value', 'Thr', ...
    'Position', [790 275 820 305]);
add_block('simulink/Math Operations/Gain', [blk '/invSlope'], ...
    'Gain', '1/Slope', 'Position', [940 240 970 270]);
add_block('simulink/Discontinuities/Saturation', [blk '/Sat'], ...
    'UpperLimit', '1', 'LowerLimit', '0', 'Position', [1000 240 1030 270]);
add_block('simulink/Sinks/Out1', [blk '/V_mV'], 'Port', '1', ...
    'Position', [1070 150 1100 164]);
add_block('simulink/Sinks/Out1', [blk '/S_drive'], 'Port', '2', ...
    'Position', [1070 240 1100 254]);
% wiring
add_line(blk, 'Iapp/1', 'IappSum/1', 'autorouting', 'on');
add_line(blk, 'IappSum/1', 'dV/1', 'autorouting', 'on');
add_line(blk, 'sumG/1', 'gTimesV/1', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'gTimesV/2', 'autorouting', 'on');
add_line(blk, 'gTimesV/1', 'negGv/1', 'autorouting', 'on');
add_line(blk, 'sumE/1', 'Isyn/1', 'autorouting', 'on');
add_line(blk, 'negGv/1', 'Isyn/2', 'autorouting', 'on');
add_line(blk, 'Isyn/1', 'dV/2', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'aInf/1', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'bInf/1', 'autorouting', 'on');
add_line(blk, 'bInf/1', 'dh/1', 'autorouting', 'on');
add_line(blk, 'hInt/1', 'dh/2', 'autorouting', 'on');
add_line(blk, 'dh/1', 'invTauH/1', 'autorouting', 'on');
add_line(blk, 'invTauH/1', 'hInt/1', 'autorouting', 'on');
add_line(blk, 'eion_c/1', 'naDrive/1', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'naDrive/2', 'autorouting', 'on');
add_line(blk, 'aInf/1', 'naProd/1', 'autorouting', 'on');
add_line(blk, 'hInt/1', 'naProd/2', 'autorouting', 'on');
add_line(blk, 'naDrive/1', 'naProd/3', 'autorouting', 'on');
add_line(blk, 'naProd/1', 'gion_g/1', 'autorouting', 'on');
add_line(blk, 'gion_g/1', 'dV/3', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'dVm/1', 'autorouting', 'on');
add_line(blk, 'Vrest_c/1', 'dVm/2', 'autorouting', 'on');
add_line(blk, 'dVm/1', 'Gm/1', 'autorouting', 'on');
add_line(blk, 'Gm/1', 'dV/4', 'autorouting', 'on');
add_line(blk, 'dV/1', 'membrane/1', 'autorouting', 'on');
add_line(blk, 'membrane/1', 'membraneRC/1', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'V_mV/1', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'thrSub/1', 'autorouting', 'on');
add_line(blk, 'Thr_c/1', 'thrSub/2', 'autorouting', 'on');
add_line(blk, 'thrSub/1', 'invSlope/1', 'autorouting', 'on');
add_line(blk, 'invSlope/1', 'Sat/1', 'autorouting', 'on');
add_line(blk, 'Sat/1', 'S_drive/1', 'autorouting', 'on');
% mask (SNS_Library conventions + the NaP channel parameters)
m = Simulink.Mask.create(blk);
m.Type = 'SNS NapNeuron (RG soma)';
m.Description = ['Persistent-Na non-spiking neuron, fixed tau_h ' ...
    '(sns_toolbox NonSpikingNeuronWithPersistentSodiumChannel under ' ...
    'SNS_NumpyFixedTau = the production RG stepper). Same port ' ...
    'conventions as SNS_Library/NonSpikingNeuron: Iapp [nA] + syn1..6, ' ...
    'outputs V [mV] and S(V).'];
m.Display = strjoin({ ...
    't_ = linspace(0, 2*pi, 73);' ...
    'patch(0.95*cos(t_), 0.95*sin(t_), [0 0 0]);' ...
    'patch(0.84*cos(t_), 0.84*sin(t_), [0.99 0.96 0.81]);' ...
    'x_ = linspace(-0.55, 0.55, 61);' ...
    'y_ = -0.20 + 0.40*exp(-((x_+0.08)/0.26).^2);' ...
    'color(''black'');' ...
    'plot(x_, y_);' ...
    'plot([-0.55 0.55], [-0.20 -0.20]);' ...
    'disp(''NaP'');'}, newline);
P = {'Vrest', 'Resting potential Vrest (mV)', ...
     num2str(nrec.Vrest_mV, 12); ...
     'Gm', 'Membrane conductance Gm (uS)', num2str(nrec.Gm_uS, 12); ...
     'Cm', 'Membrane capacitance Cm (nF)', num2str(nrec.Cm_nF, 12); ...
     'Thr', 'Saturation threshold Thr (mV)', '0'; ...
     'Slope', 'Saturation slope (1/mV)', '5'; ...
     'gIon', 'NaP conductance g_ion (uS)', num2str(np.g_ion, 12); ...
     'eIon', 'NaP reversal e_ion (mV)', num2str(np.e_ion, 12); ...
     'kM', 'm-gate k_m', num2str(np.k_m, 12); ...
     'sM', 'm-gate slope_m (1/mV)', num2str(np.slope_m, 12); ...
     'eM', 'm-gate midpoint e_m (mV)', num2str(np.e_m, 12); ...
     'kH', 'h-gate k_h', num2str(np.k_h, 12); ...
     'sH', 'h-gate slope_h (1/mV)', num2str(np.slope_h, 12); ...
     'eH', 'h-gate midpoint e_h (mV)', num2str(np.e_h, 12); ...
     'tauH', 'h-gate time constant tau_h (ms), FIXED', ...
     num2str(np.tau_h_s*1000, 12); ...
     'h0', 'initial h', num2str(h0, 12)};
for k = 1:size(P, 1)
    m.addParameter('Name', P{k, 1}, 'Type', 'edit', ...
        'Prompt', P{k, 2}, 'Value', P{k, 3});
end
set_param(blk, 'MaskIconFrame', 'off', 'MaskIconUnits', 'autoscale', ...
    'MaskIconOpaque', 'on', 'MaskIconRotate', 'none');
end
