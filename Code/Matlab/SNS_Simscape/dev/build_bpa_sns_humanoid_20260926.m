function build_bpa_sns_humanoid_20260926()
%% build_bpa_sns_humanoid_20260926  Replicate the rig recipe on both
%% humanoid knees: reduce each four-bar to a z hinge at the crank midpoint,
%% weld the brackets, add EXT/FLX BPA muscles + an SNS antagonistic pair.
%%
%% Per leg (S=R uses 04_01_KT_R/04_02_KB_R/04_05_BL_001_1/04_06_FL_001_1;
%% L uses 04_03_KT_L/04_04_KB_L/04_05_BL_001_2/04_06_FL_001_2):
%%   joints: RevoluteA(BL<->KB) weld -> BL rides with KB
%%           RevoluteB(FL<->KT) weld -> FL rides with KT
%%           Cylindrical(KB<->FL) delete bare
%%           RevoluteC(BL<->KT) delete bare
%%   hinge KT<->KB: KT-side anchor = old RevoluteB KT port (x-flipped source,
%%           counter-rotated), KB-side anchor = KB<->TI weld port (identity)
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
addpath(sns);
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));

legs = struct( ...
    'S', 'R', 'kt', 'x04_01_KT_R_001_1_RIGID', 'kb', 'x04_02_KB_R_001_1_RIGID', ...
    'ti', 'x05_01_TI_R_001_1_RIGID', 'bl', 'x04_05_BL_001_1_RIGID', ...
    'fl', 'x04_06_FL_001_1_RIGID', ...
    'jBLKB', 'Revolute', 'jFLKT', 'Revolute5', 'jKBFL', 'Cylindrical1', 'jBLKT', 'Revolute4');
legs(2) = struct( ...
    'S', 'L', 'kt', 'x04_03_KT_L_001_1_RIGID', 'kb', 'x04_04_KB_L_001_1_RIGID', ...
    'ti', 'x05_02_TI_001_1_RIGID', 'bl', 'x04_05_BL_001_2_RIGID', ...
    'fl', 'x04_06_FL_001_2_RIGID', ...
    'jBLKB', 'Revolute1', 'jFLKT', 'Revolute3', 'jKBFL', 'Cylindrical', 'jBLKT', 'Revolute2');

for lg = 1:2
    S = legs(lg).S;
    P = @(n) [sub '/' n];
    bodyOf = @(h) get_param(get_param(h, 'Parent'), 'Name');

    % ---- record far ports of the four-bar joints -------------------------
    jn = {legs(lg).jBLKB, legs(lg).jFLKT, legs(lg).jKBFL, legs(lg).jBLKT};
    far = struct();
    for w = 1:numel(jn)
        ph = get_param(P(jn{w}), 'PortHandles');
        ports = [ph.RConn ph.LConn];
        fh = zeros(0, 1);
        for pp = 1:numel(ports)
            fh = [fh; far_ports_of(ports(pp))]; %#ok<AGROW>
        end
        far.(jn{w}) = fh;
    end
    h_KB   = pick_single(far.(legs(lg).jBLKB), legs(lg).kb);
    h_BL   = pick_single(far.(legs(lg).jBLKB), legs(lg).bl);
    h_KTf  = pick_single(far.(legs(lg).jFLKT), legs(lg).kt);
    h_FL   = pick_single(far.(legs(lg).jFLKT), legs(lg).fl);
    h_KTb  = pick_single(far.(legs(lg).jBLKT), legs(lg).kt);

    % KB<->TI weld port: port 1 on both KB bodies (verified via plain-line map)
    phB = get_param(P(legs(lg).kb), 'PortHandles');
    bP = [phB.RConn phB.LConn];
    h_KBti = bP(1);

    % ---- surgery ---------------------------------------------------------
    delete_block(P(legs(lg).jBLKB));  add_line(sub, h_BL, h_KB);    % BL rides with KB
    delete_block(P(legs(lg).jFLKT));  add_line(sub, h_KTf, h_FL);   % FL rides with KT
    delete_block(P(legs(lg).jKBFL));                                % bare
    delete_block(P(legs(lg).jBLKT));                                % bare

    % hinge frames
    deal_frame(['pinT_KT_' S], h_KTf, [-0.0061 0.0077 0.0745], true, sub);
    deal_frame(['pinT_KB_' S], h_KBti, [0.0281 -0.0242 -0.0005], false, sub);
    add_block('sm_lib/Joints/Revolute Joint', P(['knee_hinge_' S]), 'Position', [200 560+200*lg 260 620+200*lg]);
    phK = get_param(P(['pinT_KT_' S]), 'PortHandles'); pK = [phK.RConn phK.LConn];
    phB = get_param(P(['pinT_KB_' S]), 'PortHandles'); pB = [phB.RConn phB.LConn];
    nh = get_param(P(['knee_hinge_' S]), 'PortHandles'); np = [nh.RConn nh.LConn];
    add_line(sub, pK(2), np(1));
    add_line(sub, pB(2), np(2));

    % ---- muscles ---------------------------------------------------------
    % EXT on +x side, FLX on -x side (pin-frame axes = part axes)
    deal_frame(['EXT_' S '_origT'], pK(2), [0.028 0.070 0.0], false, sub);
    deal_frame(['EXT_' S '_insT'],  pB(2), [0.033 -0.058 0.0], false, sub);
    deal_frame(['FLX_' S '_origT'], pK(2), [-0.028 0.070 0.0], false, sub);
    deal_frame(['FLX_' S '_insT'],  pB(2), [-0.033 -0.058 0.0], false, sub);

    ts = 'sm_lib/Frames and Transforms/Transform Sensor';
    for m = {'EXT', 'FLX'}
        pfx = [m{1} '_' S];
        add_block(ts, P([pfx '_TS']), 'Position', [160 40 220 120]);
        set_param(P([pfx '_TS']), 'SenseDist', 'on');
        ot = get_param(P([pfx '_origT']), 'PortHandles'); oP = [ot.RConn ot.LConn];
        it = get_param(P([pfx '_insT']), 'PortHandles'); iP = [it.RConn it.LConn];
        tsph = get_param(P([pfx '_TS']), 'PortHandles'); tsP = [tsph.LConn tsph.RConn];
        add_line(sub, oP(2), tsP(1));
        add_line(sub, iP(2), tsP(2));
        add_block('nesl_utility/PS-Simulink Converter', P([pfx '_c4']), 'Position', [260 40 310 60]);
        cph = get_param(P([pfx '_c4']), 'PortHandles');
        add_line(sub, tsP(3), cph.LConn(1));
        add_block('SNS_Library/BPA_20mm', P([pfx '_BPA']), 'Position', [560 40 640 110]);
        set_param(P([pfx '_BPA']), 'Rest', '0.1281', 'Kmax', '0.1057');
        add_line(sub, [pfx '_c4/1'], [pfx '_BPA/2'], 'autorouting', 'on');
        add_block('sm_lib/Forces and Torques/Internal Force', P([pfx '_IF']), 'Position', [1020 40 1080 100]);
        add_block('nesl_utility/Simulink-PS Converter', P([pfx '_sps']), 'Position', [860 40 910 60]);
        ifh = get_param(P([pfx '_IF']), 'PortHandles'); ifP = [ifh.LConn ifh.RConn];
        add_line(sub, oP(2), ifP(1));
        add_line(sub, iP(2), ifP(3));
        sph = get_param(P([pfx '_sps']), 'PortHandles');
        add_line(sub, [pfx '_BPA/1'], [pfx '_sps/1'], 'autorouting', 'on');
        add_line(sub, sph.RConn(1), ifP(2));
        add_block('simulink/Sinks/To Workspace', P([pfx '_logL']), 'Position', [340 320 400 350]);
        set_param(P([pfx '_logL']), 'VariableName', ['L_' pfx], 'SaveFormat', 'Structure With Time');
        add_line(sub, [pfx '_c4/1'], [pfx '_logL/1'], 'autorouting', 'on');
    end
    fprintf('leg %s: knee reduced + muscles built\n', S);
end

% ---- SNS: two antagonistic pairs (one per knee) --------------------------------
if ~bdIsLoaded('SNS_Library'), load_system('SNS_Library'); end
vrest = str2double(get_param('SNS_Library/NonSpikingNeuron', 'Vrest'));
for s = {'R', 'L'}
    B = @(src, n, pos) add_block(src, [sub '/' n], 'Position', pos);
    B('SNS_Library/NonSpikingNeuron', ['N_E_' s{1}], [560 480+200*(s{1}=='L') 620 540+200*(s{1}=='L')]);
    B('SNS_Library/NonSpikingNeuron', ['N_F_' s{1}], [560 600+200*(s{1}=='L') 620 660+200*(s{1}=='L')]);
    B('SNS_Library/NonSpikingSynapse', ['S_EF_' s{1}], [660 560+200*(s{1}=='L') 700 590+200*(s{1}=='L')]);
    B('SNS_Library/NonSpikingSynapse', ['S_FE_' s{1}], [660 620+200*(s{1}=='L') 700 650+200*(s{1}=='L')]);
    for sy = {'S_EF_', 'S_FE_'}
        set_param([sub '/' sy{1} s{1}], 'gmax', '1.5', 'Esyn', num2str(vrest - 25));
    end
    B('simulink/Sources/Sine Wave', ['drv_E_' s{1}], [440 480+200*(s{1}=='L') 470 510+200*(s{1}=='L')]);
    set_param([sub '/' ['drv_E_' s{1}]], 'Amplitude', '3', 'Bias', '3.2', 'Frequency', '3.1416');
    B('simulink/Sources/Sine Wave', ['drv_F_' s{1}], [440 600+200*(s{1}=='L') 470 630+200*(s{1}=='L')]);
    set_param([sub '/' ['drv_F_' s{1}]], 'Amplitude', '3', 'Bias', '3.2', 'Frequency', '3.1416', 'Phase', '3.1416');
    neh = get_param([sub '/' ['N_E_' s{1}]], 'PortHandles');
    nfh = get_param([sub '/' ['N_F_' s{1}]], 'PortHandles');
    deh = get_param([sub '/' ['drv_E_' s{1}]], 'PortHandles');
    dfh = get_param([sub '/' ['drv_F_' s{1}]], 'PortHandles');
    seh = get_param([sub '/' ['S_EF_' s{1}]], 'PortHandles');
    sfh = get_param([sub '/' ['S_FE_' s{1}]], 'PortHandles');
    add_line(sub, deh.Outport(1), neh.Inport(1));
    add_line(sub, dfh.Outport(1), nfh.Inport(1));
    add_line(sub, neh.Outport(1), seh.Inport(1));
    add_line(sub, nfh.Outport(1), sfh.Inport(1));
    add_line(sub, seh.Outport(1), nfh.Inport(2));
    add_line(sub, sfh.Outport(1), neh.Inport(2));
    % pressure maps
    for k = 1:2
        if k == 1, ss = 'E'; nn = ['N_E_' s{1}]; mp = ['EXT_' s{1}];
        else, ss = 'F'; nn = ['N_F_' s{1}]; mp = ['FLX_' s{1}]; end
        B('simulink/Math Operations/Sum', ['dV_' ss '_' s{1}], [860 480+120*k+200*(s{1}=='L') 890 510+120*k+200*(s{1}=='L')]);
        set_param([sub '/' ['dV_' ss '_' s{1}]], 'Inputs', '+-');
        B('simulink/Sources/Constant', ['thr_' ss '_' s{1}], [860 540+120*k+200*(s{1}=='L') 890 560+120*k+200*(s{1}=='L')]);
        set_param([sub '/' ['thr_' ss '_' s{1}]], 'Value', num2str(vrest + 8));
        B('simulink/Discontinuities/Saturation', ['sat_' ss '_' s{1}], [920 480+120*k+200*(s{1}=='L') 950 510+120*k+200*(s{1}=='L')]);
        set_param([sub '/' ['sat_' ss '_' s{1}]], 'UpperLimit', '25', 'LowerLimit', '0');
        B('simulink/Math Operations/Gain', ['pmap_' ss '_' s{1}], [980 480+120*k+200*(s{1}=='L') 1010 510+120*k+200*(s{1}=='L')]);
        set_param([sub '/' ['pmap_' ss '_' s{1}]], 'Gain', '620/25');
        nh2 = get_param([sub '/' nn], 'PortHandles');
        dh = get_param([sub '/' ['dV_' ss '_' s{1}]], 'PortHandles');
        th = get_param([sub '/' ['thr_' ss '_' s{1}]], 'PortHandles');
        sh = get_param([sub '/' ['sat_' ss '_' s{1}]], 'PortHandles');
        gh = get_param([sub '/' ['pmap_' ss '_' s{1}]], 'PortHandles');
        add_line(sub, nh2.Outport(1), dh.Inport(1));
        add_line(sub, th.Outport(1), dh.Inport(2));
        add_line(sub, dh.Outport(1), sh.Inport(1));
        add_line(sub, sh.Outport(1), gh.Inport(1));
        add_line(sub, ['pmap_' ss '_' s{1} '/1'], [mp '_BPA/1'], 'autorouting', 'on');
    end
end
fprintf('SNS pairs wired\n');

try
    set_param(mdl, 'StopTime', '3');
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('UPDATE OK\n');
    save_system(mdl);
    out = sim(mdl);
    LR = out.L_EXT_R.signals(1).values; LL = out.L_EXT_L.signals(1).values;
    fprintf('R knee EXT length: %.4f -> %.4f..%.4f m\n', LR(1), min(LR), max(LR));
    fprintf('L knee EXT length: %.4f -> %.4f..%.4f m\n', LL(1), min(LL), max(LL));
    save(fullfile(here, '..', 'results', 'humanoid_bpa_sns_20260926.mat'), 'out');
catch ME
    fprintf('UPDATE/SIM FAILED: %s\n', ME.message);
    for c = 1:min(numel(ME.cause), 8)
        fprintf('  CAUSE %d: %s\n', c, ME.cause{c}.message(1:min(end, 170)));
    end
end
save_system(mdl);
close_system(mdl, 0);
fprintf('=== build_bpa_sns_humanoid DONE ===\n');
end

function deal_frame(nm, srcPort, dT, counterRotate, sub)
add_block('sm_lib/Frames and Transforms/Rigid Transform', [sub '/' nm], 'Position', [40 700 100 750]);
if counterRotate
    set_param([sub '/' nm], 'TranslationMethod', 'Cartesian', ...
        'TranslationCartesianOffset', mat2str(dT, 6), 'RotationMethod', 'ArbitraryAxis', ...
        'RotationAngleUnits', 'rad', 'RotationAngle', '3.14159265358979', 'RotationArbitraryAxis', '[1 0 0]');
else
    set_param([sub '/' nm], 'TranslationMethod', 'Cartesian', ...
        'TranslationCartesianOffset', mat2str(dT, 6), 'RotationMethod', 'None');
end
ph = get_param([sub '/' nm], 'PortHandles');
add_line(sub, srcPort, ph.RConn(1));
end

function fp = far_ports_of(porth)
fp = zeros(0, 1);
ln = get_param(porth, 'Line');
if ln < 0, return; end
phs = [];
try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
br = -1;
try, br = get_param(porth, 'Branch'); catch, end
if isscalar(br) && br > 0
    try, phs = [phs get_param(br, 'BranchHandles')]; catch, end
end
seen = [];
for q = 1:numel(phs)
    h = phs(q);
    if ~isscalar(h) || h <= 0 || h == porth, continue; end
    if any(seen == h), continue; end
    seen(end+1) = h; %#ok<AGROW>
    try
        pn = get_param(get_param(h, 'Parent'), 'Name');
    catch
        continue;
    end
    if ~isempty(regexpi(pn, '_RIGID$', 'once'))
        fp(end+1, 1) = h; %#ok<AGROW>
    end
end
end

function h = pick_single(fh, bodyA)
h = [];
for k = 1:size(fh, 1)
    pn = get_param(get_param(fh(k, 1), 'Parent'), 'Name');
    if strcmp(pn, bodyA), h = fh(k, 1); end
end
assert(~isempty(h), 'pick_single failed for %s', bodyA);
end
