function build_bpa_sns_rig_v2_20260926()
%% build_bpa_sns_rig_v2_20260926  Full BPA + SNS layer on the leg rig.
%% Port layouts (probed 2026-09-26):
%%   Transform Sensor: LConn(1)=base frame, RConn(1)=follower frame,
%%     RConn(2..5)=x,y,z,dist PS outputs   (base->follower, base axes)
%%   Internal Force:   LConn(1)=frame, LConn(2)=f PS input, RConn(1)=frame
%%   SNS blocks are plain Simulink (Iapp in, V out); BPA_20mm P[kPa],L[m]->F[N]

here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
addpath(sns);
mdl = 'mdl_leg_rig_ba003_imported';
if bdIsLoaded(mdl), close_system(mdl, 0); end
load_system(fullfile(sns, [mdl '.slx']));
if ~bdIsLoaded('SNS_Library'), load_system('SNS_Library'); end

V = @(n) [mdl '/' n];
B = @(src, n, pos) add_block(src, V(n), 'Position', pos);

% ---------------- branch source ports (frames) -------------------------------
kt = [mdl '/x04_01_KT_R_003_1_RIGID'];
kb = [mdl '/x04_02_KB_R_003_1_RIGID'];
ktPort = line_partner_port([mdl '/Revolute2'], kt);   % KT-side frame at knee
kbPort = line_partner_port([mdl '/Revolute1'], kb);   % KB-side frame at knee

% neuron defaults for the pressure map
vrest = str2double(get_param('SNS_Library/NonSpikingNeuron', 'Vrest'));
fprintf('neuron Vrest = %g mV\n', vrest);

muscles = struct('name', {'EXT', 'FLX'}, ...
    'dO', {{[0.0288 0.0447 0.0370], [-0.0212 0.0447 0.0370]}}, ...
    'dI', {{[0.0104 -0.0394 0.0380], [-0.0496 -0.0394 0.0380]}});

for m = 1:2
    pfx = muscles(m).name;
    % frames
    B('sm_lib/Frames and Transforms/Rigid Transform', [pfx '_origT'], [40 40*m 100 90*m]);
    set_param(V([pfx '_origT']), 'TranslationMethod', 'Cartesian', ...
        'TranslationCartesianOffset', mat2str(muscles(m).dO{1}, 6), 'RotationMethod', 'None');
    B('sm_lib/Frames and Transforms/Rigid Transform', [pfx '_insT'], [40 240+40*m 100 290+40*m]);
    set_param(V([pfx '_insT']), 'TranslationMethod', 'Cartesian', ...
        'TranslationCartesianOffset', mat2str(muscles(m).dI{1}, 6), 'RotationMethod', 'None');
    % sensor + force element
    B('sm_lib/Frames and Transforms/Transform Sensor', [pfx '_TS'], [160 40*m 220 120*m]);
    set_param(V([pfx '_TS']), 'SenseX', 'on', 'SenseY', 'on', 'SenseZ', 'on', 'SenseDist', 'on');
    B('sm_lib/Forces and Torques/Internal Force', [pfx '_IF'], [1020 40*m 1080 100*m]);
    % PS->Simulink converters for x,y,z,dist
    for c = 1:4
        B('nesl_utility/PS-Simulink Converter', [pfx '_c' num2str(c)], [260 40*m+30*c 310 50*m+30*c]);
    end
    % BPA force law
    B('SNS_Library/BPA_20mm', [pfx '_BPA'], [560 40*m 640 110*m]);
    set_param(V([pfx '_BPA']), 'Rest', '0.1301', 'Kmax', '0.1073');
    % force vector: Fx = F*x/dist etc
    B('simulink/Math Operations/Divide', [pfx '_invD'], [360 40*m 390 70*m]);
    set_param(V([pfx '_invD']), 'Inputs', '*/');
    B('simulink/Sources/Constant', [pfx '_one'], [300 40*m 330 60*m]);
    set_param(V([pfx '_one']), 'Value', '1');
    for c = 1:3
        B('simulink/Math Operations/Product', [pfx '_F' num2str(c)], [700 40*m+34*c 730 60*m+34*c]);
        set_param(V([pfx '_F' num2str(c)]), 'Inputs', '3');
    end
    B('simulink/Signal Routing/Mux', [pfx '_Mux'], [800 40*m 810 40*m+110]);
    set_param(V([pfx '_Mux']), 'Inputs', '3');
    B('nesl_utility/Simulink-PS Converter', [pfx '_sps'], [860 40*m 910 60*m]);
    % logging
    B('simulink/Sinks/To Workspace', [pfx '_logF'], [760 320+40*m 820 350+40*m]);
    set_param(V([pfx '_logF']), 'VariableName', ['F_' pfx], 'SaveFormat', 'Structure With Time');
    B('simulink/Sinks/To Workspace', [pfx '_logL'], [340 320+40*m 400 350+40*m]);
    set_param(V([pfx '_logL']), 'VariableName', ['L_' pfx], 'SaveFormat', 'Structure With Time');
end

% ---------------- wire frames + sensors + forces ------------------------------
for m = 1:2
    pfx = muscles(m).name;
    ot = get_param(V([pfx '_origT']), 'PortHandles'); otP = [ot.RConn ot.LConn];
    it = get_param(V([pfx '_insT']), 'PortHandles');  itP = [it.RConn it.LConn];
    ts = get_param(V([pfx '_TS']), 'PortHandles');    tsP = [ts.LConn ts.RConn];   % 1=base 2=follower 3..6=x y z dist
    ifh = get_param(V([pfx '_IF']), 'PortHandles');   ifP = [ifh.LConn ifh.RConn]; % 1=frame 2=PS f 3=frame
    add_line(mdl, ktPort, otP(1));
    add_line(mdl, kbPort, itP(1));
    add_line(mdl, otP(2), tsP(1));                    % base = origin frame
    add_line(mdl, itP(2), tsP(2));                    % follower = insertion frame
    add_line(mdl, otP(2), ifP(1));                    % force base  (branch)
    add_line(mdl, itP(2), ifP(3));                    % force follower (branch)
    % x,y,z,dist -> converters
    for c = 1:4
        cph = get_param(V([pfx '_c' num2str(c)]), 'PortHandles');
        add_line(mdl, tsP(2+c), cph.LConn(1));
    end
    % dist -> BPA.L ; dist -> invD
    c4 = get_param(V([pfx '_c4']), 'PortHandles');
    add_line(mdl, [pfx '_c4/1'], [pfx '_BPA/2'], 'autorouting', 'on');
    add_line(mdl, [pfx '_c4/1'], [pfx '_invD/2'], 'autorouting', 'on');
    add_line(mdl, [pfx '_one/1'], [pfx '_invD/1'], 'autorouting', 'on');
    % x,y,z -> F components
    for c = 1:3
        cc = get_param(V([pfx '_c' num2str(c)]), 'PortHandles');
        add_line(mdl, sprintf('%s_c%d/1', pfx, c), sprintf('%s_F%d/2', pfx, c), 'autorouting', 'on');
    end
    % F components -> Mux -> Simulink-PS -> IF.f
    for c = 1:3
        add_line(mdl, sprintf('%s_F%d/1', pfx, c), sprintf('%s_Mux/%d', pfx, c), 'autorouting', 'on');
    end
    mph = get_param(V([pfx '_Mux']), 'PortHandles');
    sph = get_param(V([pfx '_sps']), 'PortHandles');
    add_line(mdl, mph.Outport(1), sph.Inport(1));
    add_line(mdl, sph.RConn(1), ifP(2));
    % logging taps
    add_line(mdl, [pfx '_c4/1'], [pfx '_logL/1'], 'autorouting', 'on');
    add_line(mdl, [pfx '_BPA/1'], [pfx '_logF/1'], 'autorouting', 'on');
    fprintf('%s wired\n', pfx);
end

% ---------------- SNS antagonistic pair ---------------------------------------
B('SNS_Library/NonSpikingNeuron', 'N_E', [560 480 620 540]);
B('SNS_Library/NonSpikingNeuron', 'N_F', [560 600 620 660]);
B('SNS_Library/NonSpikingSynapse', 'S_EF', [660 560 700 590]);   % N_E inhibits N_F
B('SNS_Library/NonSpikingSynapse', 'S_FE', [660 620 700 650]);   % N_F inhibits N_E
for s = {'S_EF', 'S_FE'}
    try
        set_param(V(s{1}), 'gmax', '1.5', 'Esyn', num2str(vrest - 25));
    catch
        dp = fieldnames(get_param(V(s{1}), 'DialogParameters'));
        fprintf('synapse %s params: %s\n', s{1}, strjoin(dp', ' | '));
    end
end
% drives: antiphase sines (nA) into Iapp
B('simulink/Sources/Sine Wave', 'drv_E', [440 480 470 510]);
set_param(V('drv_E'), 'Amplitude', '3', 'Bias', '3.2', 'Frequency', '3.1416');  % 0.5 Hz
B('simulink/Sources/Sine Wave', 'drv_F', [440 600 470 630]);
set_param(V('drv_F'), 'Amplitude', '3', 'Bias', '3.2', 'Frequency', '3.1416', 'Phase', '3.1416');
B('simulink/Sinks/To Workspace', 'logV', [760 520 820 560]);
set_param(V('logV'), 'VariableName', 'V_E', 'SaveFormat', 'Structure With Time');

neh = get_param(V('N_E'), 'PortHandles');
nfh = get_param(V('N_F'), 'PortHandles');
deh = get_param(V('drv_E'), 'PortHandles');
dfh = get_param(V('drv_F'), 'PortHandles');
seh = get_param(V('S_EF'), 'PortHandles');
sfh = get_param(V('S_FE'), 'PortHandles');
% neuron ports: In1 = Iapp, In2.. = syn; Out1 = V  (verify by counts)
fprintf('N_E ports: in=%d out=%d | S_EF: in=%d out=%d\n', ...
    numel(neh.Inport), numel(neh.Outport), numel(seh.Inport), numel(seh.Outport));
add_line(mdl, deh.Outport(1), neh.Inport(1));
add_line(mdl, dfh.Outport(1), nfh.Inport(1));
add_line(mdl, neh.Outport(1), seh.Inport(1));      % V_E -> synapse pre
add_line(mdl, nfh.Outport(1), sfh.Inport(1));
add_line(mdl, seh.Outport(1), nfh.Inport(2));      % [g gE] -> N_F syn1
add_line(mdl, sfh.Outport(1), neh.Inport(2));      % [g gE] -> N_E syn1
lph = get_param(V('logV'), 'PortHandles');
add_line(mdl, neh.Outport(1), lph.Inport(1));

% pressure maps: P = 620 * sat((V - (Vrest+8))/25, 0..1)
for s = {'E', 'F'}
    B('simulink/Math Operations/Sum', ['dV_' s{1}], [860 480+120*(s{1}=='E') 890 510+120*(s{1}=='E')]);
    set_param(V(['dV_' s{1}]), 'Inputs', '+-');
    B('simulink/Sources/Constant', ['thr_' s{1}], [860 540+120*(s{1}=='E') 890 560+120*(s{1}=='E')]);
    set_param(V(['thr_' s{1}]), 'Value', num2str(vrest + 8));
    B('simulink/Discontinuities/Saturation', ['sat_' s{1}], [920 480+120*(s{1}=='E') 950 510+120*(s{1}=='E')]);
    set_param(V(['sat_' s{1}]), 'UpperLimit', '25', 'LowerLimit', '0');
    B('simulink/Math Operations/Gain', ['pmap_' s{1}], [980 480+120*(s{1}=='E') 1010 510+120*(s{1}=='E')]);
    set_param(V(['pmap_' s{1}]), 'Gain', '620/25');
    nx = ['N_' s{1}];
    nh = get_param(V(nx), 'PortHandles');
    dh = get_param(V(['dV_' s{1}]), 'PortHandles');
    th = get_param(V(['thr_' s{1}]), 'PortHandles');
    sh = get_param(V(['sat_' s{1}]), 'PortHandles');
    gh = get_param(V(['pmap_' s{1}]), 'PortHandles');
    add_line(mdl, nh.Outport(1), dh.Inport(1));
    add_line(mdl, th.Outport(1), dh.Inport(2));
    add_line(mdl, dh.Outport(1), sh.Inport(1));
    add_line(mdl, sh.Outport(1), gh.Inport(1));
    % P -> BPA.P of the matching muscle
    if strcmp(s{1}, 'E'), mp = 'EXT'; else, mp = 'FLX'; end
    add_line(mdl, ['pmap_' s{1} '/1'], [mp '_BPA/1'], 'autorouting', 'on');
end
fprintf('SNS wired\n');

% ---------------- compile + sim ------------------------------------------------
set_param(mdl, 'StopTime', '4');
try
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('UPDATE OK\n');
catch ME
    fprintf('UPDATE FAILED: %s\n', ME.message);
    for c = 1:min(numel(ME.cause), 8)
        fprintf('  CAUSE %d: %s\n', c, ME.cause{c}.message(1:min(end, 180)));
    end
    save_system(mdl);
    close_system(mdl, 0);
    return;
end
save_system(mdl);
out = sim(mdl);
L_E = out.L_EXT; L_F = out.L_FLX; F_E = out.F_EXT; F_F = out.F_FLX; VE = out.V_E;
fprintf('EXT length  %.4f -> %.4f -> %.4f m (t=0, 2, 4 s)\n', ...
    L_E.signals(1).values(1), interp1(L_E.time, L_E.signals(1).values, 2), ...
    interp1(L_E.time, L_E.signals(1).values, 4));
fprintf('FLX length  %.4f -> %.4f -> %.4f m\n', ...
    L_F.signals(1).values(1), interp1(L_F.time, L_F.signals(1).values, 2), ...
    interp1(L_F.time, L_F.signals(1).values, 4));
fprintf('EXT force   max %.1f N | FLX force max %.1f N | V_E range %.1f..%.1f mV\n', ...
    max(abs(F_E.signals(1).values)), max(abs(F_F.signals(1).values)), min(VE.signals(1).values), max(VE.signals(1).values));
save(fullfile(here, '..', 'results', 'rig_bpa_sns_20260926.mat'), 'out');
close_system(mdl, 0);
fprintf('=== build_bpa_sns_rig_v2 DONE ===\n');
end

function port = line_partner_port(jointBlk, bodyBlk)
ph = get_param(jointBlk, 'PortHandles');
ports = [ph.RConn ph.LConn];
hb = get_param(bodyBlk, 'Handle');
port = [];
for pp = 1:numel(ports)
    ln = get_param(ports(pp), 'Line');
    if ln < 0, continue; end
    phs = [];
    try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
    for q = 1:numel(phs)
        h = phs(q);
        if ~isscalar(h) || h <= 0 || h == ports(pp), continue; end
        try
            if get_param(get_param(h, 'Parent'), 'Handle') == hb
                port = ports(pp);
                return;
            end
        catch
        end
    end
end
assert(~isempty(port), 'no line from %s to %s', jointBlk, bodyBlk);
end
