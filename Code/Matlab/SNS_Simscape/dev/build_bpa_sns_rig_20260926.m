function build_bpa_sns_rig_20260926()
%% build_bpa_sns_rig_20260926  Add BPA actuation + SNS control to the leg rig:
%%   - EXT + FLX point-muscles (Internal Force elements) between KT (thigh,
%%     merged femur) and KB (shank, merged tibia), force law = SNS_Library
%%     BPA_20mm (Ben's Festo equations), length from a Transform Sensor
%%   - SNS: NonSpikingNeuron antagonistic pair with mutual inhibition,
%%     sine drive, membrane-voltage -> pressure map
%%   - log muscle length + force + pressures, 4 s sim, save model + results

here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
addpath(sns);                       % SNS_Library.slx on path
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

% ---------- discover SNS library block names --------------------------------
if ~bdIsLoaded('SNS_Library'), load_system('SNS_Library'); end
libBlks = find_system('SNS_Library', 'SearchDepth', 1, 'Type', 'Block');
neuBlk = ''; synBlk = ''; bpaBlk = '';
for k = 1:numel(libBlks)
    nm = get_param(libBlks{k}, 'Name');
    if strcmpi(nm, 'NonSpikingNeuron'), neuBlk = libBlks{k}; end
    if strcmpi(nm, 'NonSpikingSynapse'), synBlk = libBlks{k}; end
    if strcmpi(nm, 'BPA_20mm'), bpaBlk = libBlks{k}; end
end
assert(~isempty(neuBlk), 'NonSpikingNeuron not found in SNS_Library');
assert(~isempty(bpaBlk), 'BPA_20mm not found in SNS_Library');
fprintf('SNS blocks: neuron=%s synapse=%s bpa=%s\n', neuBlk, synBlk, bpaBlk);

% neuron defaults (for the pressure map)
nv = fieldnames(get_param(neuBlk, 'DialogParameters'));
fprintf('neuron params: %s\n', strjoin(nv', ', '));

% ---------- frame branch sources ---------------------------------------------
kt = [mdl '/x04_01_KT_R_003_1_RIGID'];   % thigh (femur merged)
kb = [mdl '/x04_02_KB_R_003_1_RIGID'];   % shank (tibia merged)

% muscle attachment points, host part coords (m), identity part rotations:
% EXT femur (0.0207, 0.0593, 0.0000) | EXT shank (0.0343, -0.0293, 0.0001)
% FLX femur (-0.0293, 0.0593, 0.0000) | FLX shank (-0.0257, -0.0293, 0.0001)
% branch frames (smiData, identity rotation):
%   KT-side entry25: (-0.00806, 0.01461, -0.0370) | KB-side entry11: (0.02389, 0.01014, -0.03787)
mkMuscle = @(nm, ptKT, ptKB) struct('name', nm, 'ptKT', ptKT, 'ptKB', ptKB);
muscles = [ ...
    mkMuscle('EXT', [0.0207 0.0593 0.0000], [0.0343 -0.0293 0.0001]); ...
    mkMuscle('FLX', [-0.0293 0.0593 0.0000], [-0.0257 -0.0293 0.0001])];
fKT = [-0.00806 0.01461 -0.0370];
fKB = [0.02389 0.01014 -0.03787];

% ---------- find the Revolute lines to branch --------------------------------
% Revolute2: FLX002 <-> KT ; Revolute1: FLX002 <-> KB  (from topology dump)
r2 = [mdl '/Revolute2'];
r1 = [mdl '/Revolute1'];
ktPort = line_partner_port(r2, kt);
kbPort = line_partner_port(r1, kb);
fprintf('branch ports resolved\n');

% ---------- build per muscle --------------------------------------------------
for m = 1:numel(muscles)
    mu = muscles(m);
    pfx = mu.name;
    dT = mu.ptKT - fKT;
    dI = mu.ptKB - fKB;
    Rest = 0.1301; Kmax = 0.1073;

    add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/' pfx '_origT'], 'Position', [60 60*m 140 120*m]);
    set_param([mdl '/' pfx '_origT'], 'TranslationMethod', 'Cartesian', ...
        'TranslationCartesianOffset', mat2str(dT, 6), 'RotationMethod', 'None');
    add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/' pfx '_insT'], 'Position', [60 40+60*m 140 100+60*m]);
    set_param([mdl '/' pfx '_insT'], 'TranslationMethod', 'Cartesian', ...
        'TranslationCartesianOffset', mat2str(dI, 6), 'RotationMethod', 'None');

    % transform sensor: base = origin frame, follower = insertion frame
    ts = 'sm_lib/Frames and Transforms/Transform Sensor';
    add_block(ts, [mdl '/' pfx '_TS'], 'Position', [220 40*m 300 120*m]);
    enable_sensor_xyz([mdl '/' pfx '_TS']);

    % BPA force law (Simulink)
    add_block(bpaBlk, [mdl '/' pfx '_BPA'], 'Position', [560 40*m 640 120*m], ...
        'Rest', num2str(Rest), 'Kmax', num2str(Kmax));

    % internal force element
    add_block('sm_lib/Forces and Torques/Internal Force', [mdl '/' pfx '_IF'], 'Position', [900 40*m 980 120*m]);
    print_block_params([mdl '/' pfx '_IF'], pfx);

    % wire frames
    ot = get_param([mdl '/' pfx '_origT'], 'PortHandles'); otP = [ot.RConn ot.LConn];
    it = get_param([mdl '/' pfx '_insT'], 'PortHandles');  itP = [it.RConn it.LConn];
    tsp = get_param([mdl '/' pfx '_TS'], 'PortHandles');   tsP = [tsp.RConn tsp.LConn];
    ifp = get_param([mdl '/' pfx '_IF'], 'PortHandles');   ifP = [ifp.RConn ifp.LConn];
    add_line(mdl, ktPort, otP(1));        % branch thigh frame
    add_line(mdl, otP(2), tsP(1));        % sensor base
    add_line(mdl, otP(2), ifP(1));        % force element base  (branch)
    add_line(mdl, kbPort, itP(1));        % branch shank frame
    add_line(mdl, itP(2), tsP(2));        % sensor follower
    add_line(mdl, itP(2), ifP(2));        % force element follower
    fprintf('%s frames wired\n', pfx);
end

save_system(mdl);
fprintf('=== stage 1 (frames+sensors+forces) saved ===\n');
close_system(mdl, 0);
end

function port = line_partner_port(jointBlk, bodyBlk)
% the port handle ON jointBlk whose line reaches bodyBlk
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

function enable_sensor_xyz(tsBlk)
% enable x/y/z + distance measurement outputs on a Transform Sensor
set_param(tsBlk, 'SenseX', 'on', 'SenseY', 'on', 'SenseZ', 'on', 'SenseDist', 'on');
ph = get_param(tsBlk, 'PortHandles');
allp = [ph.RConn ph.LConn];
nm = strings(0, 1);
for k = 1:numel(allp)
    s = '';
    try, s = get_param(allp(k), 'Name'); catch, end
    nm(k, 1) = string(s);
end
fprintf('TS %s ports: %s\n', strrep(tsBlk, [bdroot '/'], ''), strjoin(nm, ' | '));
end

function print_block_params(blk, tag)
dp = fieldnames(get_param(blk, 'DialogParameters'));
fprintf('%s params: %s\n', tag, strjoin(dp', ' | '));
end
