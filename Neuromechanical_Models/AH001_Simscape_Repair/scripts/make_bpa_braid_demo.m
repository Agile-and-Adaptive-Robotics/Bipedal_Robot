% make_bpa_braid_demo.m  (AH001_Simscape_Repair, 2026-11-01, Ben's direction)
% Frame-mountable, diameter-expanding BPA model + pendulum demo.
%   braid kinematics (inextensible braid, N turns): b^2 = L^2 + (pi*N*D)^2
%     D(L) = sqrt(b^2 - L^2)/(pi*N)      % diameter GROWS as it contracts
%   static force (virtual work):  F = P/(4*pi*N^2) * (3L^2 - b^2)
%     tension only, active while Lmin <= L <= L0 (contraction <= Kmax)
% The BPA assembly is two frame ports (B/F) + Internal Force + Transform
% Sensor + braid math -- the same pattern used for the lab's SNS BPA_20mm,
% extended with the diameter state. Mount = connect B/F to any two robot
% frames.

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'bpa_braid_demo_20261101.txt');
fid = fopen(logFile, 'w');

%% demo parameters
P    = 300e3;      % Pa
D0   = 0.020;      % m braid diameter at rest
L0   = 0.200;      % m rest length (active length ceiling in demo)
Kmax = 0.25;
N    = 1.5;        % braid turns (ASSUMPTION -- replace with measured braid)
piN  = pi*N;
b2   = L0^2 + (piN*D0)^2;     % inextensible braid length squared
kF   = P/(4*piN^2);
Lmin = L0*(1-Kmax);
assignin('base', 'b2', b2); assignin('base', 'kF', kF);
assignin('base', 'piN', piN);
assignin('base', 'L0v', L0); assignin('base', 'Lminv', Lmin);

mdl = 'bpa_braid_demo';
new_system(mdl);

%% world + configs
add_block('sm_lib/Frames and Transforms/World Frame', [mdl '/World']);
add_block('nesl_utility/Solver Configuration',        [mdl '/SolverCfg']);
add_block('sm_lib/Utilities/Mechanism Configuration', [mdl '/MechCfg']);
set_param([mdl '/MechCfg'], 'GravityVector', '[0 -9.80665 0]');
phW = get_param([mdl '/World'], 'PortHandles'); wcon = [phW.LConn phW.RConn];
phS = get_param([mdl '/SolverCfg'], 'PortHandles'); scon = [phS.LConn phS.RConn];
phM = get_param([mdl '/MechCfg'], 'PortHandles'); mcon = [phM.LConn phM.RConn];
add_line(mdl, wcon(1), scon(1));
add_line(mdl, wcon(1), mcon(1));   % branch

%% pendulum: pivot revolute + arm brick
add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/T_piv']);
add_block('sm_lib/Joints/Revolute Joint', [mdl '/pivot']);
set_param([mdl '/pivot'], 'SensePosition', 'on');
set_param([mdl '/pivot'], 'VelocityTargetSpecify', 'on');
set_param([mdl '/pivot'], 'VelocityTargetValue', '-45');
add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/T_armMid'], ...
    'TranslationMethod', 'Cartesian', 'TranslationCartesianOffset', '[125 0 0]', ...
    'TranslationCartesianOffsetUnits', 'mm', 'RotationMethod', 'None');
add_block('sm_lib/Body Elements/Brick Solid', [mdl '/arm']);
set_param([mdl '/arm'], 'BrickDimensions', '[0.25 0.02 0.02]');
set_param([mdl '/arm'], 'BrickDimensionsUnits', 'm');
set_param([mdl '/arm'], 'Density', '7800');
add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/T_Fmount'], ...
    'TranslationMethod', 'Cartesian', 'TranslationCartesianOffset', '[250 0 0]', ...
    'TranslationCartesianOffsetUnits', 'mm', 'RotationMethod', 'None');
add_line(mdl, wcon(1), get_param([mdl '/T_piv'], 'PortHandles').LConn(1));
add_line(mdl, get_param([mdl '/T_piv'],'PortHandles').RConn(1), get_param([mdl '/pivot'],'PortHandles').LConn(1));
add_line(mdl, get_param([mdl '/pivot'],'PortHandles').RConn(1), get_param([mdl '/T_armMid'],'PortHandles').LConn(1));
add_line(mdl, get_param([mdl '/T_armMid'],'PortHandles').RConn(1), solidPort([mdl '/arm']));
add_line(mdl, get_param([mdl '/pivot'],'PortHandles').RConn(1), get_param([mdl '/T_Fmount'],'PortHandles').LConn(1));

%% BPA: B mount above pivot, F mount at arm end
add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/T_Bmount'], ...
    'TranslationMethod', 'Cartesian', 'TranslationCartesianOffset', '[0 50 0]', ...
    'TranslationCartesianOffsetUnits', 'mm', 'RotationMethod', 'None');
add_line(mdl, wcon(1), get_param([mdl '/T_Bmount'], 'PortHandles').LConn(1));
Bport = get_param([mdl '/T_Bmount'], 'PortHandles').RConn(1);
Fport = get_param([mdl '/T_Fmount'], 'PortHandles').RConn(1);

% Internal Force along B-F
add_block('sm_lib/Forces and Torques/Internal Force', [mdl '/BPA_force']);
phIF = get_param([mdl '/BPA_force'], 'PortHandles');
add_line(mdl, Bport, phIF.LConn(1));
add_line(mdl, Fport, phIF.RConn(1));

% Transform Sensor B->F (x,y,z PS outputs)
add_block('sm_lib/Frames and Transforms/Transform Sensor', [mdl '/BPA_sensor']);
phTS = get_param([mdl '/BPA_sensor'], 'PortHandles');
add_line(mdl, Bport, phTS.LConn(1));
add_line(mdl, Fport, phTS.RConn(1));
% TS signal outputs are PS outports (x,y,z + maybe rotation) - find via handles:
outPS = phTS.Outport;
for k = 1:numel(outPS)
    try
        nm = get_param(outPS(k), 'PortType');
        p(fid, 'TS outport %d: %s\n', k, nm);
    catch
    end
end

% PS-Simulink converters for x, y, z
cvt = cell(1,3);
for s = 1:3
    cvt{s} = sprintf('%s/cvt_xyz%d', mdl, s); %#ok<AGROW>
    add_block('nesl_utility/PS-Simulink Converter', cvt{s});
end
add_line(mdl, outPS(1), get_param(cvt{1}, 'PortHandles').Inport);
add_line(mdl, outPS(2), get_param(cvt{2}, 'PortHandles').Inport);
add_line(mdl, outPS(3), get_param(cvt{3}, 'PortHandles').Inport);

% L^2 = x^2+y^2+z^2   (Product with Inputs='**' squares its single input)
add_block('simulink/Math Operations/Product', [mdl '/sq_x'], 'Inputs', '**');
add_block('simulink/Math Operations/Product', [mdl '/sq_y'], 'Inputs', '**');
add_block('simulink/Math Operations/Product', [mdl '/sq_z'], 'Inputs', '**');
add_block('simulink/Math Operations/Sum', [mdl '/sum3'], 'Inputs', '3');
add_line(mdl, out_path(cvt{1}), in_path([mdl '/sq_x'], 1));
add_line(mdl, out_path(cvt{2}), in_path([mdl '/sq_y'], 1));
add_line(mdl, out_path(cvt{3}), in_path([mdl '/sq_z'], 1));
add_line(mdl, out_path([mdl '/sq_x']), in_path([mdl '/sum3'], 1));
add_line(mdl, out_path([mdl '/sq_y']), in_path([mdl '/sum3'], 2));
add_line(mdl, out_path([mdl '/sq_z']), in_path([mdl '/sum3'], 3));

% L = sqrt(L2)
add_block('simulink/Math Operations/Sqrt', [mdl '/L']);
add_line(mdl, out_path([mdl '/sum3']), in_path([mdl '/L'], 1));

% F_raw = kF*(3*L2 - b2)
add_block('simulink/Math Operations/Gain', [mdl '/g3'], 'Gain', '3');
add_block('simulink/Math Operations/Sum', [mdl '/sumF'], 'Inputs', '+-');
add_block('simulink/Math Operations/Gain', [mdl '/gkF'], 'Gain', 'kF');
add_block('simulink/Sources/Constant', [mdl '/b2c'], 'Value', 'b2');
add_line(mdl, out_path([mdl '/sum3']), in_path([mdl '/g3'], 1));
add_line(mdl, out_path([mdl '/g3']), in_path([mdl '/sumF'], 1));
add_line(mdl, out_path([mdl '/b2c']), in_path([mdl '/sumF'], 2));
add_line(mdl, out_path([mdl '/sumF']), in_path([mdl '/gkF'], 1));

% F_pos = max(F_raw, 0)
add_block('simulink/Signal Routing/Switch', [mdl '/posSw'], 'Criteria', 'u2 > Threshold', 'Threshold', '0');
add_block('simulink/Sources/Constant', [mdl '/zero'], 'Value', '0');
add_line(mdl, out_path([mdl '/gkF']), in_path([mdl '/posSw'], 1));
add_line(mdl, out_path([mdl '/gkF']), in_path([mdl '/posSw'], 2));
add_line(mdl, out_path([mdl '/zero']), in_path([mdl '/posSw'], 3));

% gates: L < L0 and L > Lmin
add_block('simulink/Logic and Bit Operations/Compare To Constant', [mdl '/cmpLo'], ...
    'relop', '<', 'Value', 'L0v');
add_block('simulink/Logic and Bit Operations/Compare To Constant', [mdl '/cmpHi'], ...
    'relop', '>', 'Value', 'Lminv');
add_block('simulink/Math Operations/Product', [mdl '/gate'], 'Inputs', '3');
add_line(mdl, out_path([mdl '/L']), in_path([mdl '/cmpLo'], 1));
add_line(mdl, out_path([mdl '/L']), in_path([mdl '/cmpHi'], 1));
add_line(mdl, out_path([mdl '/cmpLo']), in_path([mdl '/gate'], 1));
add_line(mdl, out_path([mdl '/cmpHi']), in_path([mdl '/gate'], 2));
add_line(mdl, out_path([mdl '/posSw']), in_path([mdl '/gate'], 3));

% F -> PS -> Internal Force input (PS ports live in LConn/RConn, not Inport)
add_block('nesl_utility/Simulink-PS Converter', [mdl '/F2ps']);
add_line(mdl, out_path([mdl '/gate']), get_param([mdl '/F2ps'], 'PortHandles').Inport(1));
phF2 = get_param([mdl '/F2ps'], 'PortHandles');
% find the Internal Force PS input: a conserving port that is NOT B or F
phIFx = get_param([mdl '/BPA_force'], 'PortHandles');
psIn = -1;
usedL = get_param(phIFx.LConn(1), 'Line');
usedR = get_param(phIFx.RConn(1), 'Line');
for cnd = [phIFx.LConn phIFx.RConn]
    if cnd ~= -1 && get_param(cnd, 'Line') == -1
        psIn = cnd;
    end
end
if psIn == -1
    % B/F wired first; try indexing: LConn(1)=B, RConn(1)=F, extra RConn(2)=PS
    psIn = phIFx.RConn(min(2, numel(phIFx.RConn)));
end
psOut = [phF2.LConn phF2.RConn];
add_line(mdl, psOut(1), psIn);

% D = sqrt(max(b2 - L2,0))/piN
add_block('simulink/Math Operations/Sum', [mdl '/sumD'], 'Inputs', '+-');
add_block('simulink/Signal Routing/Switch', [mdl '/dSw'], 'Criteria', 'u2 > Threshold', 'Threshold', '0');
add_block('simulink/Math Operations/Sqrt', [mdl '/sqrtD']);
add_block('simulink/Math Operations/Gain', [mdl '/gD'], 'Gain', '1/piN');
add_line(mdl, out_path([mdl '/b2c']), in_path([mdl '/sumD'], 1));
add_line(mdl, out_path([mdl '/sum3']), in_path([mdl '/sumD'], 2));
add_line(mdl, out_path([mdl '/sumD']), in_path([mdl '/dSw'], 1));
add_line(mdl, out_path([mdl '/sumD']), in_path([mdl '/dSw'], 2));
add_line(mdl, out_path([mdl '/zero']), in_path([mdl '/dSw'], 3));
add_line(mdl, out_path([mdl '/dSw']), in_path([mdl '/sqrtD'], 1));
add_line(mdl, out_path([mdl '/sqrtD']), in_path([mdl '/gD'], 1));

% logging
logw(mdl, [mdl '/L'], 'log_L');
logw(mdl, [mdl '/gD'], 'log_D');
logw(mdl, [mdl '/gate'], 'log_F');
% joint angle: pivot PS position out -> converter -> workspace
qPS = phIF; %#ok<NASGU>
phP = get_param([mdl '/pivot'], 'PortHandles');
add_block('nesl_utility/PS-Simulink Converter', [mdl '/cvt_q']);
add_line(mdl, phP.Outport(1), get_param([mdl '/cvt_q'], 'PortHandles').Inport);
logw(mdl, [mdl '/cvt_q'], 'log_th');

set_param(mdl, 'StopTime', '2.0', 'Solver', 'ode23t');
save_system(mdl, fullfile(repairRoot, 'working', [mdl '.slx']));

%% run + report
out = sim(mdl, 'StopTime', '2.0');
L = evalin('base', 'log_L'); D = evalin('base', 'log_D');
F = evalin('base', 'log_F'); th = evalin('base', 'log_th');
Lv = squeeze(L.signals.values); Dv = squeeze(D.signals.values);
Fv = squeeze(F.signals.values); tv = squeeze(th.signals.values);
p(fid, 'L range: %.4f .. %.4f m\n', min(Lv), max(Lv));
p(fid, 'D range: %.5f .. %.5f m\n', min(Dv), max(Dv));
p(fid, 'F range: %.2f .. %.2f N\n', min(Fv), max(Fv));
% correlation check: where L shrinks, D must grow
[~, imin] = min(Lv); [~, imax] = max(Lv);
p(fid, 'at min L=%.4f -> D=%.5f ; at max L=%.4f -> D=%.5f (must be inverse)\n', ...
    Lv(imin), Dv(imin), Lv(imax), Dv(imax));
p(fid, 'pivot sweep: %.1f -> %.1f deg\n', tv(1), tv(end));
close_system(mdl, 0);
p(fid, 'BPA BRAID DEMO DONE\n');
fclose(fid);

function h = solidPort(blk)
ph = get_param(blk, 'PortHandles');
c = [ph.LConn ph.RConn];
h = c(1);
end

function wire(mdl, srcBlk, outIdx, dstBlk)
add_line(mdl, out_path(srcBlk), in_path(dstBlk, outIdx));
end

function o = out_path(b)
o = [b '/1'];
end

function i = in_path(b, k)
i = sprintf('%s/%d', b, k);
end

function logw(mdl, blk, varname)
add_block('simulink/Sinks/To Workspace', [mdl '/' varname], 'VariableName', varname, 'SaveFormat', 'Structure With Time');
add_line(mdl, out_path(blk), in_path([mdl '/' varname], 1));
end

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
