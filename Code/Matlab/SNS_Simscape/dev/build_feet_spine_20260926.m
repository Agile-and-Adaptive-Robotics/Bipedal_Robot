function build_feet_spine_20260926()
%% build_feet_spine_20260926  Add simple block feet (welded at the tibia
%% bottom) and a spine (cylinder column + lumped mass) to the humanoid so the
%% total mass is 73 kg. Feet: 0.24x0.035x0.09 m Onyx bricks; spine: r=0.06,
%% h=0.45 m Onyx cylinder + 61.09 kg point mass (trunk+head+arms surrogate).

sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));

ONYX = 1200;   % kg/m^3 (Markforged datasheet; solid-print assumption)

% ---- branch frames ----------------------------------------------------------
kbR = [sub '/x04_02_KB_R_001_1_RIGID'];
kbL = [sub '/x04_04_KB_L_001_1_RIGID'];
pe  = [sub '/x02_01_PE_001_1_RIGID'];
phR = get_param(kbR, 'PortHandles'); kbRPort = [phR.RConn phR.LConn]; kbRPort = kbRPort(1);  % -> TI bolt line
phL = get_param(kbL, 'PortHandles'); kbLPort = [phL.RConn phL.LConn]; kbLPort = kbLPort(1);
phP = get_param(pe,  'PortHandles'); pePort  = [phP.RConn phP.LConn]; pePort  = pePort(2);   % -> Spherical (KNOB1)

% ---- feet -------------------------------------------------------------------
add_block('sm_lib/Frames and Transforms/Rigid Transform', [sub '/foot_T_R'], 'Position', [660 40 740 100]);
set_param([sub '/foot_T_R'], 'TranslationMethod', 'Cartesian', ...
    'TranslationCartesianOffset', '[0.03992 -0.3323 0.00391]', ...
    'RotationMethod', 'None');
add_block('sm_lib/Body Elements/Brick Solid', [sub '/foot_R'], 'Position', [780 40 860 100]);

add_block('sm_lib/Frames and Transforms/Rigid Transform', [sub '/foot_T_L'], 'Position', [660 140 740 200]);
set_param([sub '/foot_T_L'], 'TranslationMethod', 'Cartesian', ...
    'TranslationCartesianOffset', '[0.03992 0.3323 0.00391]', ...
    'RotationMethod', 'ArbitraryAxis', 'RotationAngleUnits', 'rad', ...
    'RotationAngle', '3.14159265358979', 'RotationArbitraryAxis', '[1 0 0]');
add_block('sm_lib/Body Elements/Brick Solid', [sub '/foot_L'], 'Position', [780 140 860 200]);

for s = {'R', 'L'}
    blk = [sub '/foot_' s{1}];
    dp = fieldnames(get_param(blk, 'DialogParameters'));
    fprintf('foot_%s dialog params: %s\n', s{1}, strjoin(dp', ', '));
end

% ---- spine ------------------------------------------------------------------
add_block('sm_lib/Frames and Transforms/Rigid Transform', [sub '/spine_T'], 'Position', [660 260 740 320]);
set_param([sub '/spine_T'], 'TranslationMethod', 'Cartesian', ...
    'TranslationCartesianOffset', '[-0.017956 0.374509 -0.065]', ...
    'RotationMethod', 'ArbitraryAxis', 'RotationAngleUnits', 'rad', ...
    'RotationAngle', '-1.5708', 'RotationArbitraryAxis', '[1 0 0]');
add_block('sm_lib/Body Elements/Cylindrical Solid', [sub '/spine_col'], 'Position', [780 260 860 320]);

add_block('sm_lib/Frames and Transforms/Rigid Transform', [sub '/lump_T'], 'Position', [660 380 740 440]);
set_param([sub '/lump_T'], 'TranslationMethod', 'Cartesian', ...
    'TranslationCartesianOffset', '[-0.017956 0.399509 -0.065]', ...
    'RotationMethod', 'None');
add_block('sm_lib/Body Elements/Brick Solid', [sub '/spine_lump'], 'Position', [780 380 860 440]);

% ---- set solid parameters (schema discovered at runtime) --------------------
% Brick: size + density
try, set_param([sub '/foot_R'], 'Size', '[0.24 0.035 0.09]'); catch ME, fprintf('brick size: %s\n', ME.message); end
try, set_param([sub '/foot_R'], 'Density', num2str(ONYX)); catch ME, fprintf('brick density: %s\n', ME.message); end
try, set_param([sub '/foot_L'], 'Size', '[0.24 0.035 0.09]'); catch, end
try, set_param([sub '/foot_L'], 'Density', num2str(ONYX)); catch, end
% Cylinder: radius/length + density
try, set_param([sub '/spine_col'], 'Radius', '0.06'); catch ME, fprintf('cyl radius: %s\n', ME.message); end
try, set_param([sub '/spine_col'], 'Length', '0.45'); catch, end
try, set_param([sub '/spine_col'], 'Density', num2str(ONYX)); catch, end
% Point Mass: mass
try, set_param([sub '/spine_lump'], 'Mass', '61.0906438814214'); catch ME, fprintf('pm mass: %s\n', ME.message); end

% ---- wire: branch off the existing bolt/spherical lines ---------------------
ftR = get_param([sub '/foot_T_R'], 'PortHandles'); ftRP = [ftR.RConn ftR.LConn];
ftL = get_param([sub '/foot_T_L'], 'PortHandles'); ftLP = [ftL.RConn ftL.LConn];
spT = get_param([sub '/spine_T'],  'PortHandles'); spTP = [spT.RConn spT.LConn];
lpT = get_param([sub '/lump_T'],   'PortHandles'); lpTP = [lpT.RConn lpT.LConn];
fR  = get_param([sub '/foot_R'],   'PortHandles'); fRP  = [fR.RConn fR.LConn];
fL  = get_param([sub '/foot_L'],   'PortHandles'); fLP  = [fL.RConn fL.LConn];
sc  = get_param([sub '/spine_col'],'PortHandles'); scP  = [sc.RConn sc.LConn];
sl  = get_param([sub '/spine_lump'],'PortHandles'); slP = [sl.RConn sl.LConn];

add_line(sub, kbRPort, ftRP(1));   % branch: KB<->TI bolt  -> foot_R transform
add_line(sub, ftRP(2), fRP(1));
add_line(sub, kbLPort, ftLP(1));
add_line(sub, ftLP(2), fLP(1));
add_line(sub, pePort,  spTP(1));   % branch: PE<->KNOB spherical frame -> spine
add_line(sub, spTP(2), scP(1));
add_line(sub, pePort,  lpTP(1));
add_line(sub, lpTP(2), slP(1));
fprintf('feet + spine wired\n');

% ---- compile + save ----------------------------------------------------------
try
    set_param(mdl, 'StopTime', '0.01');
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('%s: UPDATE OK\n', mdl);
catch ME
    fprintf('%s: UPDATE FAILED: %s\n', mdl, ME.message);
    for c = 1:min(numel(ME.cause), 8)
        fprintf('  CAUSE %d: %s\n', c, ME.cause{c}.message(1:min(end, 180)));
    end
end
save_system(mdl);
close_system(mdl, 0);
fprintf('=== build_feet_spine DONE ===\n');
end
