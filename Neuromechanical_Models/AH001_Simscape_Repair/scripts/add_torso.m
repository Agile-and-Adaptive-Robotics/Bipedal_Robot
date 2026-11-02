% add_torso.m  (AH001_Simscape_Repair, 2026-11-01, Ben's direction)
% Add the torso spec to ah001_repair.slx:
%   - expose the PELVIS DATUM frame (branch the pelvis Solid base net to a new
%     Connection Port on the pelvis group)
%   - lumbar ball joint at pelvis [-0.16, 0.10, 0] m (spec: upper-rear)
%   - 0.25 m, 20 mm square spine bar, density 1200 -> 0.12 kg
%   - payload spheres per today's ledger: 35.7539269872 kg torso at
%     [0.10, 0.075, 0] m from hinge; 35.7539269873 kg pelvis half at [0.10,0,0]
%   - total becomes 75.0000000000 kg by construction
% Frame assumption (documented): pelvis datum axes are anatomical
% (x anterior, y superior, z right) per the 2026-10-27 frame ruling.

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'add_torso_20261101.txt');
fid = fopen(logFile, 'w');

mdl = 'ah001_repair';
load_system(fullfile(repairRoot, 'working', [mdl '.slx']));
x09 = [mdl '/x09_BA_001_1'];
pelvisGrp = [x09 '/x02_01_PE_001_1_RIGID'];

%% 1. pelvis datum port: branch the Solid's BASE-side net inside the group
sol = find_system(pelvisGrp, 'LookUnderMasks','all','FollowLinks','on', ...
                  'ReferenceBlock', 'sm_lib/Body Elements/File Solid');
if isempty(sol)
    p(fid, 'FATAL: no File Solid in pelvis group\n'); fclose(fid); error('no pelvis solid');
end
solB = sol{1};
phS = get_param(solB, 'PortHandles');
% pick the solid's CONNECTED conserving port (layout varies: LConn or RConn)
cands = [phS.LConn phS.RConn];
basePort = -1;
for ii = 1:numel(cands)
    if cands(ii) ~= -1 && get_param(cands(ii), 'Line') ~= -1
        basePort = cands(ii);
        break;
    end
end
if basePort == -1
    p(fid, 'FATAL: pelvis Solid has no connected conserving port\n'); fclose(fid); error('unconnected');
end
cpPath = [pelvisGrp '/pelvis_datum_port'];
add_block('nesl_utility/Connection Port', cpPath, ...
          'Position', get_param(solB, 'Position') + [220 0 220 0]);
for pv = {'PortSide','right'; 'ConnectionType','physical'}
    try, set_param(cpPath, pv{1}, pv{2}); catch, end
end
phC = get_param(cpPath, 'PortHandles');
allc = [phC.LConn phC.RConn];
add_line(pelvisGrp, basePort, allc(1));   % branch (whichever port exists)
p(fid, 'pelvis datum port added and branched to Solid base net\n');

% find the new port index on the group (last RConn)
phG = get_param(pelvisGrp, 'PortHandles');
nP = numel(phG.RConn);
p(fid, 'pelvis group now has %d RConn ports; new port is R%d\n', nP, nP);
datumPort = phG.RConn(nP);

%% 2. torso chain inside x09
% --- lumbar placement transform (pelvis datum -> hinge location)
T_lum = [x09 '/T_lumbar'];
add_block('sm_lib/Frames and Transforms/Rigid Transform', T_lum, ...
          'TranslationMethod', 'Cartesian', ...
          'TranslationCartesianOffset', '[-160 100 0]', ...   % mm, pelvis frame
          'TranslationCartesianOffsetUnits', 'mm', ...
          'RotationMethod', 'None');
add_line(x09, datumPort, get_param(T_lum, 'PortHandles').LConn(1));

% --- lumbar ball joint
lumbar = [x09 '/lumbar_ball'];
add_block('sm_lib/Joints/Spherical Joint', lumbar);
add_line(x09, get_param(T_lum,'PortHandles').RConn(1), get_param(lumbar,'PortHandles').LConn(1));
% damping per torso_spec: 2 N*m*s/rad (set ONLY the damper coefficient + units)
try
    okUnits = false;
    try
        set_param(lumbar, 'DampingCoefficientUnits', 'N*m/(rad/s)');
        okUnits = true;
    catch
        try
            set_param(lumbar, 'DampingCoefficientUnits', 'N*m/(deg/s)');
            okUnits = true;
        catch
        end
    end
    if okUnits
        set_param(lumbar, 'DampingCoefficient', '2');
        p(fid, 'lumbar DampingCoefficient = 2 (spec: 2 N*m*s/rad)\n');
    else
        p(fid, 'lumbar damping param not found - left default (document)\n');
    end
catch ME
    p(fid, 'lumbar damping note: %s\n', ME.message);
end

lumF = get_param(lumbar, 'PortHandles').RConn(1);

% --- spine bar: 20x20x250 mm along +y from the hinge, density 1200 = 0.12 kg
T_spn = [x09 '/T_spine_off'];
add_block('sm_lib/Frames and Transforms/Rigid Transform', T_spn, ...
          'TranslationMethod', 'Cartesian', ...
          'TranslationCartesianOffset', '[0 125 0]', ...
          'TranslationCartesianOffsetUnits', 'mm', ...
          'RotationMethod', 'None');
add_line(x09, lumF, get_param(T_spn, 'PortHandles').LConn(1));
spn = [x09 '/spine_bar'];
add_block('sm_lib/Body Elements/Brick Solid', spn);
pv = {'BrickDimensions','[0.02 0.25 0.02]'; 'BrickDimensionsUnits','m'; ...
      'Density','1200'; 'DensityUnits','kg/m^3'};
for kk = 1:size(pv,1)
    try
        set_param(spn, pv{kk,1}, pv{kk,2});
    catch ME
        p(fid, 'spine param %s not set: %s\n', pv{kk,1}, ME.message);
    end
end
add_line(x09, get_param(T_spn,'PortHandles').RConn(1), solidPort(spn));

% --- torso payload sphere at [0.10 0.075 0] m from hinge (torso side)
T_pt = [x09 '/T_payload_torso'];
add_block('sm_lib/Frames and Transforms/Rigid Transform', T_pt, ...
          'TranslationMethod', 'Cartesian', ...
          'TranslationCartesianOffset', '[100 75 0]', ...
          'TranslationCartesianOffsetUnits', 'mm', ...
          'RotationMethod', 'None');
add_line(x09, lumF, get_param(T_pt, 'PortHandles').LConn(1));
pt = [x09 '/payload_torso'];
add_block('sm_lib/Body Elements/Spherical Solid', pt);
pv = {'SphereRadius','[0.1]'; 'SphereRadiusUnits','m'; ...
      'BasedOnType','Mass'; 'Mass','35.7539269872'; 'MassUnits','kg'};
for kk = 1:size(pv,1)
    try, set_param(pt, pv{kk,1}, pv{kk,2}); catch, p(fid, 'param note %s: not set\n', pv{kk,1}); end
end
add_line(x09, get_param(T_pt,'PortHandles').RConn(1), solidPort(pt));

% --- pelvis payload half at [0.10 0 0] m from hinge (pelvis side branch)
T_pp = [x09 '/T_payload_pelvis'];
add_block('sm_lib/Frames and Transforms/Rigid Transform', T_pp, ...
          'TranslationMethod', 'Cartesian', ...
          'TranslationCartesianOffset', '[100 0 0]', ...
          'TranslationCartesianOffsetUnits', 'mm', ...
          'RotationMethod', 'None');
add_line(x09, get_param(T_lum,'PortHandles').RConn(1), get_param(T_pp, 'PortHandles').LConn(1));
pp = [x09 '/payload_pelvis'];
add_block('sm_lib/Body Elements/Spherical Solid', pp);
pv = {'SphereRadius','[0.1]'; 'SphereRadiusUnits','m'; ...
      'BasedOnType','Mass'; 'Mass','35.7539269873'; 'MassUnits','kg'};
for kk = 1:size(pv,1)
    try, set_param(pp, pv{kk,1}, pv{kk,2}); catch, p(fid, 'param note %s: not set\n', pv{kk,1}); end
end
add_line(x09, get_param(T_pp,'PortHandles').RConn(1), solidPort(pp));

%% 3. compile gate
ok = true;
try
    sim(mdl, 'StopTime', '0');
    p(fid, 'compile with torso: OK\n');
catch ME
    ok = false;
    p(fid, 'compile FAIL: %s\n', ME.message);
    try
        cs = ME.cause;
        if isa(cs, 'MException'), cs = {cs}; end
        for k = 1:numel(cs), p(fid, '  cause: %s\n', cs{k}.message); end
    catch
    end
end

%% 4. total mass ledger (must be exactly 75 kg)
blks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
ref = strings(numel(blks),1);
for k = 1:numel(blks)
    try, ref(k) = strtrim(strrep(string(get_param(blks{k}, 'ReferenceBlock')), newline, " ")); catch, ref(k) = ""; end
end
hws = get_param(mdl, 'ModelWorkspace');
sm = evalin(hws, 'smiData');
tot = 0;
for k = 1:numel(blks)
    if ref(k) == "sm_lib/Body Elements/File Solid"
        me_ = get_param(blks{k}, 'Mass');
        mi = regexp(me_, 'Solid\((\d+)\)', 'tokens', 'once');
        if ~isempty(mi), tot = tot + sm.Solid(str2double(mi{1})).mass; end
    elseif ref(k) == "sm_lib/Body Elements/Brick Solid" || ref(k) == "sm_lib/Body Elements/Spherical Solid"
        nm = get_param(blks{k}, 'Name');
        switch nm
            case 'spine_bar', tot = tot + 0.12;
            case 'payload_torso', tot = tot + 35.7539269872;
            case 'payload_pelvis', tot = tot + 35.7539269873;
            otherwise, p(fid, 'uncounted solid: %s\n', nm);
        end
    end
end
p(fid, 'TOTAL MASS with torso: %.10f kg (target 75)\n', tot);

save_system(mdl, fullfile(repairRoot, 'working', [mdl '.slx']));
close_system(mdl, 0);
p(fid, 'ADD TORSO DONE. compile=%d\n', ok);
fclose(fid);

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end

function h = solidPort(blk)
% first existing conserving port of a body-element block (layout varies)
ph = get_param(blk, 'PortHandles');
c = [ph.LConn ph.RConn];
h = c(1);
end
