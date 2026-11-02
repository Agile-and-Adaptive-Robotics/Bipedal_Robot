% add_feet.m  (AH001_Simscape_Repair, 2026-11-02 -- Ben's overnight assignment)
% Attach BOTH articulated feet to ah001_repair.slx:
%   per foot: 2-DOF subtalar-calcaneus (Exterior -rev1- Interior -rev2- BottomJoint)
%             1-DOF MTP (mainplate -rev- axle1 -rev- frontfoot)
%   toes rigid on forefoot, mid foot rigid on rearfoot. NO welded feet.
% Frames = the CAD's own mate transforms (T26/T25, T11/T12, T39/T40, T21/T22);
% rigid placements = InstanceTree poses relative to each segment's anchor part.
% Mount = Ankle Exterior frame on the tibia part datum, IDENTITY provisional
% transform (single block to adjust once the true mount datum is measured).

function add_feet()
global poses mateT STEP massOf repairRoot fid plog

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'add_feet_20261102.txt');
fid = fopen(logFile, 'w');
plog = @pfun;

mdl = 'ah001_repair';
load_system(fullfile(repairRoot, 'working', [mdl '.slx']));
x09 = [mdl '/x09_BA_001_1'];
stepDir = fullfile(repairRoot, 'cad', 'export_foot');

%% data
D = jsondecode(fileread(fullfile(stepDir, 'foot_build_data.json')));
poses = D.poses;

massOf = containers.Map('KeyType','char','ValueType','double');
mk = { 'Ankle Middle Joint Shell',0.012631857024606435; 'Ankle Exterior',0.12192384008158326;
       'Toe',0.0037914944747527227; 'foot axle',0.0034815476923724066;
       'Ankle Bottom Joint',0.013168929260936995; 'New Front foot',0.033435102009916902;
       'Ankle Interior Joint',0.01794923410936473; 'New foot - Copy - Copy',0.11544451399879704;
       'Left mid foot',0.028945071427893498; 'small axle lock',0.00025608685660253341;
       'small axle',0.002487204744708915; 'large axle lock',0.00044007561614053812;
       'Foot T Lock',0.00015959290630236151; 'Ankle Lock 1',0.0005912645598564407;
       'small axle - Copy',0.0016673998288551384 };
for k = 1:size(mk,1), massOf(mk{k,1}) = mk{k,2}; end

STEP = containers.Map('KeyType','char','ValueType','char');
dd = dir(fullfile(stepDir, '*_Default_sldprt.STEP'));
for k = 1:numel(dd)
    STEP(strrep(dd(k).name, '_Default_sldprt.STEP', '')) = fullfile(stepDir, dd(k).name);
end

mateT = containers.Map('KeyType','char','ValueType','any');
df = fileread(fullfile(stepDir, 'foot_ankle_DataFile.m'));
toks = regexp(df, 'smiData\.RigidTransform\((\d+)\)\.(translation|angle|axis) = (.+?);', 'tokens');
for k = 1:numel(toks)
    key = ['T' toks{k}{1}];
    if ~isKey(mateT, key), mateT(key) = struct('translation',[0 0 0],'angle',0,'axis',[0 0 0]); end
    s = mateT(key);
    s.(toks{k}{2}) = str2num(toks{k}{3}); %#ok<ST2NM>
    mateT(key) = s;
end
pfun('mate transforms: %d\n', mateT.Count);

HOUS = {'Ankle Middle Joint Shell-1','Ankle Middle Joint Shell'; 'Ankle Lock 1-1','Ankle Lock 1'};
REAR = {'New foot - Copy - Copy-1','New foot - Copy - Copy'; 'Left mid foot-1','Left mid foot'; ...
        'Foot T Lock-1','Foot T Lock'; 'foot axle-2','foot axle'; 'small axle-1','small axle'; ...
        'small axle - Copy-2','small axle - Copy'; 'small axle lock-1','small axle lock'; ...
        'small axle lock-2','small axle lock'; 'large axle lock-1','large axle lock'; ...
        'large axle lock-2','large axle lock'};
FORE = {'Toe-1','Toe'; 'Toe-2','Toe'; 'Toe-3','Toe'; 'Toe-4','Toe'; 'Toe-5','Toe'};

% ---- anatomical basis of the leg CAD (verified: det+1, FL anterior >0) ----
% superior = mean tibia axis (knee<-distal); lateral = hip-ball line (knobs);
% anterior = superior x lateral. KT/KB are TOP/BOTTOM brackets (Ben) -- the
% anterior sign comes from the four-bar FRONT link FL vs back link BL.
uSup = [0.7784 -0.3333 0.5320];
zlat = [0.6108  0.2063 -0.7644];
uSup = uSup / norm(uSup);
zlat = zlat - dot(zlat, uSup) * uSup; zlat = zlat / norm(zlat);
xant = cross(uSup, zlat);            % right-handed ISB: x = y x z
Q = [xant; uSup; zlat];              % CAD anatomical axes as ISB rows
assert(abs(det(Q) - 1) < 1e-9, 'Q not proper');

% reorient the WHOLE assembly (Ben: align x-forward, y-up once). The rotation
% must sit AFTER the 6-DOF follower: rotating the root Transform would rotate
% robot AND joint frame together, leaving the readings unchanged.
set_param([mdl '/Transform'], 'RotationMethod', 'None');   % root back to identity
six = [mdl '/SixDOF'];
ph6 = get_param(six, 'PortHandles');
lF = get_param(ph6.RConn(1), 'Line');
farF = [];
srcF = get_param(lF, 'SrcPortHandle'); dstF = reshape(get_param(lF, 'DstPortHandle'), 1, []);
for xx = [srcF dstF]
    if xx ~= -1 && xx ~= ph6.RConn(1), farF = xx; end
end
delete_line(lF);
axQ = rotm2axang(Q);
twc = [mdl '/T_world_cad'];
add_block('sm_lib/Frames and Transforms/Rigid Transform', twc);
set_param(twc, 'RotationMethod', 'ArbitraryAxis', ...
    'RotationArbitraryAxis', mat2str(axQ(1:3), 12), ...
    'RotationAngle', mat2str(axQ(4), 12), 'RotationAngleUnits', 'rad', ...
    'TranslationMethod', 'Cartesian', 'TranslationCartesianOffset', '[0 0 0]', ...
    'TranslationCartesianOffsetUnits', 'm');
phT = get_param(twc, 'PortHandles');
add_line(mdl, ph6.RConn(1), phT.LConn(1));
add_line(mdl, phT.RConn(1), farF);
pfun('world frame reoriented via T_world_cad: x fwd, y up, z right\n');

% foot frame (from instance map): x = width, y = UP, z = ANTERIOR.
% foot->ISB: y->+Y, z->+X, x->-Z (right-handed; LEFT-foot design => +x = left)
Rf2i = [0 0 1; 0 1 0; -1 0 0];
assert(abs(det(Rf2i) - 1) < 1e-9, 'Rf2i not proper');
Rleg = Q' * Rf2i;    % foot -> ISB -> CAD (mount lives in the CAD tibia frame)
pfun('mount rotation (foot->CAD tibia frame):\n'); disp(Rleg);

buildFoot('fr', x09, 'x05_01_TI_R_001_1_RIGID', HOUS, REAR, FORE, Rleg);
buildFoot('fl', x09, 'x05_02_TI_001_1_RIGID', HOUS, REAR, FORE, Rleg);

%% payload: keep total exactly 75.0000000000
footTot = massOf('Ankle Exterior')+massOf('Ankle Middle Joint Shell')+massOf('Ankle Lock 1') ...
        + massOf('Ankle Interior Joint') ...
        + massOf('Ankle Bottom Joint')+massOf('New foot - Copy - Copy')+massOf('Left mid foot') ...
        + massOf('Foot T Lock')+2*massOf('foot axle')+massOf('small axle')+massOf('small axle - Copy') ...
        + 2*massOf('small axle lock')+2*massOf('large axle lock') ...
        + massOf('New Front foot')+5*massOf('Toe');
pfun('foot mass each: %.10f kg\n', footTot);
half = (75 - 3.3721460255 - 0.12 - 2*footTot) / 2;
set_param([x09 '/payload_torso'],  'Mass', num2str(half, '%.10f'));
set_param([x09 '/payload_pelvis'], 'Mass', num2str(half, '%.10f'));
pfun('payload spheres -> %.10f each (total exactly 75)\n', half);

%% compile + census
ok = true;
try
    sim(mdl, 'StopTime', '0');
    pfun('compile with feet: OK\n');
catch ME
    ok = false;
    pfun('compile FAIL: %s\n', ME.message);
    try
        cs = ME.cause; if isa(cs,'MException'), cs = {cs}; end
        for k = 1:numel(cs), pfun('  cause: %s\n', cs{k}.message); end
    catch
    end
end
blks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
ref = strings(numel(blks),1);
for k = 1:numel(blks)
    try, ref(k) = strtrim(strrep(string(get_param(blks{k}, 'ReferenceBlock')), newline, " ")); catch, ref(k) = ""; end
end
pfun('joint census:\n');
uj = unique(ref(startsWith(ref, "sm_lib/Joints/")));
for t = uj'
    pfun('  %2d x %s\n', nnz(ref == t), t);
end
save_system(mdl, fullfile(repairRoot, 'working', [mdl '.slx']));
close_system(mdl, 0);
pfun('ADD FEET DONE. compile=%d\n', ok);
fclose(fid);
end

function buildFoot(pre, x09, tibiaGrpName, HOUS, REAR, FORE, Rleg)
global poses mateT STEP massOf fid plog
grp = [x09 '/' tibiaGrpName];
sol = find_system(grp, 'LookUnderMasks','all','FollowLinks','on', ...
                  'ReferenceBlock', 'sm_lib/Body Elements/File Solid');
phS = get_param(sol{1}, 'PortHandles');
cands = [phS.LConn phS.RConn];
basePort = -1;
for ii = 1:numel(cands)
    if cands(ii) ~= -1 && get_param(cands(ii), 'Line') ~= -1, basePort = cands(ii); break; end
end
cp = [grp '/' pre '_datum_port'];
add_block('nesl_utility/Connection Port', cp);
phC = get_param(cp, 'PortHandles'); ccl = [phC.LConn phC.RConn];
add_line(grp, basePort, ccl(1));
phG = get_param(grp, 'PortHandles');
datum = -1;
for pp = [phG.LConn phG.RConn]
    if pp ~= -1 && get_param(pp, 'Line') == -1, datum = pp; end
end

addRT([pre '_mount'], x09, [0 0 0]', Rleg);   % anatomical basis; position PROVISIONAL at tibia datum
add_line(x09, datum, outP(x09, [pre '_mount']));
addFS([pre '_ext'], x09, 'Ankle Exterior');
add_line(x09, outP(x09, [pre '_mount']), solP(x09, [pre '_ext']));

% housing rigid parts (rel Exterior)
rigidSeg(x09, pre, 'h', HOUS, 'Ankle Exterior-1', [pre '_ext']);

% ankle axis 1
addMateRT([pre '_a1F'], x09, 'T26');
addJ([pre '_ank1'], x09);
addMateRT([pre '_a1B'], x09, 'T25');
add_line(x09, solP(x09, [pre '_ext']), inP(x09, [pre '_a1F']));
add_line(x09, outP(x09, [pre '_a1F']), inP(x09, [pre '_ank1']));
add_line(x09, outP(x09, [pre '_ank1']), inP(x09, [pre '_a1B']));
addFS([pre '_spider'], x09, 'Ankle Interior Joint');
add_line(x09, outP(x09, [pre '_a1B']), solP(x09, [pre '_spider']));

% ankle axis 2 (T11 hangs on the spider frame net)
addMateRT([pre '_a2B'], x09, 'T11');
addJ([pre '_ank2'], x09);
addMateRT([pre '_a2F'], x09, 'T12');
add_line(x09, solP(x09, [pre '_spider']), inP(x09, [pre '_a2B']));
add_line(x09, outP(x09, [pre '_a2B']), inP(x09, [pre '_ank2']));
add_line(x09, outP(x09, [pre '_ank2']), inP(x09, [pre '_a2F']));
addFS([pre '_bot'], x09, 'Ankle Bottom Joint');
add_line(x09, outP(x09, [pre '_a2F']), solP(x09, [pre '_bot']));

% rear segment rigid parts (rel Bottom Joint), incl. MAIN PLATE
rigidSeg(x09, pre, 'r', REAR, 'Ankle Bottom Joint-1', [pre '_bot']);

% MTP: T39 must hang on the MAIN PLATE frame (block name = first rear part)
mainBlk = [pre '_r1'];   % rigidSeg adds parts as pre_r1.. with main first
addMateRT([pre '_m1B'], x09, 'T39');
addJ([pre '_mtpA'], x09);
addMateRT([pre '_m1F'], x09, 'T40');
add_line(x09, solP(x09, mainBlk), inP(x09, [pre '_m1B']));
add_line(x09, outP(x09, [pre '_m1B']), inP(x09, [pre '_mtpA']));
add_line(x09, outP(x09, [pre '_mtpA']), inP(x09, [pre '_m1F']));
addFS([pre '_axle1'], x09, 'foot axle');
add_line(x09, outP(x09, [pre '_m1F']), solP(x09, [pre '_axle1']));

% forefoot through the second axle revolute
addMateRT([pre '_m2F'], x09, 'T22');
addJ([pre '_mtpB'], x09);
addMateRT([pre '_m2B'], x09, 'T21');
addFS([pre '_front'], x09, 'New Front foot');
add_line(x09, solP(x09, [pre '_axle1']), inP(x09, [pre '_m2F']));
add_line(x09, outP(x09, [pre '_m2F']), inP(x09, [pre '_mtpB']));
add_line(x09, outP(x09, [pre '_mtpB']), inP(x09, [pre '_m2B']));
add_line(x09, outP(x09, [pre '_m2B']), solP(x09, [pre '_front']));

% toes rigid on forefoot (rel New Front foot)
rigidSeg(x09, pre, 'f', FORE, 'New Front foot-1', [pre '_front']);
pfun('%s built\n', pre);
end

function rigidSeg(x09, pre, tag, list, anchorInst, anchorSolidBlk)
global poses mateT STEP massOf
[Ra, ta] = poseOf(anchorInst);
for k = 1:size(list,1)
    inst = list{k,1}; part = list{k,2};
    [Rp, tp] = poseOf(inst);
    Rrel = Ra' * Rp; trel = Ra' * (tp - ta);
    nm = sprintf('%s_%s%d', pre, tag, k);
    addRT(nm, x09, trel, Rrel);
    addFS([nm '_s'], x09, part);
    add_line(x09, solP(x09, anchorSolidBlk), inP(x09, nm));
    add_line(x09, outP(x09, nm), solP(x09, [nm '_s']));
end
end

function [R, t] = poseOf(name)
global poses
for k = 1:numel(poses)
    if strcmp(poses(k).inst, name)
        R = poses(k).R; t = poses(k).t; return;
    end
end
error('pose not found: %s', name);
end

function addRT(name, sys, t_m, R)
add_block('sm_lib/Frames and Transforms/Rigid Transform', [sys '/' name]);
if norm(R - eye(3)) < 1e-12
    set_param([sys '/' name], 'RotationMethod', 'None');
else
    ax = rotm2axang(R);
    set_param([sys '/' name], 'RotationMethod', 'ArbitraryAxis', ...
        'RotationArbitraryAxis', mat2str(ax(1:3), 12), ...
        'RotationAngle', mat2str(ax(4), 12), 'RotationAngleUnits', 'rad');
end
set_param([sys '/' name], 'TranslationMethod', 'Cartesian', ...
    'TranslationCartesianOffset', mat2str(t_m(:)', 12), ...
    'TranslationCartesianOffsetUnits', 'm');
end

function addMateRT(name, sys, tkey)
global mateT
s = mateT(tkey);
add_block('sm_lib/Frames and Transforms/Rigid Transform', [sys '/' name]);
if abs(s.angle) < 1e-15
    set_param([sys '/' name], 'RotationMethod', 'None');
else
    set_param([sys '/' name], 'RotationMethod', 'ArbitraryAxis', ...
        'RotationArbitraryAxis', mat2str(s.axis, 12), ...
        'RotationAngle', num2str(s.angle, '%.12g'), 'RotationAngleUnits', 'rad');
end
set_param([sys '/' name], 'TranslationMethod', 'Cartesian', ...
    'TranslationCartesianOffset', mat2str(s.translation/1000, 12), ...
    'TranslationCartesianOffsetUnits', 'm');
end

function addJ(name, sys)
add_block('sm_lib/Joints/Revolute Joint', [sys '/' name]);
end

function addFS(name, sys, part)
global STEP massOf
add_block('sm_lib/Body Elements/File Solid', [sys '/' name]);
set_param([sys '/' name], 'ExtGeomFileName', STEP(part));
set_param([sys '/' name], 'BasedOnType', 'Mass');
set_param([sys '/' name], 'Mass', num2str(massOf(part), '%.12g'));
set_param([sys '/' name], 'MassUnits', 'kg');
end

function h = solP(sys, name)
ph = get_param([sys '/' name], 'PortHandles');
c = [ph.LConn ph.RConn]; h = c(1);
end
function h = outP(sys, name)
ph = get_param([sys '/' name], 'PortHandles');
c = [ph.RConn ph.LConn]; h = c(1);
end
function h = inP(sys, name)
ph = get_param([sys '/' name], 'PortHandles');
c = [ph.LConn ph.RConn]; h = c(1);
end

function pfun(fmt, varargin)
global fid
s = sprintf(fmt, varargin{:});
fprintf('%s', s); fprintf(fid, '%s', s);
end

