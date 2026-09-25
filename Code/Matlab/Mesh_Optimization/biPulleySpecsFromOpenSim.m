% Bi Pulley Specs From OpenSim
% Author: Ben Bolen
% Date: 2026-09-24
% Description: Parses the in-repo OpenSim robot-body model
% (Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_robotbody.osim, plain
% XML via matlab's xmlread -- no toolbox needed) and builds BiPam_pulley
% muscle-spec structs for a named list of muscles. Each spec carries the
% path points grouped by body frame in the BiPamData convention (rows before
% Cross(1) in the proximal frame, rows Cross(1)..Cross(2)-1 in the middle
% frame, rows Cross(2)..end in the distal frame), the two-hinge kinematic
% data needed to build the 4-D TransformationMat, and max_isometric_force.
% The driver (Opt_run_BiPulley) owns the transforms; this builder only
% extracts and documents the data.
%
% DOCUMENTED SIMPLIFICATIONS:
%   1. WRAPS REFUSED: stock wrap objects (PathWrapSet) are NOT approximated.
%      A wrap-carrying muscle errors out (insertWrapMidpoints) because the
%      chord-midpoint draft mixed body frames and appended the via past the
%      insertion. Supply wrap muscles with an explicit via point instead.
%      The default target muscles carry zero wraps and are unaffected.
%   2. HINGE JOINTS: each crossed OpenSim CustomJoint is reduced to a hinge
%      at its parent offset frame using the joint's FIRST coordinate's
%      rotation axis. The knee keeps its SimmSpline translation functions
%      (rolling knee; kin.transTheta/transX/transY); all other joints are
%      pure rotations. Secondary coordinates (hip adduction/rotation,
%      subtalar, mtp) are held at their default (zero) values.
%   3. FRAME FOLDS: the class computes moment arms about the CHILD-frame
%      origin, so every frame slot must be joint-centered. The hip and knee
%      child frames are already joint-centered in this model (their offset
%      frames are zero and the child body origin sits at the joint).
%      Calcaneus path points are re-expressed into an ANKLE-centered frame
%      by adding the subtalar offset (-0.04877, -0.04195, 0.00792), the
%      exact default-pose fold of the tibia->talus->calcn offset chain; this
%      reproduces the master's own VecTrans(T_a*T_s, ...) convention
%      (kin.pointFolds carries the offset; applied by the driver).
%   4. GASTROCNEMIUS: gait2392 has no muscle named "gastrocnemius"; the
%      default list maps it to med_gas_r (medial gastrocnemius, right).
%      Right-side muscles are used because the robot-body model carries the
%      real moved paths on the right leg (left-leg paths are stock/
%      placeholder in this model).
%   5. PATHS ENDING PROXIMALLY: robot-body variants of some muscles end
%      their path before the distal body (rect_fem_r ends on the femur).
%      Such a muscle gets a rigid distal joint (identity transform); if
%      only one row remains past the real crossing, a duplicate of the last
%      row is appended so the two-crossing bookkeeping stays well formed.
%      The rigid crossing's outputs are inert (zero-length segment) and
%      documented in spec.notes.
%   6. ANKLE ALIAS: the ankle CustomJoint connects tibia_r to talus_r, but
%      the muscles insert on calcn_r; the hinge data is aliased to the
%      tibia_r>calcn_r pair (subtalar frozen at zero, header note 3).
%
% TRANSMISSION CONFIG + BUNDLE PLACEMENT (Ben, 2026-09-24): the optional
% third argument carries the reverse-pulley transmission config
% (nPulleyBPA, tackleLineParts, routingMode; defaults 1 / 1 /
% 'moving_via') attached to each spec as spec.pulley and consumed by
% Opt_run_BiPulley / BiPam_pulley. When nPulleyBPA > 1 the builder also
% computes the BPA-bundle placement (spec.bundleOffsets, n x 3, expressed
% in the frame of the LAST path point): the bundle sits SYMMETRICALLY
% about the ORIGINAL OpenSim line -- offsets perpendicular to the line,
% spacing >= BPA diameter + margin (margin default 5 mm), bundle centroid
% ON the line -- so the first-order line of action and moment arm are
% preserved. DOCUMENTED COST: the footprint widens tangent to the local
% bone anatomy to footprintWidth = (n-1)*spacing + diameter (echoed in
% the batch summary CSV); ASYMMETRIC placement would shift the line of
% action and change the moment arm. spec.bundleMargin carries the margin
% so consumers can re-verify the spacing.
%
% Usage:
%   specs = biPulleySpecsFromOpenSim();     % default 5-muscle list
%   specs = biPulleySpecsFromOpenSim({'med_gas_r', 'soleus_r'});
%   specs = biPulleySpecsFromOpenSim({}, osimFile);
%   specs = biPulleySpecsFromOpenSim([], [], struct('nPulleyBPA', 2));

function specs = biPulleySpecsFromOpenSim(names, osimFile, pulleyConfig)

%% Inputs and model location
if nargin < 1 || isempty(names)
    % Keep the default list SHORT (per the batch contract): one per major
    % functional group. "Gastrocnemius" maps to med_gas_r (header note 4).
    names = {'med_gas_r', 'soleus_r', 'tib_ant_r', 'rect_fem_r', 'bifemsh_r'};
end
if nargin < 2 || isempty(osimFile)
    osimFile = fullfile(repoRoot(), 'Solid_Models', 'OpenSim', ...
        'Gait2392_Robotbody', 'gait2392_robotbody.osim');
end
if nargin < 3
    pulleyConfig = struct();
end
if ~exist(osimFile, 'file')
    error('biPulleySpecsFromOpenSim:missingModel', ...
        'OpenSim model not found at %s', osimFile)
end

%% Parse the model XML
doc = xmlread(osimFile);
model = doc.getElementsByTagName('Model').item(0);

kin = osimKinematics(model);

%% Body-depth map for slot assignment (right chain of gait2392)
% pelvis(0) -> femur_r(1) -> tibia_r(2) -> talus_r(3) -> calcn_r(4).
frameDepth = containers.Map( ...
    {'pelvis', 'femur_r', 'tibia_r', 'talus_r', 'calcn_r', 'toes_r', 'torso'}, ...
    {0, 1, 2, 3, 4, 5, 0});

%% Build one spec per requested name
muscleNodes = model.getElementsByTagName('Thelen2003Muscle');
nMuscle = muscleNodes.getLength();

specs = struct([]);   %#ok<STRNU>
for k = 1:numel(names)
    target = char(names{k});

    node = [];
    for m = 0:nMuscle - 1
        cand = muscleNodes.item(m);
        if strcmp(char(cand.getAttribute('name')), target)
            node = cand;
            break
        end
    end
    if isempty(node)
        error('biPulleySpecsFromOpenSim:unknownMuscle', ...
            'Muscle %s not found in %s', target, osimFile)
    end

    fmax = str2double(char(node.getElementsByTagName( ...
        'max_isometric_force').item(0).getFirstChild().getData()));

    pts0 = osimPathPoints(node);            % ordered {frameName, 1x3 loc}
    pts0 = insertWrapMidpoints(node, pts0); % wrap-as-via (header note 1)

    [points, cross, frames, rigidSlot2, notes] = ...
        groupSlots(target, pts0, frameDepth);

    spec = assembleSpec(target, points, cross, frames, rigidSlot2, ...
        fmax, notes, kin, pulleyConfig);
    spec.name = target;

    if isempty(specs)
        specs = spec;
    else
        specs(end + 1) = spec;   %#ok<AGROW>
    end
end

end

%% ---------------------------------------------------------------------
function rootDir = repoRoot()
% Self-locate the Bipedal_Robot repo root from this file's location.
rootDir = fileparts(mfilename('fullpath'));
for k = 1:8
    [parent, name] = fileparts(rootDir);
    if strcmpi(name, 'Bipedal_Robot')
        return
    end
    if strcmp(parent, rootDir)
        error('biPulleySpecsFromOpenSim:noRepo', ...
            'Could not locate the Bipedal_Robot repo root')
    end
    rootDir = parent;
end
end

%% ---------------------------------------------------------------------
function val = cfgField(cfg, name, default)
% Config field with default (used for the transmission config).
if isfield(cfg, name)
    val = cfg.(name);
else
    val = default;
end
end

%% ---------------------------------------------------------------------
function childText = childText(node, tagName)
% First child text of the named tag, or '' when absent/empty (the osim
% carries plenty of empty elements, e.g. <coordinates></coordinates>).
childText = '';
nodes = node.getElementsByTagName(tagName);
if nodes.getLength() == 0
    return
end
n0 = nodes.item(0);
if ~n0.hasChildNodes()
    return
end
childText = strtrim(char(n0.getFirstChild().getData()));
end

%% ---------------------------------------------------------------------
function pts = osimPathPoints(muscleNode)
% Ordered {frameName, 1x3 location} rows from the GeometryPath
% (first = origin, last = insertion, between = via).
pts = {};
gp = muscleNode.getElementsByTagName('GeometryPath').item(0);
ppNodes = gp.getElementsByTagName('PathPoint');
n = ppNodes.getLength();
pts = cell(n, 2);
for i = 0:n - 1
    p = ppNodes.item(i);
    sock = childText(p, 'socket_parent_frame');
    pts{i + 1, 1} = regexprep(sock, '.*/', '');   % strip /bodyset/
    pts{i + 1, 2} = sscanf(childText(p, 'location'), '%f')';
end
if n < 2
    error('biPulleySpecsFromOpenSim:shortPath', ...
        'Path of %s has fewer than two points', ...
        char(muscleNode.getAttribute('name')))
end
end

%% ---------------------------------------------------------------------
function pts = insertWrapMidpoints(muscleNode, pts)
% Wrap handling (header note 1): a real PathWrap needs its wrap surface
% modeled, not approximated by an appended via. The earlier chord-midpoint
% draft mixed body frames (the last two path points can live in different
% frames) and landed the synthetic via PAST the insertion, silently moving
% the route end. Wrap-carrying muscles are therefore refused outright:
% model them with an explicit via point, or extend this builder with real
% wrap-surface geometry. Muscles without wraps pass through unchanged.
gp = muscleNode.getElementsByTagName('GeometryPath').item(0);
wrapNodes = gp.getElementsByTagName('PathWrap');
nWrap = wrapNodes.getLength();
if nWrap == 0
    return
end
error('biPulleySpecsFromOpenSim:WrapNotSupported', ...
    ['Muscle %s carries %d PathWrap(s). This builder does not ' ...
     'approximate wraps: supply the muscle with an explicit via point ' ...
     'at the wrap, or extend biPulleySpecsFromOpenSim with wrap-surface ' ...
     'geometry first.'], ...
    char(muscleNode.getAttribute('name')), nWrap)
end

%% ---------------------------------------------------------------------
function [points, cross, frames, rigidSlot2, notes] = groupSlots(name, pts, frameDepth)
% Assign path points to up to three frame slots ordered by body depth.
% Slot 1 = proximal frame, slot 2 = middle frame, slot 3 = distal frame.
% When the found frames span an intermediate body (med_gas: femur -> calcn
% skips the tibia), an EMPTY middle slot is kept and both crossings fall on
% the first slot-3 row (BiPamData duplicate-crossing convention; the class
% chains the two transforms child-to-parent in that case).
chainBodies = {'pelvis', 'femur_r', 'tibia_r', 'talus_r', 'calcn_r', 'toes_r'};

n = size(pts, 1);
depth = zeros(n, 1);
for i = 1:n
    key = pts{i, 1};
    if ~isKey(frameDepth, key)
        error('biPulleySpecsFromOpenSim:unknownFrame', ...
            '%s: path point on unhandled body %s', name, key)
    end
    depth(i) = frameDepth(key);
end

ud = unique(depth);
if numel(ud) < 2
    error('biPulleySpecsFromOpenSim:singleFrame', ...
        '%s: path crosses no joint (all points on one body)', name)
end

% Decide the three slot depth levels.
if numel(ud) >= 3
    slotDepths = ud(1:3);
    emptyMiddle = false;
elseif ud(2) - ud(1) > 1 && ud(1) <= 1
    % An intermediate body is crossed without path points on it AND the
    % proximal frame is pelvis/femur (the knee lies inside the span):
    % keep it as an empty middle slot so both joints stay represented
    % (med_gas: femur -> [empty tibia] -> calcn, duplicate crossing).
    slotDepths = [ud(1), ud(1) + 1, ud(2)];
    emptyMiddle = true;
else
    % Tibia-proximal spans (soleus/tib_ant: tibia -> calcn) get the ankle
    % alias at crossing 1 and a rigid clone at crossing 2 (header note 5).
    slotDepths = [ud(1), ud(2), ud(2)];
    emptyMiddle = false;
end

% Slot index per row.
slot = zeros(n, 1);
points = zeros(n, 3);
for i = 1:n
    s = find(depth(i) == slotDepths, 1);
    if isempty(s)
        error('biPulleySpecsFromOpenSim:depthJump', ...
            '%s: point depth %d outside slots %s', ...
            name, depth(i), mat2str(slotDepths))
    end
    slot(i) = s;
    points(i, :) = pts{i, 2};
end

% Frame labels per slot (first body seen in the slot, else the chain name).
frames = cell(1, 3);
for s = 1:3
    hit = find(slot == s, 1);
    if ~isempty(hit)
        frames{s} = pts{hit, 1};
    else
        frames{s} = chainBodies{slotDepths(s) + 1};
        if s == 2
            frames{s} = [frames{s}, ' (no path points)'];
        end
    end
end

% Cross rows: first row of slot 2, first row of slot 3 (BiPamData rule).
% Slot 2 can be empty (duplicate-crossing case); firstRow guards find().
cross = zeros(1, 2);
cross(1) = firstRow(slot == 2);
cross(2) = firstRow(slot == 3);
rigidSlot2 = false;
notes = {};

if emptyMiddle
    % Duplicate-crossing convention: both crossings at the slot-3 row.
    cross(1) = cross(2);
    notes{end + 1} = ...   %#ok<AGROW>
        'empty middle frame: duplicate crossing, transforms chained';
elseif nSlotOrZero(slot) == 2
    % One real crossing. Crossing 1 = slot1->slot2 (the real joint);
    % crossing 2 is made rigid (header note 5).
    rigidSlot2 = true;
    if cross(1) < n
        cross(2) = cross(1) + 1;
        frames{3} = [frames{2}, ' (rigid)'];
        notes{end + 1} = ...   %#ok<AGROW>
            'slot 3 = rigid clone; crossing-2 outputs inert';
    else
        % Single row past the crossing: append a duplicate of it.
        points(end + 1, :) = points(end, :);
        slot(end + 1) = 3;
        cross(2) = n + 1;
        frames{3} = [frames{2}, ' (rigid, duplicate row)'];
        notes{end + 1} = ...   %#ok<AGROW>
            'duplicate insertion row appended; crossing-2 outputs inert';
    end
end
end

%% ---------------------------------------------------------------------
function nSlots = nSlotOrZero(slot)
nSlots = numel(unique(slot));
end

%% ---------------------------------------------------------------------
function idx = firstRow(mask)
% First true index of a logical mask, or 0 when none.
hit = find(mask, 1);
if isempty(hit)
    idx = 0;
else
    idx = hit;
end
end

%% ---------------------------------------------------------------------
function spec = assembleSpec(name, points, crossRows, frames, rigidSlot2, ...
    fmax, notes, kin, pulleyConfig)
% Attach the two-hinge kinematics to the grouped points. Joint i drives
% crossing i: joint 1 = slot1->slot2, joint 2 = slot2->slot3. (The crossing
% argument is named crossRows so MATLAB's cross() stays callable below.)
spec = struct();
spec.name = name;
spec.source = 'opensim';
spec.points = points;
spec.cross = crossRows;
spec.frames = frames;
spec.diameter = 20;

% Transmission config (defaults 1 / 1 / 'moving_via'); passed through to
% BiPam_pulley's 15th argument by the driver.
spec.pulley = struct( ...
    'nPulleyBPA', cfgField(pulleyConfig, 'nPulleyBPA', 1), ...
    'tackleLineParts', cfgField(pulleyConfig, 'tackleLineParts', 1), ...
    'routingMode', cfgField(pulleyConfig, 'routingMode', 'moving_via'));

% Symmetric BPA-bundle placement when nPulleyBPA > 1 (see header):
% offsets perpendicular to the ORIGINAL line, spacing >= diameter +
% margin, centroid ON the line, first-order line of action preserved.
nBPA = spec.pulley.nPulleyBPA;
bundleMargin = 0.005;   % m between adjacent BPA skins
spacing = spec.diameter * 1e-3 + bundleMargin;
if nBPA > 1
    % Original-line direction: the last NON-ZERO-LENGTH path segment
    % (walking back over the rigid-clone duplicate rows, header note 5).
    u = [];
    for r = size(points, 1) - 1:-1:1
        cand = points(r + 1, :) - points(r, :);
        if norm(cand) > 1e-12
            u = cand / norm(cand);
            break
        end
    end
    if isempty(u)
        error('biPulleySpecsFromOpenSim:degeneratePath', ...
            '%s: every path segment is zero-length; no bundle line', name)
    end
    % A perpendicular to the line (any unit vector orthogonal to u).
    if abs(u(3)) < 0.9
        perp = cross(u, [0, 0, 1]);
    else
        perp = cross(u, [0, 1, 0]);
    end
    perp = perp / norm(perp);
    offs = ((1:nBPA)' - (nBPA + 1) / 2) * spacing;   % centroid at zero
    spec.bundleOffsets = offs .* perp;               % n x 3, distal frame
    spec.bundleMargin = bundleMargin;
    spec.footprintWidth = (nBPA - 1) * spacing + spec.diameter * 1e-3;
    notes{end + 1} = sprintf( ...
        ['%d-BPA bundle symmetric about the original line (spacing ' ...
         '%.1f mm, footprint %.1f mm); asymmetric placement would shift ' ...
         'the line of action'], nBPA, 1000 * spacing, ...
        1000 * spec.footprintWidth); %#ok<AGROW>
else
    spec.bundleOffsets = zeros(0, 3);
    spec.bundleMargin = bundleMargin;
    spec.footprintWidth = spec.diameter * 1e-3;
end

spec.fmax = fmax;
spec.fmaxSource = sprintf('gait2392_robotbody.osim max_isometric_force (%g N)', fmax);

body1 = cleanBody(frames{1});
body2 = cleanBody(frames{2});
body3 = cleanBody(frames{3});

% Build field-by-field: struct() with cell values would create a struct
% ARRAY, not the scalar we need.
joints = struct();
joints.pivots = zeros(2, 3);
joints.axes = zeros(2, 3);
joints.types = {'fixed', 'fixed'};
joints.transTheta = [];
joints.transX = [];
joints.transY = [];
joints.pointFolds = {[0; 0; 0], [0; 0; 0]};
joints.kneeSlot = 0;                    % set below when the knee is present
joints.thetaRangesDeg = [0, 0; 0, 0];   % grid sweep ranges per crossing

% Joint 1: between body1 and body2 (always a real pair).
[joints, note1] = applyJoint(joints, 1, body1, body2, kin);

% Joint 2: between body2 and body3 (rigid clone case: skip).
if rigidSlot2 || strcmp(body2, body3)
    note2 = 'crossing 2 rigid (identity transform)';
else
    [joints, note2] = applyJoint(joints, 2, body2, body3, kin);
end

% Which joint is the knee (for the driver's rolling-knee pivot override):
% the knee joint is the femur_r<->tibia_r pair.
if strcmp(body1, 'femur_r') && strcmp(body2, 'tibia_r')
    joints.kneeSlot = 1;    % knee at crossing 1 (med_gas duplicate chain)
elseif strcmp(body2, 'femur_r') && strcmp(body3, 'tibia_r')
    joints.kneeSlot = 2;    % knee at crossing 2
else
    joints.kneeSlot = 0;    % no knee in this muscle's chain
end

% Fold calcaneus points into ankle-centered coordinates (header note 3).
if strcmp(body3, 'calcn_r') && ~strcmp(body2, 'calcn_r')
    joints.pointFolds{2} = kin.subtalarOffset;
    notes{end + 1} = 'calcn points shifted by the subtalar offset'; ...
        %#ok<AGROW>
end

    spec.kin = joints;
    spec.notes = sprintf('osim route %s: %d points, cross = [%d %d]; %s', ...
    name, size(points, 1), crossRows(1), crossRows(2), ...
    strjoin([notes, {note1, note2}], '; '));
end

%% ---------------------------------------------------------------------
function body = cleanBody(frameLabel)
body = regexprep(strtrim(frameLabel), '\s*\(.*$', '');   % strip "(rigid...)"
end

%% ---------------------------------------------------------------------
function [joints, note] = applyJoint(joints, slot, bodyParent, bodyChild, kin)
% Look up the hinge data for the parent>child joint parsed from JointSet.
key = sprintf('%s_to_%s', bodyParent, bodyChild);
if isKey(kin.joints, key)
    J = kin.joints(key);
    joints.types{slot} = 'hinge';
    joints.pivots(slot, :) = J.pivot;
    joints.axes(slot, :) = J.axis;
    joints.thetaRangesDeg(slot, :) = jointRange(J.name);
    if ~isempty(J.transTheta)
        joints.transTheta = J.transTheta;
        joints.transX = J.transX;
        joints.transY = J.transY;
    end
    note = sprintf('crossing %d = %s hinge (axis [%g %g %g], pivot [%g %g %g])', ...
        slot, J.name, J.axis, J.pivot);
else
    joints.types{slot} = 'fixed';
    note = sprintf('crossing %d: no hinge for %s, held rigid', slot, key);
end
end

%% ---------------------------------------------------------------------
function range = jointRange(jname)
% Sweep ranges (deg) from the stock gait2392 coordinate ranges (the model
% carries [-120,10] for the knee, etc.; the plantar/dorsi ankle range is
% narrowed to the functional -50..20 window).
switch jname
    case 'hip_r'
        range = [-15, 110];
    case 'knee_r'
        range = [-120, 10];
    case 'ankle_r'
        range = [-50, 20];
    case 'subtalar_r'
        range = [-25, 35];
    otherwise
        range = [0, 0];
end
end

%% ---------------------------------------------------------------------
function kin = osimKinematics(model)
% Parse JointSet CustomJoints into hinge approximations (header note 2).
% Pivots are in the PARENT body frame (matching RpToTrans child->parent).
% kin.joints is a containers.Map keyed '<parentBody>_to_<childBody>'.
kin.joints = containers.Map('KeyType', 'char', 'ValueType', 'any');
kin.subtalarOffset = [-0.04877; -0.04195; 0.00792];  % talus->calcn offset

js = model.getElementsByTagName('JointSet').item(0);
jNodes = js.getElementsByTagName('CustomJoint');
for i = 0:jNodes.getLength() - 1
    j = jNodes.item(i);
    jname = char(j.getAttribute('name'));

    parentOff = offsetFrame(j, 'socket_parent_frame');
    childOff = offsetFrame(j, 'socket_child_frame');
    if isempty(parentOff.frame) || isempty(childOff.frame)
        continue
    end

    % First coordinate with a non-empty axis line drives the hinge.
    coordName = '';
    axis = [0, 0, 1];
    taNodes = j.getElementsByTagName('TransformAxis');
    nTA = taNodes.getLength();
    for t = 0:nTA - 1
        ta = taNodes.item(t);
        cn = childText(ta, 'coordinates');
        if isempty(cn)
            continue
        end
        axis = sscanf(childText(ta, 'axis'), '%f')';
        coordName = regexprep(cn, ',', ' ');   % multi-DoF: first wins
        break
    end

    % Translation functions bound to the same coordinate (rolling knee).
    transTheta = [];
    transX = [];
    transY = [];
    for t = 0:nTA - 1
        ta = taNodes.item(t);
        cn = childText(ta, 'coordinates');
        if isempty(cn) || ~strcmp(regexprep(cn, ',', ' '), coordName)
            continue
        end
        ax = sscanf(childText(ta, 'axis'), '%f')';
        fn = ta.getElementsByTagName('SimmSpline');
        if fn.getLength() == 0
            continue
        end
        xy = sscanf(char(fn.item(0).getFirstChild().getData()), '%f')';
        if mod(numel(xy), 2) ~= 0
            continue
        end
        xy = reshape(xy, 2, []).';
        if all(abs(ax - [1, 0, 0]) < 1e-9)
            transTheta = xy(:, 1)';
            transX = xy(:, 2)';
        elseif all(abs(ax - [0, 1, 0]) < 1e-9)
            transTheta = xy(:, 1)';
            transY = xy(:, 2)';
        end
    end

    parentBody = regexprep(parentOff.frame, '_offset$', '');
    childBody = regexprep(childOff.frame, '_offset$', '');
    J = struct('pivot', parentOff.t, 'axis', axis, ...
        'transTheta', transTheta, 'transX', transX, 'transY', transY, ...
        'name', jname);
    kin.joints(sprintf('%s_to_%s', parentBody, childBody)) = J;

    % Ankle alias (header note 6): tibia_r>talus_r data also serves the
    % tibia_r>calcn_r pair (subtalar frozen at zero).
    if strcmp(jname, 'ankle_r')
        kin.joints('tibia_r_to_calcn_r') = J;
    end
end
end

%% ---------------------------------------------------------------------
function off = offsetFrame(jointNode, socketTag)
% Resolve a joint socket to its PhysicalOffsetFrame (frame name + t).
off = struct('frame', '', 't', [0; 0; 0]);
socks = jointNode.getElementsByTagName(socketTag);
if socks.getLength() == 0
    return
end
target = regexprep(strtrim(char(socks.item(0).getFirstChild().getData())), ...
    '.*/', '');

frames = jointNode.getElementsByTagName('PhysicalOffsetFrame');
for i = 0:frames.getLength() - 1
    f = frames.item(i);
    if strcmp(char(f.getAttribute('name')), target)
        off.frame = target;
        off.t = sscanf(childText(f, 'translation'), '%f')';
        return
    end
end
end
