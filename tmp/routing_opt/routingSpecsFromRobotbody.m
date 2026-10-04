function [specs, report, chain] = routingSpecsFromRobotbody(mapStruct, osimFile)
% Routing Specs From Robotbody
% Campaign: routing27 (2026-10-03)
% Description: Builds BiPam_pulley muscle-spec structs for the 27-actuator
% routing campaign from gait2392_robotbody.osim (the post-09-26 robot-route
% model of record; gait2327.osim and the mjc trees are one path-point stale
% and must not be used). Adapted from biPulleySpecsFromOpenSim.m (the
% committed 5-muscle builder) with these campaign-specific extensions,
% every one documented in the returned report:
%
%   1. FULL RIGHT CHAIN incl. torso and toes. Frame depths: pelvis 0,
%      femur_r 1, tibia_r 2, talus_r 3, calcn_r 4, toes_r 5, torso 6.
%   2. MOVING / CONDITIONAL PATH POINTS: the committed builder silently
%      dropped them (it reads only 'PathPoint' tags; rect_fem_r lost its
%      tibial point that way). Here MovingPathPoints are frozen at their
%      coordinate's DEFAULT value by evaluating their SimmSpline knots
%      (pchip), and ConditionalPathPoints are included only when active
%      at the default coordinate value (range check in radians). A probe
%      (tmp/routing_opt/probe_robot_routes.py, run 2026-10-03) verified
%      none of the 27 seed muscles carries a PathWrap, so wraps remain a
%      hard error if ever encountered.
%   3. TRUNK: pelvis->torso routes cross the back joint (hinge on
%      lumbar_extension, pivot at the pelvis_offset translation). Psoas'
%      route starts on the torso BEFORE crossing to the pelvis; those
%      torso rows are folded into the pelvis frame through the back joint
%      at its default pose (lumbar frozen; documented simplification).
%   4. TIBIA->CALCN SPANS: the committed builder aliased tibia->calcn to
%      the ankle hinge with a rigid clone, which cannot express subtalar
%      torque (tib_post and the peronei are subtalar-primary). Here the
%      empty TALUS middle frame is inserted instead, so crossing 1 is the
%      true ankle (tibia->talus) and crossing 2 the true subtalar
%      (talus->calcn), no offset fold needed. Femur->calcn spans (med/lat
%      gas) keep the committed empty-tibia middle + ankle alias + calcn
%      subtalar-offset fold, unchanged from the committed builder.
%   5. EXPLICIT COORDINATE SELECTION: each crossing's hinge axis is the
%      TransformAxis of the campaign coordinate (hip_flexion_r,
%      knee_angle_r, ankle_angle_r, subtalar_angle_r, mtp_angle_r,
%      lumbar_extension), not whatever axis comes first in the XML.
%   6. SEED PREFERENCE with fallback: each actuator lists candidate
%      robotbody muscles; the first that resolves to a well-formed spec
%      wins. appliesForce flags are ignored (they belong to Connor's
%      analysis, not the routing; the campaign re-enables what it uses).
%   7. targetCrossing: the crossing whose coordinate equals the
%      actuator's primary DOF, verified at build time.
%
% Usage (campaign driver):
%   [specs, report] = routingSpecsFromRobotbody(mapStruct)   % default model
%   [specs, report] = routingSpecsFromRobotbody(mapStruct, osimFile)

if nargin < 2 || isempty(osimFile)
    osimFile = fullfile(repoRoot(), 'Solid_Models', 'OpenSim', ...
        'Gait2392_Robotbody', 'gait2392_robotbody.osim');
end
if ~exist(osimFile, 'file')
    error('routingSpecsFromRobotbody:missingModel', ...
        'OpenSim model not found at %s', osimFile)
end

doc = xmlread(osimFile);
model = doc.getElementsByTagName('Model').item(0);
kin = osimKinematics(model);
coords = osimCoordinates(model);

frameDepth = containers.Map( ...
    {'pelvis', 'femur_r', 'tibia_r', 'talus_r', 'calcn_r', 'toes_r', 'torso'}, ...
    {0, 1, 2, 3, 4, 5, 6});

muscleNodes = model.getElementsByTagName('Thelen2003Muscle');
muscleByName = containers.Map('KeyType', 'char', 'ValueType', 'any');
for m = 0:muscleNodes.getLength() - 1
    cand = muscleNodes.item(m);
    muscleByName(char(cand.getAttribute('name'))) = cand;
end

specs = struct([]);
report = struct('actuatorId', {}, 'masterName', {}, 'seedTried', {}, ...
    'seedChosen', {}, 'note', {});

for k = 1:numel(mapStruct.actuators)
    a = mapStruct.actuators(k);
    chosen = [];
    tries = {};
    for s = 1:numel(a.robot_seed_preference)
        cname = char(a.robot_seed_preference(s));
        try
            node = muscleByName(cname);   % errors if unknown muscle
            spec = buildOneSpec(a, node, frameDepth, kin, coords, mapStruct);
            chosen = spec;
            tries{end + 1} = sprintf('%s: OK', cname); %#ok<AGROW>
            break
        catch err
            tries{end + 1} = sprintf('%s: %s', cname, err.message); %#ok<AGROW>
        end
    end
    rep = struct();
    rep.actuatorId = a.id;
    rep.masterName = char(a.master_name);
    rep.seedTried = tries;
    if isempty(chosen)
        rep.seedChosen = '';
        rep.note = 'ALL SEEDS FAILED';
        error('routingSpecsFromRobotbody:noSeed', ...
            'Actuator %d (%s): every seed failed: %s', a.id, ...
            char(a.master_name), strjoin(tries, ' | '))
    end
    rep.seedChosen = chosen.name;
    rep.note = chosen.notes;
    if isempty(specs)
        specs = chosen;
    else
        specs(end + 1) = chosen; %#ok<AGROW>
    end
    report(end + 1) = rep; %#ok<AGROW>
end

% Global chain for cross-actuator geometry (bone clouds, pelvis-frame
% segment mapping): every chain edge with its hinge data and its
% DEFAULT-pose transform. Built once from the parsed JointSet.
edgeDefs = { ...
    'pelvis>femur_r', 'pelvis_to_femur_r', 'hip_flexion_r'; ...
    'femur_r>tibia_r', 'femur_r_to_tibia_r', 'knee_angle_r'; ...
    'tibia_r>talus_r', 'tibia_r_to_talus_r', 'ankle_angle_r'; ...
    'talus_r>calcn_r', 'talus_r_to_calcn_r', 'subtalar_angle_r'; ...
    'calcn_r>toes_r', 'calcn_r_to_toes_r', 'mtp_angle_r'; ...
    'pelvis>torso', 'pelvis_to_torso', 'lumbar_extension'};
chain = struct();
chain.edges = containers.Map('KeyType', 'char', 'ValueType', 'any');
for e = 1:size(edgeDefs, 1)
    key = edgeDefs{e, 1};
    jkey = edgeDefs{e, 2};
    cn = edgeDefs{e, 3};
    if ~isKey(kin.joints, jkey) || ~isKey(kin.joints(jkey).axes, cn)
        error('routingSpecsFromRobotbody:noEdge', ...
            'Chain edge %s (%s on %s) missing from the model', key, cn, jkey)
    end
    J = kin.joints(jkey);
    A = J.axes(cn);
    edge = struct('axis', A.axis, 'pivot', J.pivot, 'name', J.name, ...
        'transThetaX', [], 'transX', [], 'transThetaY', [], ...
        'transY', []);
    if isfield(A, 'transThetaX')
        edge.transThetaX = A.transThetaX;
        edge.transX = A.transX;
    end
    if isfield(A, 'transThetaY')
        edge.transThetaY = A.transThetaY;
        edge.transY = A.transY;
    end
    th0 = 0;
    if isKey(coords, cn)
        th0 = coords(cn).default_rad;
    end
    edge.Tdefault = routingChainLib('edgeTrans', edge, th0, true);
    chain.edges(key) = edge;
end
chain.coords = coords;
end

%% ---------------------------------------------------------------------
function rootDir = repoRoot()
rootDir = fileparts(mfilename('fullpath'));
for k = 1:8
    [parent, name] = fileparts(rootDir);
    if strcmpi(name, 'Bipedal_Robot')
        return
    end
    if strcmp(parent, rootDir)
        error('routingSpecsFromRobotbody:noRepo', ...
            'Could not locate the Bipedal_Robot repo root')
    end
    rootDir = parent;
end
end

%% ---------------------------------------------------------------------
function childText = childText(node, tagName)
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
function c = osimCoordinates(model)
% Coordinate name -> struct(default_rad, range_deg). OpenSim 4.x nests
% <Coordinate> elements inside the joints (no model-level CoordinateSet).
c = containers.Map('KeyType', 'char', 'ValueType', 'any');
js = model.getElementsByTagName('JointSet').item(0);
nodes = js.getElementsByTagName('Coordinate');
for i = 0:nodes.getLength() - 1
    n = nodes.item(i);
    nm = char(n.getAttribute('name'));
    s = struct();
    d = str2double(childText(n, 'default_value'));
    if isnan(d)
        d = 0;
    end
    s.default_rad = d;
    lo = str2double(childText(n, 'range_min'));
    hi = str2double(childText(n, 'range_max'));
    if ~isnan(lo) && ~isnan(hi)
        s.range_deg = [lo, hi] * 180 / pi;
    else
        s.range_deg = [0, 0];
    end
    c(nm) = s;
end
end

%% ---------------------------------------------------------------------
function pts = osimPathPoints(node, coords)
% Ordered {frameName, 1x3 loc, pointType} rows. PathPoint = fixed;
% MovingPathPoint frozen at the coordinate default (spline knots, pchip);
% ConditionalPathPoint kept only when active at the coordinate default.
pts = {};
gp = node.getElementsByTagName('GeometryPath').item(0);
pps = gp.getElementsByTagName('PathPointSet').item(0);
% OpenSim 4.x nests the points under PathPointSet/objects; fall back to
% direct children for flat variants.
holder = pps;
objs = pps.getElementsByTagName('objects');
if objs.getLength() > 0
    holder = objs.item(0);
end
kids = holder.getChildNodes();
for i = 0:kids.getLength() - 1
    p = kids.item(i);
    tag = char(p.getNodeName());
    if ~any(strcmp(tag, {'PathPoint', 'MovingPathPoint', ...
            'ConditionalPathPoint'}))
        continue
    end
    sock = childText(p, 'socket_parent_frame');
    frame = regexprep(sock, '.*/', '');
    switch tag
        case 'PathPoint'
            loc = sscanf(childText(p, 'location'), '%f')';
            typ = 'fixed';
        case 'MovingPathPoint'
            loc = zeros(1, 3);
            for c = 1:3
                ax = char('x' + c - 1);
                sockc = childText(p, sprintf('socket_%s_coordinate', ax));
                cn = regexprep(sockc, '.*/', '');
                if ~isKey(coords, cn)
                    error('routingSpecsFromRobotbody:movingCoord', ...
                        'Unknown coordinate %s on a moving point', cn)
                end
                theta0 = coords(cn).default_rad;
                loc(c) = splineEval(p, sprintf('%s_location', ax), theta0);
            end
            typ = 'moving-at-default';
        case 'ConditionalPathPoint'
            rng = sscanf(childText(p, 'range'), '%f');
            sockc = childText(p, 'socket_coordinate');
            cn = regexprep(sockc, '.*/', '');
            if ~isKey(coords, cn)
                error('routingSpecsFromRobotbody:condCoord', ...
                    'Unknown coordinate %s on a conditional point', cn)
            end
            theta0 = coords(cn).default_rad;
            if theta0 < rng(1) || theta0 > rng(2)
                continue    % inactive at default: dropped (documented)
            end
            loc = sscanf(childText(p, 'location'), '%f')';
            typ = 'conditional-active';
    end
    pts(end + 1, :) = {frame, loc, typ}; %#ok<AGROW>
end
gpTags = gp.getElementsByTagName('PathWrap');
if gpTags.getLength() > 0
    error('routingSpecsFromRobotbody:wrap', ...
        'Muscle %s carries %d PathWrap(s); the probe found none on the 27 seeds, so this is unexpected. Model the wrap explicitly before including this muscle.', ...
        char(node.getAttribute('name')), gpTags.getLength())
end
if size(pts, 1) < 2
    error('routingSpecsFromRobotbody:shortPath', ...
        'Path of %s has fewer than two usable points', ...
        char(node.getAttribute('name')))
end
end

%% ---------------------------------------------------------------------
function v = splineEval(pointNode, locTag, theta)
% Evaluate a SimmSpline knot list (or Constant) at theta. SimmSpline
% stores its knots in explicit <x> and <y> children (reading the first
% child alone yields only the x knots; the committed builder's odd-count
% guard silently dropped the rolling-knee splines that way).
fn = pointNode.getElementsByTagName(locTag).item(0);
sp = fn.getElementsByTagName('SimmSpline');
if sp.getLength() > 0
    n0 = sp.item(0);
    xs = sscanf(childText(n0, 'x'), '%f')';
    ys = sscanf(childText(n0, 'y'), '%f')';
    if numel(xs) ~= numel(ys) || numel(xs) < 2
        error('routingSpecsFromRobotbody:badSpline', ...
            'Malformed SimmSpline in %s', locTag)
    end
    v = interp1(xs, ys, theta, 'pchip', 'extrap');
else
    v = str2double(char(fn.getFirstChild().getData()));
    if isnan(v)
        v = 0;
    end
end
end

%% ---------------------------------------------------------------------
function spec = buildOneSpec(a, node, frameDepth, kin, coords, mapStruct)
name = char(node.getAttribute('name'));
fmax = str2double(childText(node, 'max_isometric_force'));

pts = osimPathPoints(node, coords);

% Torso-before-pelvis fold (psoas): torso rows occurring before any
% pelvis/below row are folded into the pelvis frame through the back
% joint at its default pose (rotation = default angle, translation =
% child-offset chain). At lumbar default 0 the rotation is identity, so
% the fold is the pure offset sum (documented simplification).
folded = false;
depths = zeros(size(pts, 1), 1);
for i = 1:size(pts, 1)
    key = pts{i, 1};
    if ~isKey(frameDepth, key)
        error('routingSpecsFromRobotbody:unknownFrame', ...
            '%s: path point on unhandled body %s', name, key)
    end
    depths(i) = frameDepth(key);
end
if any(depths == 6) && depths(1) == 6 && any(depths < 6)
    Tfold = kin.joints('pelvis_to_torso');
    for i = 1:size(pts, 1)
        if depths(i) == 6
            th0 = coords('lumbar_extension').default_rad;
            R0 = axisRot(Tfold.axes('lumbar_extension').axis, th0);
            pts{i, 2} = (R0 * pts{i, 2}.' + Tfold.pivot.').';
            pts{i, 1} = 'pelvis';
            depths(i) = 0;
        end
    end
    folded = true;
end

[points, cross, frames, emptyMiddle, rigidSlot2, notes] = ...
    groupSlots(name, pts, depths);

spec = struct();
spec.name = name;
spec.masterName = char(a.master_name);
spec.actuatorId = a.id;
spec.source = 'robotbody';
spec.points = points;
spec.cross = cross;
spec.frames = frames;
spec.diameter = 20;
spec.humanGroup = a.human_group;
spec.primaryDof = char(a.primary_dof);
spec.fmax = fmax;
spec.fmaxSource = sprintf('gait2392_robotbody.osim max_isometric_force (%g N)', fmax);
if folded
    notes{end + 1} = 'torso rows folded to pelvis via back joint at default'; %#ok<AGROW>
end

% Kinematics + explicit coordinate per crossing.
joints = struct();
joints.pivots = zeros(2, 3);
joints.axes = zeros(2, 3);
joints.types = {'fixed', 'fixed'};
joints.transThetaX = [];
joints.transX = [];
joints.transThetaY = [];
joints.transY = [];
joints.kneeSlot = 0;
joints.thetaRangesDeg = [0, 0; 0, 0];
joints.coordNames = {'', ''};

chainCoord = chainCoordinateFor(cleanBody(frames{1}), cleanBody(frames{2}), ...
    emptyMiddle, rigidSlot2, char(a.primary_dof));
[joints, note1] = applyJoint(joints, 1, cleanBody(frames{1}), ...
    cleanBody(frames{2}), chainCoord, kin, coords, mapStruct);
joints.coordNames{1} = chainCoord;
if rigidSlot2
    note2 = 'crossing 2 rigid (identity transform)';
    joints.coordNames{2} = '';
else
    chainCoord2 = chainCoordinateFor2(cleanBody(frames{2}), ...
        cleanBody(frames{3}), char(a.primary_dof), joints.coordNames{1});
    [joints, note2] = applyJoint(joints, 2, cleanBody(frames{2}), ...
        cleanBody(frames{3}), chainCoord2, kin, coords, mapStruct, false, false);
    joints.coordNames{2} = chainCoord2;
end

% Which crossing carries the primary DOF (verified).
spec.targetCrossing = find(strcmp(joints.coordNames, char(a.primary_dof)), 1);
if isempty(spec.targetCrossing)
    spec.targetCrossing = 1;
    notes{end + 1} = sprintf(['primary DOF %s not on a crossing (frozen ' ...
        'at default); target crossing 1'], char(a.primary_dof)); %#ok<AGROW>
end

% Calcn fold for the femur->calcn alias geometry (committed convention).
if strcmp(cleanBody(frames{3}), 'calcn_r') && ...
        strcmp(cleanBody(frames{2}), 'tibia_r') && emptyMiddle
    spec.points(cross(2):end, :) = spec.points(cross(2):end, :) ...
        + kin.subtalarOffset.';
    notes{end + 1} = 'calcn points shifted by the subtalar offset (alias)'; %#ok<AGROW>
end

spec.kin = joints;
spec.notes = sprintf(['robotbody route %s (actuator %d %s): %d points, ' ...
    'cross = [%d %d], primary %s at crossing %d; %s'], name, a.id, ...
    char(a.master_name), size(points, 1), cross(1), cross(2), ...
    char(a.primary_dof), spec.targetCrossing, ...
    strjoin([notes, {note1, note2}], '; '));

% Sanity: nonzero RoM on the primary crossing.
pr = joints.thetaRangesDeg(spec.targetCrossing, :);
if pr(2) <= pr(1)
    error('routingSpecsFromRobotbody:zeroRoM', ...
        '%s: zero RoM on the primary crossing', name)
end
end

%% ---------------------------------------------------------------------
function body = cleanBody(frameLabel)
body = regexprep(strtrim(frameLabel), '\s*\(.*$', '');
end

%% ---------------------------------------------------------------------
function cc = chainCoordinateFor(body1, body2, emptyMiddle, rigidSlot2, primary)
% Coordinate of crossing 1 for the body pair (campaign chain table).
if rigidSlot2
    cc = primary;   % single real crossing: force it to the primary DOF's joint when the pair matches
    return
end
key = sprintf('%s>%s', body1, body2);
switch key
    case 'pelvis>torso'
        cc = 'lumbar_extension';
    case 'pelvis>femur_r'
        cc = 'hip_flexion_r';
    case 'femur_r>tibia_r'
        cc = 'knee_angle_r';
    case 'tibia_r>talus_r'
        cc = 'ankle_angle_r';
    case 'tibia_r>calcn_r'   % empty-tibia alias (med/lat gas geometry)
        cc = 'ankle_angle_r';
    otherwise
        if emptyMiddle && strcmp(body2, 'calcn_r') && ~strcmp(body1, 'talus_r')
            cc = 'ankle_angle_r';
        else
            cc = primary;
        end
end
end

function cc = chainCoordinateFor2(body2, body3, primary, coord1)
key = sprintf('%s>%s', body2, body3);
switch key
    case 'femur_r>tibia_r'
        cc = 'knee_angle_r';
    case 'talus_r>calcn_r'
        cc = 'subtalar_angle_r';
    case 'calcn_r>toes_r'
        cc = 'mtp_angle_r';
    case 'tibia_r>calcn_r'
        cc = 'ankle_angle_r';
    otherwise
        cc = primary;
end
if strcmp(cc, coord1)
    cc = primary;   % never duplicate the same coordinate on both crossings
end
end

%% ---------------------------------------------------------------------
function [points, cross, frames, emptyMiddle, rigidSlot2, notes] = ...
    groupSlots(name, pts, depths)
% Frame-slot assignment (BiPamData 3-slot convention). Depths are the
% chain depths pelvis 0 < femur 1 < tibia 2 < talus 3 < calcn 4 < toes 5,
% torso 6 (direct child of pelvis, see the trunk special case).
n = numel(depths);
ud = unique(depths);
if numel(ud) < 2
    error('routingSpecsFromRobotbody:singleFrame', ...
        '%s: path crosses no joint', name)
end

chainBodies = {'pelvis', 'femur_r', 'tibia_r', 'talus_r', 'calcn_r', 'toes_r'};
notes = {};
emptyMiddle = false;
rigidSlot2 = false;
if numel(ud) >= 3
    slotDepths = ud(1:3);
elseif ud(1) == 0 && ud(2) == 6
    % Trunk: pelvis -> torso is a DIRECT chain edge (back joint).
    slotDepths = [0, 6, 6];
    rigidSlot2 = true;
elseif ud(2) - ud(1) > 1
    % Skip-frame span: insert the next chain body as an empty middle so
    % BOTH joints stay expressible (hamstrings femur; gas tibia via the
    % ankle alias; tib->calcn the true talus/subtalar).
    slotDepths = [ud(1), ud(1) + 1, ud(2)];
    emptyMiddle = true;
else
    slotDepths = [ud(1), ud(2), ud(2)];
    rigidSlot2 = true;
end

slot = zeros(n, 1);
points = zeros(n, 3);
for i = 1:n
    s = find(depths(i) == slotDepths, 1);
    if isempty(s)
        error('routingSpecsFromRobotbody:depthJump', ...
            '%s: point depth %d outside slots %s', name, depths(i), ...
            mat2str(slotDepths))
    end
    slot(i) = s;
    points(i, :) = pts{i, 2};
end

frames = cell(1, 3);
for s = 1:3
    hit = find(slot == s, 1);
    if ~isempty(hit)
        frames{s} = pts{hit, 1};
    elseif slotDepths(s) <= 5
        frames{s} = chainBodies{slotDepths(s) + 1};
        frames{s} = [frames{s}, ' (no path points)'];
    else
        frames{s} = 'torso (no path points)';
    end
end

cross = zeros(1, 2);
cross(1) = firstRowOrZero(slot == 2);
cross(2) = firstRowOrZero(slot == 3);
if emptyMiddle
    cross(1) = cross(2);    % duplicate-crossing convention
    notes{end + 1} = 'empty middle frame: duplicate crossing, transforms chained';
elseif rigidSlot2
    if cross(1) < n
        cross(2) = cross(1) + 1;
        frames{3} = [frames{2}, ' (rigid)'];
        notes{end + 1} = 'slot 3 = rigid clone; crossing-2 outputs inert';
    else
        points(end + 1, :) = points(end, :);
        cross(2) = n + 1;
        frames{3} = [frames{2}, ' (rigid, duplicate row)'];
        notes{end + 1} = 'duplicate insertion row appended; crossing-2 outputs inert';
    end
end
end

function idx = firstRowOrZero(mask)
hit = find(mask, 1);
if isempty(hit)
    idx = 0;
else
    idx = hit;
end
end

%% ---------------------------------------------------------------------
function R = axisRot(axis, theta)
a = axis(:) / norm(axis);
Km = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * Km + (1 - cos(theta)) * (Km * Km);
end

%% ---------------------------------------------------------------------
function [joints, note] = applyJoint(joints, slot, bodyParent, bodyChild, ...
    coordName, kin, coords, mapStruct, varargin) %#ok<INUSD>
key = sprintf('%s_to_%s', bodyParent, bodyChild);
if strcmp(key, 'tibia_r_to_calcn_r') && ~isKey(kin.joints, key)
    key = 'tibia_r_to_talus_r';    % true-chain lookup
end
if isKey(kin.joints, key) && isKey(kin.joints(key).axes, coordName)
    J = kin.joints(key);
    A = J.axes(coordName);
    joints.types{slot} = 'hinge';
    joints.pivots(slot, :) = J.pivot;
    joints.axes(slot, :) = A.axis;
    if isfield(A, 'transThetaX') && ~isempty(A.transThetaX)
        joints.transThetaX = A.transThetaX;
        joints.transX = A.transX;
    end
    if isfield(A, 'transThetaY') && ~isempty(A.transThetaY)
        joints.transThetaY = A.transThetaY;
        joints.transY = A.transY;
    end
    joints.thetaRangesDeg(slot, :) = campaignRange(coordName, coords, ...
        mapStruct);
    if strcmp(coordName, 'knee_angle_r')
        joints.kneeSlot = slot;   % the rolling-knee lives at this crossing
    end
    note = sprintf('crossing %d = %s hinge on %s (axis [%g %g %g], pivot [%g %g %g])', ...
        slot, J.name, coordName, A.axis, J.pivot);
else
    joints.types{slot} = 'fixed';
    note = sprintf('crossing %d: no %s hinge for %s, held rigid', slot, ...
        coordName, key);
end
end

%% ---------------------------------------------------------------------
function r = campaignRange(coordName, coords, mapStruct)
% Campaign RoM: the map's rom_table_deg first (it carries the documented
% judgment values for mtp/lumbar, whose model ranges are EMPTY), then the
% committed biPulleySpecsFromOpenSim convention, then the model range.
if isfield(mapStruct, 'rom_table_deg') && ...
        isfield(mapStruct.rom_table_deg, coordName) && ...
        ~isempty(mapStruct.rom_table_deg.(coordName))
    r = double(mapStruct.rom_table_deg.(coordName));
    if numel(r) == 2
        return
    end
end
switch coordName
    case 'hip_flexion_r'
        r = [-15, 110];
    case 'knee_angle_r'
        r = [-120, 10];
    case 'ankle_angle_r'
        r = [-50, 20];
    case 'subtalar_angle_r'
        r = [-25, 35];
    otherwise
        if isKey(coords, coordName)
            r = coords(coordName).range_deg;
        else
            r = [0, 0];
        end
end
end

%% ---------------------------------------------------------------------
function kin = osimKinematics(model)
% Parse JointSet CustomJoints: per joint, the pivot (parent-offset
% translation) and a per-coordinate map of TransformAxis data (axis +
% SimmSpline translations bound to the same coordinate).
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

    axes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    taNodes = j.getElementsByTagName('TransformAxis');
    for t = 0:taNodes.getLength() - 1
        ta = taNodes.item(t);
        cn = childText(ta, 'coordinates');
        if isempty(cn)
            continue
        end
        coordParts = regexp(cn, ',', 'split');
        coordName = strtrim(coordParts{1});
        ax = sscanf(childText(ta, 'axis'), '%f')';
        A = struct('axis', ax, 'transTheta', [], 'transX', [], 'transY', []);
        % Translation splines bound to the same coordinate (rolling knee).
        for t2 = 0:taNodes.getLength() - 1
            ta2 = taNodes.item(t2);
            cn2 = childText(ta2, 'coordinates');
            if isempty(cn2)
                continue
            end
            parts = regexp(regexprep(cn2, ',', ' '), '\s+', 'split');
            if ~any(strcmp(parts, coordName))
                continue
            end
            ax2 = sscanf(childText(ta2, 'axis'), '%f')';
            fn = ta2.getElementsByTagName('SimmSpline');
            if fn.getLength() == 0
                continue
            end
            fnode = fn.item(0);
            xs = sscanf(childText(fnode, 'x'), '%f')';
            ys = sscanf(childText(fnode, 'y'), '%f')';
            if numel(xs) ~= numel(ys) || numel(xs) < 2
                continue
            end
            if all(abs(ax2 - [1, 0, 0]) < 1e-9)
                A.transThetaX = xs;    % x translation spline knots
                A.transX = ys;
            elseif all(abs(ax2 - [0, 1, 0]) < 1e-9)
                A.transThetaY = xs;    % y translation spline knots
                A.transY = ys;
            end
        end
        if ~isKey(axes, coordName)
            axes(coordName) = A;
        end
    end

    parentBody = regexprep(parentOff.frame, '_offset$', '');
    childBody = regexprep(childOff.frame, '_offset$', '');
    J = struct('pivot', parentOff.t, 'name', jname, 'axes', axes);
    kin.joints(sprintf('%s_to_%s', parentBody, childBody)) = J;
    if strcmp(jname, 'ankle_r')
        % Committed ankle alias (subtalar frozen) for tibia->calcn pairs.
        kin.joints('tibia_r_to_calcn_r') = J;
    end
end
end

%% ---------------------------------------------------------------------
function off = offsetFrame(jointNode, socketTag)
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
