function kin = fr_kinematics(osimPath)
%FR_KINEMATICS Parse a gait2392-family .osim into a continuous kinematics model.
%
% Reads CustomJoint definitions (PhysicalOffsetFrames, SpatialTransform axes
% with LinearFunction / Constant / SimmSpline functions), Coordinate
% defaults and ranges, and body names. Transforms are composed downstream
% by fr_bodyT / fr_chainT.
%
% COMPOSITION CONVENTION (validated in this session against
% geometry_20261103/extraction_robotbody.json bodyTransformsInGround at
% knee 0 / -60 / -120 deg - exact at the spline-node knee angles; see
% fr_validateKin): X_joint(q) = [R(q), t(q)] mapping child-offset coords
% to parent-offset coords, with
%   t(q) = sum_i a_trans_i * f_i(q)  (parent-offset axes, placed FIRST)
%   R(q) = Rot(a1,q1)*Rot(a2,q2)*Rot(a3,q3) (axes through the parent-offset
%          origin, i.e. about the translated child origin),
% then X_parentBody<-childBody = X_body<-parentOffset * X_joint *
% inv(X_body<-childOffset). Because the offset frames sit AT the joints
% (hip pelvis_offset = hip pivot; femur/tibia body origins = hip/knee
% pivots; talus origin = ankle pivot; subtalar calcn_r_offset pivot; MTP
% calcn_r_offset pivot = toes origin; back pelvis_offset = lumbar pivot),
% every joint transform's translation IS the joint pivot, which is exactly
% the transform requirement of BiPam_pulley2_20261022.

kin = struct();
kin.osimPath = osimPath;
doc = xmlread(osimPath);
root = doc.getDocumentElement();

% ---- bodies -------------------------------------------------------------
bodyNodes = root.getElementsByTagName('Body');
kin.bodies = cell(1, bodyNodes.getLength());
for i = 1:bodyNodes.getLength()
    kin.bodies{i} = char(bodyNodes.item(i-1).getAttribute('name'));
end

% ---- coordinates (defaults and ranges, radians / meters) ----------------
coordNodes = root.getElementsByTagName('Coordinate');
kin.coordDefault = containers.Map();
kin.coordRange = containers.Map();
for i = 0:coordNodes.getLength()-1
    cn = coordNodes.item(i);
    nm = char(cn.getAttribute('name'));
    dv = fr_xmlChild(cn, 'default_value');
    if ~isempty(dv)
        kin.coordDefault(nm) = fr_xmlNum(dv);
    end
    rg = fr_xmlChild(cn, 'range');
    if ~isempty(rg)
        kin.coordRange(nm) = fr_xmlNums(rg);
    end
end

% ---- offset frames: PER-JOINT scope (the names collide across joints:
% ground_pelvis / hip_r / back EACH define their own 'pelvis_offset' with
% DIFFERENT translations, and knee_r vs ankle_r define different
% 'tibia_r_offset' values - the frame name is only unique within its
% owning joint) ------------------------------------------------------------
jointNodes = root.getElementsByTagName('CustomJoint');

% ---- joints -------------------------------------------------------------
kin.joints = struct();
for i = 0:jointNodes.getLength()-1
    jn = jointNodes.item(i);
    J = struct();
    J.name = char(jn.getAttribute('name'));
    J.parentFrame = fr_stripFramePath(fr_xmlChild(jn, 'socket_parent_frame'));
    J.childFrame = fr_stripFramePath(fr_xmlChild(jn, 'socket_child_frame'));
    J.frames = struct();
    offNodes = jn.getElementsByTagName('PhysicalOffsetFrame');
    for k = 0:offNodes.getLength()-1
        fn = offNodes.item(k);
        F = struct();
        F.name = char(fn.getAttribute('name'));
        F.body = fr_stripFramePath(fr_xmlChild(fn, 'socket_parent'));
        F.t = fr_xmlNums(fr_xmlChild(fn, 'translation'));
        F.R = fr_xyzExtrinsic(fr_xmlNums(fr_xmlChild(fn, 'orientation')));
        J.frames.(F.name) = F;
    end
    st = jn.getElementsByTagName('SpatialTransform');
    taNodes = st.item(0).getElementsByTagName('TransformAxis');
    J.axes = cell(1, 6);
    for k = 0:taNodes.getLength()-1
        an = taNodes.item(k);
        A = struct();
        A.axis = fr_xmlNums(fr_xmlChild(an, 'axis'));
        co = fr_xmlChild(an, 'coordinates');
        if isempty(co) || isempty(strtrim(char(co.getTextContent())))
            A.coords = {};
        else
            A.coords = strsplit(strtrim(char(co.getTextContent())));
        end
        A.location = [0 0 0];   % gait2392 rotation axes sit at the offset origin
        fnl = fr_xmlChild(an, 'LinearFunction');
        fcn = fr_xmlChild(an, 'Constant');
        fsp = fr_xmlChild(an, 'SimmSpline');
        if ~isempty(fnl)
            A.type = 'linear';
            A.coeffs = fr_xmlNums(fnl);   % [scale offset]
        elseif ~isempty(fcn)
            A.type = 'constant';
            A.value = fr_xmlNum(fcn);
        elseif ~isempty(fsp)
            A.type = 'spline';
            A.x = fr_xmlNums(fr_xmlChild(fsp, 'x'));
            A.y = fr_xmlNums(fr_xmlChild(fsp, 'y'));
        else
            error('fr_kinematics:AxisFunction', ...
                'Unsupported TransformAxis function in joint %s', J.name);
        end
        J.axes{k+1} = A;
    end
    kin.joints.(J.name) = J;
end

% ---- joint chain (child body -> parent body + joint) ---------------------
kin.chain = struct();
fn = fieldnames(kin.joints);
for i = 1:numel(fn)
    J = kin.joints.(fn{i});
    S = struct('parentBody', fr_frameOwner(kin, J.parentFrame), ...
        'jointName', J.name);
    kin.chain.(fr_frameOwner(kin, J.childFrame)) = S;
end

% ---- coordinate order used by the campaign ------------------------------
kin.coordNames = {'pelvis_tilt', 'pelvis_list', 'pelvis_rotation', ...
    'pelvis_tx', 'pelvis_ty', 'pelvis_tz', ...
    'lumbar_extension', 'lumbar_bending', 'lumbar_rotation', ...
    'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r', ...
    'knee_angle_r', 'ankle_angle_r', 'subtalar_angle_r', 'mtp_angle_r'};

kin.qDefault = struct();
for i = 1:numel(kin.coordNames)
    nm = kin.coordNames{i};
    if isKey(kin.coordDefault, nm)
        kin.qDefault.(nm) = kin.coordDefault(nm);
    else
        kin.qDefault.(nm) = 0;
    end
end
end

function b = fr_frameOwner(kin, frameName)
% frame owner across ALL joints' scoped offset frames (the owning BODY is
% the same wherever the frame name appears; only the transform differs)
fn = fieldnames(kin.joints);
for i = 1:numel(fn)
    if isfield(kin.joints.(fn{i}).frames, frameName)
        b = kin.joints.(fn{i}).frames.(frameName).body;
        return
    end
end
b = frameName;   % not an offset frame: it is a body name
end

function p = fr_stripFramePath(s)
p = strrep(strrep(strtrim(char(s.getTextContent())), '/bodyset/', ''), ...
    '/jointset/', '');
p = strrep(p, '/ground', 'ground');  % the ground body lives at /ground
end
