% Bi Pulley Specs From Masters 20261023
% Author: Ben Bolen (base builder, 2026-09-24); dated extension 2026-10-23 by
% the machinery-extension engineer per FIX_SPEC_20261023.md section 3.
%
% NEW DATED COPY of Code/Matlab/Mesh_Optimization/biPulleySpecsFromMasters.m
% (copy-and-mark-edits convention; Ben's original untouched). Everything the
% original does is KEPT VERBATIM (the 4 established routes, the patella rule,
% the field layout); this copy extends the route table to ALL 27 campaign
% actuators and adds the three hinge families the FIX_SPEC names:
%   (a) the BACK joint pelvis<->torso, center [-0.1007, 0.0815, 0] m from the
%       pelvis (Ben's value; it equals the stock gait2392_simbody.osim back
%       joint's pelvis_offset translation), lumbar AXES from the stock model
%       (lumbar_extension about [0 0 1], lumbar_bending about [1 0 0],
%       lumbar_rotation about [0 1 0] -- the back joint's TransformAxis
%       blocks), so ercspn/intobl/extobl build and grade on lumbar DOFs;
%   (b) SUBTALAR (tibia<->calcn) as a REAL hinge -- the masters' own serial
%       chain VecTrans(T_a*T_s, x) with the ankle transform T_a = [R_a, p_a]
%       and the subtalar transform T_s = [R_s, p_s] -- for tib_post and the
%       peroneals (and every other tibia->calcn route);
%   (c) MTP where the toes are the primary DOF, following the masters' own
%       toe convention (bipedal_kinematics_toe.m sweeps ankle x mtp with the
%       subtalar pinned at 0).
%
% DATA OF RECORD (FIX_SPEC section 3): the masters attachment data --
%   AttachPoints_Robot.mat, PAM_Prop.mat, MovePoints.mat -- read from the
%   campaign's read-only repo copies at tmp/old_code_copy_20261022/ by
%   default (the thumb-drive original 'D:\Bipedal humanoid\Matlab code' is
%   the fallback). This is the SAME geometry the hand-tuned row uses.
%
% KEPT VERBATIM FROM THE ORIGINAL (biPulleySpecsFromMasters.m):
%   - the 4 routes semimem_r / rect_fem_r / med_gas_r / tfl_r: identical
%     point lists, cross indices, hinge pairs (hip-flexion+knee for the
%     pelvis routes, knee+ankle for med_gas), pivots, axes, ranges;
%   - the PATELLA RULE (rect_fem's patella function point pa1 replaced by
%     the fixed tibia point t3, rigid patella, "connects patella") -- and the
%     same rule applied to vas_int's patella function point pa2 (the
%     AttachPoints_Robot.m Patella section comments "pa1 %insertion:
%     rect_fem" / "pa2 %insertion: vas_int" name the same MovePoints
%     construct for both);
%   - the twr1-as-tibia ruling, the med_gas duplicate-crossing treatment
%     (calcn frame ankle-centered: the subtalar origin offset p_s and its
%     rotation are NOT modeled on the knee+ankle routes), the transmission
%     config defaults, and every spec field name.
%
% NEW HINGE / ROUTE DECISIONS (each grounded in the masters' own files):
%   JOINT CATALOG (bipedal_kinematics_back.m lines 83-115, the file that
%   defines every transform; ranges from the campaign rom_table_deg in
%   actuator_map_27.json, which matches the masters' own RoM block):
%     hip_flexion  axis [0 0 1]   pivot p_h=[-0.0707 -0.0661 0.0835] pelvis
%     hip_adduction axis [1 0 0]  pivot p_h (masters R_h_x; rom box [-45 20],
%                  the coupled flexion envelope governs grading, FIX_SPEC 1b)
%     knee         axis [0 0 1]   rolling pivot [fcn1 fcn2 0] femur (MovePoints)
%     ankle        axis omg_a=[-0.10501355 -0.17402245 0.97912632] pivot
%                  p_a=[0 -0.43 0] tibia
%     subtalar     axis omg_s=[0.7871796 0.60474746 -0.12094949] pivot
%                  p_s=[-0.04877 -0.04195 0.00792] in the ankle-child frame
%     mtp          axis omg_m=[-0.5809544 0 0.81393611] pivot
%                  p_m=[0.1788 -0.002 0.00108] in the subtalar-child frame
%     lumbar_extension axis [0 0 1] (stock) pivot p_b=[-0.1007 0.0815 0] pelvis
%     lumbar_bending   axis [1 0 0] (stock) pivot p_b
%   SLOT-2 PIVOT RULE: the class composes frame1 <-T1- frame2 <-T2- frame3,
%   and each hinge's pivot lives in ITS chain parent frame. When the joint's
%   anatomical parent is the chain frame, the catalog pivot is used as-is.
%   Where the chain's middle frame is a JOINT-CENTERED COPY of the anatomical
%   parent (lumbar trio, hip-uniarticular group) or a SHIFTED copy (toe
%   routes' ankle-centered calcn; psoas' back-centered pelvis), dedicated
%   catalog entries carry the re-expressed pivot:
%     lumbar_extension_chain2  pivot [0 0 0]   (torsoB origin = back center)
%     hip_flexion_femurA       pivot [0 0 0]   (femurB origin = hip center)
%     mtp_calcnA               pivot p_m + p_s (mtp center in calcnA coords)
%     lumbar_extension_torsoView axis [0 0 -1] pivot [0 0 0] (psoas slot 1)
%     hip_flexion_pelvisShift  pivot p_h - p_b (psoas slot 2)
%   ROUTE GROUPS:
%   - LUMBAR TRIO ercspn/intobl/extobl: frames pelvis -> torsoA -> torsoB,
%     hinge 1 = lumbar_bending, hinge 2 = lumbar_extension (the masters'
%     compound R_b = R_b_x*R_b_y*R_b_z makes extension the innermost
%     rotation, so slot 2 = extension; lumbar_rotation pinned at 0, stated
%     in caveats). Primary lumbar_extension grades at slot 2.
%   - HAMSTRING/KNEE GROUP (bifemlh, sar, grac join the original semimem/
%     tfl/rect_fem): frames pelvis -> femur -> tibia, hip flexion + knee.
%     bifemlh has no femur rows, so it uses the duplicate-crossing chain
%     x_pelvis = T_h*(T_k*x) -- exactly the masters' VecTrans(T_h*T_k, x).
%   - KNEE-ONLY bifemsh/vas_int: frames femur -> tibia -> tibia(rigid). The
%     rigid-clone construction follows the OpenSim builder's documented
%     convention (biPulleySpecsFromOpenSim.m header note 5): the insertion
%     row is duplicated so the two-crossing bookkeeping stays well formed
%     and crossing 2's outputs are inert.
%   - GASTROCS: lat_gas joins med_gas verbatim-style: frames femur ->
%     [empty tibia] -> calcn, knee + ankle, duplicate crossing, calcn frame
%     ankle-centered (the established original med_gas convention; the
%     masters' subtalar layer dropped on these two routes only -- the class
%     models two hinges, knee+ankle here, subtalar pinned).
%   - TIBIA->CALCN GROUP soleus/tib_post/tib_ant/per_brev/per_long/per_tert:
%     frames tibia -> [empty] -> calcn(subtalar child), duplicate crossing,
%     hinge 1 = ankle [R_a p_a], hinge 2 = subtalar [R_s p_s]. This IS the
%     masters' VecTrans(T_a*T_s, x) exactly, subtalar ANIMATED. tib_post and
%     the peroneals grade subtalar at slot 2 (the FIX_SPEC fix); soleus and
%     tib_ant grade ankle at slot 1. NOTE soleus origin: the masters'
%     commented soleus formula uses t1, but the attachment file's own
%     comment on t4 reads "origin: soleus, tib_post, flex_dig, flex_hal" --
%     the attachment file (geometry of record) is followed: t4 -> a1.
%   - TOE GROUP flex_dig/flex_hal/ext_dig/ext_hal: frames tibia -> calcnA ->
%     toes, hinge 1 = ankle, hinge 2 = mtp, subtalar PINNED at 0 (the
%     masters' own toe file sweeps ankle x mtp with theta_s = 0). The calcn
%     rows are PRE-FOLDED by +p_s into ankle-centered (calcnA) coordinates
%     and the mtp pivot is carried as p_m + p_s, so at theta_s = 0 the
%     composed chain reproduces the masters' VecTrans(T_a*T_s*T_m, x)
%     EXACTLY (both the p_s translation and the point placement). The folds
%     are applied HERE (recorded in spec.kin.rowFolds); the driver-side
%     pointFolds stay zero by contract.
%   - HIP UNIARTICULAR glut_max/add_mag/iliacus: frames pelvis -> femurA ->
%     femurB, hinge 1 = hip_adduction [R_x p_h], hinge 2 = hip_flexion
%     [R_z 0] (the masters' compound R_h = R_h_x*R_h_y*R_h_z makes flexion
%     the innermost rotation; hip_rotation pinned at 0). add_mag grades
%     adduction at slot 1; glut_max/iliacus grade flexion at slot 2. At
%     adduction = 0 the chain reduces exactly to the original builder's
%     sagittal hip convention.
%   - PSOAS: the masters' own route starts on the TORSO (b1 -> psoas_wr1 ->
%     f4, uni_hip file line 176). The class chain is torso(proximal) ->
%     shifted-pelvis(middle) -> femur(distal): the pelvis row is pre-shifted
%     by -p_b so the middle frame's origin IS the back-joint center, hinge
%     1 = lumbar_extension as a PURE rotation [R_z(theta)' 0] (axis carried
%     as [0 0 -1] so the standard positive grid = positive OpenSim
%     lumbar_extension), hinge 2 = hip_flexion with pivot p_h - p_b. At the
%     neutral pose the route lengths equal the masters' formula exactly;
%     both hinge torques are about the correct pivots. CAVEAT: slot-1's
%     lumbar torque is measured in the pelvis-relative-to-torso view
%     (negate to compare with pelvis-referenced lumbar torques); psoas
%     grades hip_flexion at slot 2, which is unaffected.
%
% Usage:
%   specs = biPulleySpecsFromMasters_20261023();        % all 27 (map order)
%   specs = biPulleySpecsFromMasters_20261023({'semimem_r', 'tib_post_r'});
%   specs = biPulleySpecsFromMasters_20261023([], mastersDir);
%   specs = biPulleySpecsFromMasters_20261023([], [], struct('nPulleyBPA', 2));
% Spec fields match the original plus: spec.kin.dofs {1x2 coordinate names},
% spec.kin.jointNames, spec.kin.primarySlot, spec.kin.primaryDof,
% spec.kin.rowFolds (the exact per-row offsets this builder applied),
% spec.actuatorId / spec.masterName (actuator-map ids 1..27 by route order).

function specs = biPulleySpecsFromMasters_20261023(names, mastersDir, pulleyConfig)

%% Inputs and data location
[allNames, allIds] = masterRouteNames_20261023();
if nargin < 1 || isempty(names)
    names = allNames;
end
if nargin < 2 || isempty(mastersDir)
    % Campaign copies of record first (read-only); thumb drive fallback.
    here = fileparts(mfilename('fullpath'));
    repoRoot = iRepoRoot(here);
    repoCopy = fullfile(repoRoot, 'tmp', 'old_code_copy_20261022');
    if exist(fullfile(repoCopy, 'AttachPoints_Robot.mat'), 'file')
        mastersDir = repoCopy;
    elseif exist('D:\Bipedal humanoid\Matlab code\AttachPoints_Robot.mat', 'file')
        mastersDir = 'D:\Bipedal humanoid\Matlab code';
    else
        mastersDir = repoCopy;   % let the exist() checks below report it
    end
end
if nargin < 3
    pulleyConfig = struct();
end

attachFile = fullfile(mastersDir, 'AttachPoints_Robot.mat');
pamFile = fullfile(mastersDir, 'PAM_Prop.mat');
moveFile = fullfile(mastersDir, 'MovePoints.mat');
if ~exist(attachFile, 'file')
    error('biPulleySpecsFromMasters_20261023:missingData', ...
        'Master''s attachment data not found at %s', attachFile)
end
if ~exist(pamFile, 'file')
    error('biPulleySpecsFromMasters_20261023:missingData', ...
        'Master''s PAM properties not found at %s', pamFile)
end

S = load(attachFile);
P = load(pamFile);

%% Joint catalog (sources: bipedal_kinematics_back.m lines 83-115;
%% lumbar axes from the stock gait2392_simbody.osim back joint; ranges from
%% the campaign rom_table_deg, actuator_map_27.json). dof = campaign
%% coordinate name carried in spec.kin.dofs.
pHip = [-0.0707, -0.0661, 0.0835];
pAnkle = [0, -0.43, 0];
omgAnkle = [-0.10501355, -0.17402245, 0.97912632];
pSub = [-0.04877, -0.04195, 0.00792];
omgSub = [0.7871796, 0.60474746, -0.12094949];
pMtp = [0.1788, -0.002, 0.00108];
omgMtp = [-0.5809544, 0, 0.81393611];
pBack = [-0.1007, 0.0815, 0];          % Ben's value; = stock pelvis_offset

J = struct();
J.hip_flexion    = struct('axis', [0, 0, 1],  'pivot', pHip, ...
    'range', [-15, 110], 'dof', 'hip_flexion_r', 'rolling', false);
J.hip_adduction  = struct('axis', [1, 0, 0],  'pivot', pHip, ...
    'range', [-45, 20], 'dof', 'hip_adduction_r', 'rolling', false);
J.knee           = struct('axis', [0, 0, 1],  'pivot', [], ...
    'range', [-120, 10], 'dof', 'knee_angle_r', 'rolling', true);
J.ankle          = struct('axis', omgAnkle,   'pivot', pAnkle, ...
    'range', [-50, 20], 'dof', 'ankle_angle_r', 'rolling', false);
J.subtalar       = struct('axis', omgSub,     'pivot', pSub, ...
    'range', [-25, 35], 'dof', 'subtalar_angle_r', 'rolling', false);
J.mtp            = struct('axis', omgMtp,     'pivot', pMtp, ...
    'range', [-40, 40], 'dof', 'mtp_angle_r', 'rolling', false);
J.lumbar_extension = struct('axis', [0, 0, 1], 'pivot', pBack, ...
    'range', [-20, 15], 'dof', 'lumbar_extension', 'rolling', false);
J.lumbar_bending   = struct('axis', [1, 0, 0], 'pivot', pBack, ...
    'range', [-15, 15], 'dof', 'lumbar_bending', 'rolling', false);
% Slot-2 / shifted-view variants (see header SLOT-2 PIVOT RULE).
J.lumbar_extension_chain2 = struct('axis', [0, 0, 1], 'pivot', [0, 0, 0], ...
    'range', [-20, 15], 'dof', 'lumbar_extension', 'rolling', false);
J.hip_flexion_femurA = struct('axis', [0, 0, 1], 'pivot', [0, 0, 0], ...
    'range', [-15, 110], 'dof', 'hip_flexion_r', 'rolling', false);
J.mtp_calcnA = struct('axis', omgMtp, 'pivot', pMtp + pSub, ...
    'range', [-40, 40], 'dof', 'mtp_angle_r', 'rolling', false);
J.lumbar_extension_torsoView = struct('axis', [0, 0, -1], ...
    'pivot', [0, 0, 0], 'range', [-20, 15], 'dof', 'lumbar_extension', ...
    'rolling', false);
J.hip_flexion_pelvisShift = struct('axis', [0, 0, 1], ...
    'pivot', pHip - pBack, 'range', [-15, 110], 'dof', 'hip_flexion_r', ...
    'rolling', false);

thetaHipKnee = [-15, 110; -120, 10];     % kept name (original builder)
thetaKneeAnkle = [-120, 10; -50, 20];

%% Rolling-knee grids (knee center in femur coords vs knee angle, rad).
kneeFit = struct('thetaX', [], 'x', [], 'thetaY', [], 'y', []);
if exist(moveFile, 'file')
    M = load(moveFile, 'knee_angle_x', 'knee_x', 'knee_angle_y', 'knee_y');
    kneeFit.thetaX = M.knee_angle_x(:).';
    kneeFit.x = M.knee_x(:).';
    kneeFit.thetaY = M.knee_angle_y(:).';
    kneeFit.y = M.knee_y(:).';
else
    warning('biPulleySpecsFromMasters_20261023:noMovePoints', ...
        'MovePoints.mat not found; knee pivot held fixed at the neutral sample.')
    kneeFit.thetaX = 0;
    kneeFit.x = 0.047;
    kneeFit.thetaY = 0;
    kneeFit.y = -0.113;
end
kneeNeutral = [interp1(kneeFit.thetaX, kneeFit.x, 0), ...
    interp1(kneeFit.thetaY, kneeFit.y, 0), 0];
J.knee.pivot = kneeNeutral;

%% Route table -- ALL 27 campaign actuators (map order,
%% Routing27_v2/actuator_map_27.json). Columns:
%   {specName, mifField, frames, pointNames, nProx, nMid, hinge1, hinge2,
%    primaryDof, primarySlot, foldTag, note}
routes = { ...
    'ercspn_r',   'ercspn_mif',   {'pelvis', 'torsoA', 'torsoB'}, ...
        {'p7', 'ercspn_wr1', 'ercspn_wr2', 'ercspn_wr3', 'ercspn_wr4', ...
         'ercspn_wr5', 'ercspn_wr6', 'ercspn_wr7', 'b2'}, 2, 0, ...
        'lumbar_bending', 'lumbar_extension_chain2', 'lumbar_extension', ...
        2, 'none', ...
        'long route (back file line 266, the ACTIVE formula)'; ...
    'intobl_r',   'intobl_mif',   {'pelvis', 'torsoA', 'torsoB'}, ...
        {'p2', 'b3'}, 1, 0, ...
        'lumbar_bending', 'lumbar_extension_chain2', 'lumbar_extension', ...
        2, 'none', 'back file line 187 (p2 -> b3)'; ...
    'extobl_r',   'extobl_mif',   {'pelvis', 'torsoA', 'torsoB'}, ...
        {'p4', 'b4'}, 1, 0, ...
        'lumbar_bending', 'lumbar_extension_chain2', 'lumbar_extension', ...
        2, 'none', 'back file line 193 (p4 -> b4)'; ...
    'bifemlh_r',  'bifemlh_mif',  {'pelvis', 'femur', 'tibia'}, ...
        {'p5', 'twr2', 'twr3', 't5'}, 1, 0, ...
        'hip_flexion', 'knee', 'knee_angle_r', 2, 'none', ...
        'bi_hip line 132 (T_h*T_k chain = duplicate crossing)'; ...
    'semimem_r',  'semimem_mif',  {'pelvis', 'femur', 'tibia'}, ...
        {'p5', 'semimem_wr1', 'semimem_wr2', 'twr1', 't1', 't6'}, 1, 2, ...
        'hip_flexion', 'knee', 'knee_angle_r', 2, 'none', ...
        'ORIGINAL ROUTE VERBATIM'; ...
    'sar_r',      'sar_mif',      {'pelvis', 'femur', 'tibia'}, ...
        {'p3', 'sar_wr1', 'sar_wr2', 'sar_wr3', 'twr1', 't1', 't6'}, 1, 3, ...
        'hip_flexion', 'knee', 'knee_angle_r', 2, 'none', ...
        'bi_hip line 143'; ...
    'tfl_r',      'tfl_mif',      {'pelvis', 'femur', 'tibia'}, ...
        {'p3', 'tfl_wr1', 'tfl_wr2', 't2', 't5'}, 1, 2, ...
        'hip_flexion', 'knee', 'knee_angle_r', 2, 'none', ...
        'ORIGINAL ROUTE VERBATIM'; ...
    'rect_fem_r', 'rect_fem_mif', {'pelvis', 'femur', 'tibia'}, ...
        {'p10', 'rect_fem_wr1', 'rect_fem_wr2', 'rect_fem_wr3', 't3'}, 1, 3, ...
        'hip_flexion', 'knee', 'knee_angle_r', 2, 'none', ...
        'ORIGINAL ROUTE VERBATIM (patella rule: pa1 -> t3)'; ...
    'grac_r',     'grac_mif',     {'pelvis', 'femur', 'tibia'}, ...
        {'p8', 'grac_wr1', 'grac_wr2', 'twr1', 't1', 't6'}, 1, 2, ...
        'hip_flexion', 'knee', 'knee_angle_r', 2, 'none', ...
        'bi_hip line 163'; ...
    'bifemsh_r',  'bifemsh_mif',  {'femur', 'tibia', 'tibia_rigid'}, ...
        {'f2', 'twr2', 'twr3'}, 1, 1, ...
        'knee', 'rigid', 'knee_angle_r', 1, 'none', ...
        'bi_hip line 139 (ACTIVE); rigid clone for the 2nd crossing'; ...
    'vas_int_r',  'vas_int_mif',  {'femur', 'tibia', 'tibia_rigid'}, ...
        {'f1', 'vas_int_wr1', 'vas_int_wr2', 't3'}, 1, 2, ...
        'knee', 'rigid', 'knee_angle_r', 1, 'none', ...
        'bi_hip line 206 (ACTIVE); patella rule: pa2 -> t3; rigid clone'; ...
    'soleus_r',   'soleus_mif',   {'tibia', 'calcnA', 'calcnSub'}, ...
        {'t4', 'a1'}, 1, 0, ...
        'ankle', 'subtalar', 'ankle_angle_r', 1, 'none', ...
        't4 per AttachPoints comment (masters formula says t1; file followed)'; ...
    'med_gas_r',  'med_gas_mif',  {'femur', 'tibia', 'calcn'}, ...
        {'f4', 'med_gas_wr1', 'a1'}, 2, 0, ...
        'knee', 'ankle', 'ankle_angle_r', 2, 'none', ...
        'ORIGINAL ROUTE VERBATIM (calcn ankle-centered; subtalar dropped)'; ...
    'lat_gas_r',  'lat_gas_mif',  {'femur', 'tibia', 'calcn'}, ...
        {'f5', 'lat_gas_wr1', 'a1'}, 2, 0, ...
        'knee', 'ankle', 'ankle_angle_r', 2, 'none', ...
        'calves line 204; same convention as med_gas (subtalar dropped)'; ...
    'tib_post_r', 'tib_post_mif', {'tibia', 'calcnA', 'calcnSub'}, ...
        {'t4', 'tib_post_wr1', 'tib_post_wr2', 'a4'}, 2, 0, ...
        'ankle', 'subtalar', 'subtalar_angle_r', 2, 'none', ...
        'foot line 215 (T_a*T_s chain: subtalar REAL, slot 2)'; ...
    'tib_ant_r',  'tib_ant_mif',  {'tibia', 'calcnA', 'calcnSub'}, ...
        {'t2', 'tib_ant_wr1', 'a3'}, 2, 0, ...
        'ankle', 'subtalar', 'ankle_angle_r', 1, 'none', ...
        'foot line 232'; ...
    'per_brev_r', 'per_brev_mif', {'tibia', 'calcnA', 'calcnSub'}, ...
        {'t2', 'per_brev_wr1', 'per_brev_wr2', 'per_brev_wr3', 'a2'}, 3, 0, ...
        'ankle', 'subtalar', 'subtalar_angle_r', 2, 'none', ...
        'foot line 237'; ...
    'per_long_r', 'per_long_mif', {'tibia', 'calcnA', 'calcnSub'}, ...
        {'t2', 'per_long_wr1', 'per_long_wr2', 'per_long_wr3', ...
         'per_long_wr4', 'per_long_wr5', 'a3'}, 3, 0, ...
        'ankle', 'subtalar', 'subtalar_angle_r', 2, 'none', ...
        'foot line 242'; ...
    'per_tert_r', 'per_tert_mif', {'tibia', 'calcnA', 'calcnSub'}, ...
        {'t2', 'per_tert_wr1', 'a2'}, 2, 0, ...
        'ankle', 'subtalar', 'subtalar_angle_r', 2, 'none', ...
        'foot line 247'; ...
    'flex_dig_r', 'flex_dig_mif', {'tibia', 'calcnA', 'toes'}, ...
        {'t4', 'flex_dig_wr1', 'flex_dig_wr2', 'flex_dig_wr3', ...
         'flex_dig_wr4', 'flex_dig_wr5', 'flex_dig_wr6', 'a7'}, 2, 3, ...
        'ankle', 'mtp_calcnA', 'mtp_angle_r', 2, 'toeFolds', ...
        'toe file line 220; calcn rows folded +p_s; subtalar pinned'; ...
    'flex_hal_r', 'flex_hal_mif', {'tibia', 'calcnA', 'toes'}, ...
        {'t4', 'flex_hal_wr1', 'flex_hal_wr2', 'flex_hal_wr3', ...
         'flex_hal_wr4', 'flex_hal_wr5', 'a8'}, 2, 3, ...
        'ankle', 'mtp_calcnA', 'mtp_angle_r', 2, 'toeFolds', ...
        'toe file line 226; calcn rows folded +p_s; subtalar pinned'; ...
    'ext_dig_r',  'ext_dig_mif',  {'tibia', 'calcnA', 'toes'}, ...
        {'t2', 'ext_dig_wr1', 'ext_dig_wr2', 'ext_dig_wr3', ...
         'ext_dig_wr4', 'a5'}, 2, 2, ...
        'ankle', 'mtp_calcnA', 'mtp_angle_r', 2, 'toeFolds', ...
        'toe file line 234; calcn rows folded +p_s; subtalar pinned'; ...
    'ext_hal_r',  'ext_hal_mif',  {'tibia', 'calcnA', 'toes'}, ...
        {'t2', 'ext_hal_wr1', 'ext_hal_wr2', 'ext_hal_wr3', ...
         'ext_hal_wr4', 'ext_hal_wr5', 'a6'}, 2, 3, ...
        'ankle', 'mtp_calcnA', 'mtp_angle_r', 2, 'toeFolds', ...
        'toe file line 240; calcn rows folded +p_s; subtalar pinned'; ...
    'glut_max_r', 'glut_max_mif', {'pelvis', 'femurA', 'femurB'}, ...
        {'p6', 'glut_max_wr1', 'glut_max_wr2', 'glut_max_wr3', 'f3'}, 2, 0, ...
        'hip_adduction', 'hip_flexion_femurA', 'hip_flexion_r', 2, 'none', ...
        'uni_hip line 119; hip compound split adduction x flexion'; ...
    'add_mag_r',  'add_mag_mif',  {'pelvis', 'femurA', 'femurB'}, ...
        {'p9', 'add_mag_wr1', 'add_mag_wr2', 'f4'}, 2, 0, ...
        'hip_adduction', 'hip_flexion_femurA', 'hip_adduction_r', 1, ...
        'none', 'uni_hip line 150; adduction primary (coupled envelope)'; ...
    'iliacus_r',  'iliacus_mif',  {'pelvis', 'femurA', 'femurB'}, ...
        {'p1', 'iliacus_wr1', 'f4'}, 2, 0, ...
        'hip_adduction', 'hip_flexion_femurA', 'hip_flexion_r', 2, 'none', ...
        'uni_hip line 170'; ...
    'psoas_r',    'psoas_mif',    {'torso', 'pelvisShift', 'femur'}, ...
        {'b1', 'psoas_wr1', 'f4'}, 1, 1, ...
        'lumbar_extension_torsoView', 'hip_flexion_pelvisShift', ...
        'hip_flexion_r', 2, 'psoasFold', ...
        'uni_hip line 176; torso reference; pelvis row shifted -p_b'};

%% Actuator map (master names by id), if present
masterNameById = cell(1, 27);
mapFile = fullfile(iRepoRoot(fileparts(mfilename('fullpath'))), ...
    'Code', 'Matlab', 'Mesh_Optimization', 'Routing27_v2', 'actuator_map_27.json');
if exist(mapFile, 'file')
    try
        am = jsondecode(fileread(mapFile));
        for a = 1:numel(am.actuators)
            id = double(am.actuators(a).id);
            if id >= 1 && id <= 27
                masterNameById{id} = char(am.actuators(a).master_name);
            end
        end
    catch
        masterNameById = cell(1, 27);   % ids still carried from the route order
    end
end

%% Build one spec per requested name
nSpec = numel(names);
specs = struct([]);   %#ok<STRNU>
for k = 1:nSpec
    row = find(strcmp(routes(:, 1), names{k}), 1);
    if isempty(row)
        error('biPulleySpecsFromMasters_20261023:unknownMuscle', ...
            'No master''s route named %s. Available: %s', ...
            names{k}, strjoin(routes(:, 1), ', '))
    end

    pointNames = routes{row, 4};
    nPts = numel(pointNames);
    pts = zeros(nPts, 3);
    for j = 1:nPts
        if ~isfield(S, pointNames{j})
            error('biPulleySpecsFromMasters_20261023:missingPoint', ...
                'AttachPoints_Robot.mat has no variable %s (route %s)', ...
                pointNames{j}, names{k})
        end
        pts(j, :) = S.(pointNames{j}).';
    end

    % Optional per-row folds (toe calcn rows +p_s; psoas pelvis row -p_b).
    rowFolds = zeros(nPts, 3);
    foldNote = '';
    switch routes{row, 11}
        case 'toeFolds'
            nProx0 = routes{row, 5};
            nMid0 = routes{row, 6};
            rowFolds(nProx0 + 1:nProx0 + nMid0, :) = repmat(pSub, nMid0, 1);
            foldNote = sprintf('; calcn rows (rows %d-%d) folded by +p_s=%s', ...
                nProx0 + 1, nProx0 + nMid0, mat2str(pSub));
        case 'psoasFold'
            rowFolds(2, :) = -pBack;
            foldNote = sprintf('; pelvis row 2 shifted by -p_b=%s', ...
                mat2str(-pBack));
    end
    pts = pts + rowFolds;

    nProx = routes{row, 5};
    nMid = routes{row, 6};
    cross = [nProx + 1, nProx + nMid + 1];   % BiPamData Cross convention

    % Rigid-clone bookkeeping (knee-only routes): duplicate the insertion
    % row so crossing 2 exists; its outputs are inert (identity transform).
    isRigid2 = strcmp(routes{row, 8}, 'rigid');
    if isRigid2
        pts(end + 1, :) = pts(end, :);       %#ok<AGROW>
        rowFolds(end + 1, :) = rowFolds(end, :);
        cross(2) = size(pts, 1);
    end

    mifField = routes{row, 2};
    if ~isfield(P, mifField)
        error('biPulleySpecsFromMasters_20261023:missingMif', ...
            'PAM_Prop.mat has no variable %s (route %s)', mifField, names{k})
    end

    %% Hinge data from the catalog
    hingeKeys = {routes{row, 7}, routes{row, 8}};
    pivots = zeros(2, 3);
    axes = zeros(2, 3);
    types = {'hinge', 'hinge'};
    ranges = zeros(2, 2);
    dofs = cell(1, 2);
    kneeSlot = 0;
    for s = 1:2
        if isRigid2 && s == 2
            types{s} = 'fixed';
            pivots(s, :) = [0, 0, 0];        % identity transform
            axes(s, :) = [0, 0, 1];
            ranges(s, :) = [0, 0];
            dofs{s} = '';
            continue
        end
        h = J.(hingeKeys{s});
        axes(s, :) = h.axis;
        pivots(s, :) = h.pivot;
        ranges(s, :) = h.range;
        dofs{s} = h.dof;
        if h.rolling
            kneeSlot = s;
        end
    end

    % Per-route range blocks (kept from the original for the verbatim
    % routes; derived from the catalog otherwise).
    switch hingeKeys{1}
        case 'hip_flexion'
            thetaRanges = thetaHipKnee;       % original block, verbatim
        case 'knee'
            thetaRanges = thetaKneeAnkle;     % original block, verbatim
        otherwise
            thetaRanges = ranges;
    end
    if isRigid2
        thetaRanges(2, :) = [0, 0];          % the fixed clone never sweeps
    end

    spec = struct();
    spec.name = names{k};
    spec.source = 'masters_20261023';
    spec.points = pts;
    spec.cross = cross;
    spec.frames = routes{row, 3};
    spec.diameter = 20;

    spec.pulley = struct( ...
        'nPulleyBPA', cfgField(pulleyConfig, 'nPulleyBPA', 1), ...
        'tackleLineParts', cfgField(pulleyConfig, 'tackleLineParts', 1), ...
        'routingMode', cfgField(pulleyConfig, 'routingMode', 'moving_via'));
    spec.bundleOffsets = zeros(0, 3);
    spec.bundleMargin = 0.005;
    spec.footprintWidth = spec.diameter * 1e-3;

    spec.fmax = P.(mifField);
    spec.fmaxSource = sprintf('PAM_Prop.mat %s', mifField);
    spec.kin.pivots = pivots;             % 2x3, parent-frame coords
    spec.kin.axes = axes;                 % 2x3 unit rotation axes
    spec.kin.types = types;
    spec.kin.thetaRangesDeg = thetaRanges;
    spec.kin.dofs = dofs;                 % campaign coordinate name per hinge
    spec.kin.jointNames = hingeKeys;
    spec.kin.primarySlot = routes{row, 10};
    spec.kin.primaryDof = routes{row, 9};
    spec.kin.kneeSlot = kneeSlot;
    spec.kin.kneeFit = kneeFit;           % rolling knee center (femur coords)
    spec.kin.pointFolds = {[0; 0; 0], [0; 0; 0]};   % driver contract: none
    spec.kin.rowFolds = rowFolds;         % folds applied HERE (audit trail)
    % Gradeability self-check: the primary dof must be one of the hinges.
    if ~strcmp(dofs{spec.kin.primarySlot}, spec.kin.primaryDof)
        error('biPulleySpecsFromMasters_20261023:primarySlotMismatch', ...
            ['Route %s: primary dof %s is not modeled at slot %d ' ...
             '(hinges: %s / %s)'], names{k}, spec.kin.primaryDof, ...
            spec.kin.primarySlot, dofs{1}, dofs{2})
    end
    ix = find(strcmp(allNames, names{k}), 1);
    spec.actuatorId = allIds(ix);
    spec.masterName = masterNameById{spec.actuatorId};
    if isempty(spec.masterName)
        spec.masterName = '';
    end
    spec.notes = sprintf(['master''s route %s: points [%s]; hinges %s + %s ' ...
        '(slots 1+2); primary %s at slot %d; knee center from MovePoints ' ...
        'grids (rolling knee); %s%s%s'], names{k}, ...
        strjoin(pointNames, ' '), dofs{1}, dofs{2}, spec.kin.primaryDof, ...
        spec.kin.primarySlot, routes{row, 12}, foldNote, ...
        iIf(isRigid2, '; insertion row duplicated (rigid crossing 2)', ''));

    if isempty(specs)
        specs = spec;
    else
        specs(end + 1) = spec;   %#ok<AGROW>
    end
end

end

%% ---------------------------------------------------------------------
function [names, ids] = masterRouteNames_20261023()
% The full 27, actuator-map order (actuator ids 1..27).
names = { ...
    'ercspn_r', 'intobl_r', 'extobl_r', ...                       % 1-3
    'bifemlh_r', 'semimem_r', 'sar_r', 'tfl_r', 'rect_fem_r', 'grac_r', ...
    'bifemsh_r', 'vas_int_r', ...                                 % 4-11
    'soleus_r', 'med_gas_r', 'lat_gas_r', ...                     % 12-14
    'tib_post_r', 'tib_ant_r', 'per_brev_r', 'per_long_r', 'per_tert_r', ...  % 15-19
    'flex_dig_r', 'flex_hal_r', 'ext_dig_r', 'ext_hal_r', ...     % 20-23
    'glut_max_r', 'add_mag_r', 'iliacus_r', 'psoas_r'};           % 24-27
ids = 1:27;
end

%% ---------------------------------------------------------------------
function root = iRepoRoot(here)
root = here;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        return
    end
    if strcmp(parent, root) || isempty(parent)
        error('Could not locate the Bipedal_Robot repo root from %s', here)
    end
    root = parent;
end
end

%% ---------------------------------------------------------------------
function s = iIf(cond, a, b)
if cond
    s = a;
else
    s = b;
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
