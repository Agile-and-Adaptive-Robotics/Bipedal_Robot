% Hand-tuned BPA capability row, 2-D RE-GRADE (2026-10-23). FIX_SPEC_20261023
% section 5. Graded on eb475ws4.
%
% Re-grades the hand-tuned row (v1: ht03_grade_hand_tuned_20261022.m, 25-pt
% 1-D primary grids) on the SAME 2-D single-muscle surfaces as the optimized
% row: human_surfaces_5x5_grading_20261023.csv (single grading source of
% truth for BOTH rows). Per actuator: old-code routing (Connor-corrected
% where covered -- Connor's six section scripts run VERBATIM, exactly as v1)
% with OLD-CODE BPA sizes at 620 kPa, evaluated at each of the actuator's
% 5x5 (primary, secondary) nodes: primary node on the primary old-code
% hinge, secondary node on the crossed old-code hinge (the old code's own
% biarticular formulas; back/unihip two-axis frames use Connor's own Rx*Rz
% composition). Torque extracted about the PRIMARY dof (same MomentArm
% row/col as v1). The old code's soleus route has no knee dependence, so its
% surface is constant along the knee axis -- graded as-is (the old code's
% own formula set).
%
% METRICS per cell (identical booleans to v1 and to the optimized row's
% grade(), Opt_run_BiPulley_survey_20261022.m:524-527):
%   taut   = strain > 0 & isfinite(tauB);   strain = (LongestL-Lmt)/RestingL
%            (the old festo() contraction strain, = festo k exactly)
%   meet   = taut & (sign(tauB)==sign(tauH) | |tauH|<=1e-9)
%                  & |tauB| >= |tauH| - 1e-12
%   active = taut & |tauB| > 0
%   AREA SHARE = 100 * sum(metric & domain) / sum(domain). Cells with
%   in_graded_domain=0 (hip adduction-flexion coupling envelope, FIX_SPEC
%   1b) are excluded from the denominator for BOTH metrics and never counted
%   as failures. Flagged cells (via/step/slack/short) still count in the
%   headline percentage; the unflagged-only share is also reported
%   (FIX_SPEC section 2). Curves are never corrected.
%
% The nine v1 transcription-suspect flags + notes carry over VERBATIM from
% hand_tuned_metrics_20261022.csv; the fresh Phase D Connor-vs-Ben comparison
% must reproduce exactly the v1 mismatch set (asserted, else the run fails).
%
% GATES: (i) surfaces node grids verified vs map rom_table_deg on BOTH axes;
% (ii) Phase C T builders machine-exact (<=1e-12) vs the scripts' own
% TransformationMat (v1 verbatim); (iii) 2-D anchor: buildT2D at secondary=0
% must equal v1's 1-D buildGradedT machine-exact for every tcase; (iv)
% replay gate (FIX_SPEC 6, 0.000 tolerance): the full grading pass is
% recomputed from fresh T/Location/PamData objects and must match exactly.
%
% SOURCES (read-only): see ht03 header. OUTPUT:
% hand_tuned_metrics_20261023.csv (this folder). The repo and old-code files
% are never modified; v1 outputs stay untouched. No commits.

FAST = ~isempty(getenv('HT04_FAST'));  % HT04_FAST=1: Phase A divisions=25 (smoke only).

%% 0. Setup -------------------------------------------------------------------
host = getenv('COMPUTERNAME');
fprintf('HOST %s | %s | MATLAB %s\n', host, char(datetime('now')), version('-release'));
root  = 'D:\Github\Bipedal_Robot';
here  = fullfile(root, 'tmp', 'capability_table_20261022');
oldc  = fullfile(root, 'tmp', 'old_code_copy_20261022');
repoPOC = fullfile(root, 'Code', 'Matlab', 'Previous Optimization Code');
JOS   = fullfile(repoPOC, 'Joint_Optimization_Scripts');
addpath(fullfile(root, 'Code', 'Matlab', 'Functions'));  % maxBPAforce
addpath(repoPOC);                                        % JointData
addpath(oldc);                                           % ForceStrainTable.mat (festo)
addpath(here);                                           % staged PamData shadow FIRST
assert(~isempty(strfind(path, here)) && ~isempty(strfind(path, fullfile(root,'Code','Matlab','Functions'))));
w = which('PamData'); fprintf('PamData resolves to: %s\n', w);
w = which('PamData_20261022'); fprintf('PamData_20261022 resolves to: %s\n', w);
w = which('maxBPAforce'); fprintf('maxBPAforce resolves to: %s\n', w);
assert(strcmp(fileparts(which('PamData')), here), 'PamData shadow must win the path');

%% 1. Map + surfaces CSV + v1 metrics (carry-over) ----------------------------
mapPath = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Routing27_v2', 'actuator_map_27.json');
mj = jsondecode(fileread(mapPath));
acts = mj.actuators;
assert(numel(acts) == 27, 'map must have 27 actuators');

surfPath = fullfile(here, 'human_surfaces_5x5_grading_20261023.csv');
Sraw = readtable(surfPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char', 'CommentStyle', '#');
fprintf('Surfaces CSV: %d data rows.\n', height(Sraw));
assert(height(Sraw) == 27*25, 'surfaces CSV must carry 27x25 rows');

v1 = readtable(fullfile(here, 'hand_tuned_metrics_20261022.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
fprintf('v1 metrics CSV: %d rows (suspect flags carried over verbatim).\n', height(v1));

%% 2. PHASE A: run Connor's six section scripts verbatim ----------------------
% Sizing (Diameter/RestingL/LongestL) and the routing VALUES come straight
% from the scripts' own Muscles objects (zero transcription). v1 verbatim.
sec = struct();
divSet = 100; if FAST, divSet = 25; end

clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts Sraw v1 surfPath mapPath host
divisions = divSet;   % pre-set: scripts honor an existing 'divisions'
run(fullfile(JOS, 'RobotPAMCalculationOptimization_Back.m'));
sec.back = grabMuscles(Muscles); clear Muscles

clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts Sraw v1 surfPath mapPath host
divisions = divSet;
run(fullfile(JOS, 'RobotPAMCalculationOptimization_Bi_Hip.m'));
sec.bihip = grabMuscles(Muscles); clear Muscles

clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts Sraw v1 surfPath mapPath host
divisions = divSet;
run(fullfile(JOS, 'RobotPAMCalculationOptimization_Calves.m'));
sec.calves = grabMuscles(Muscles); clear Muscles   % {MG, LG, Soleus}

clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts Sraw v1 surfPath mapPath host
divisions = divSet;
run(fullfile(JOS, 'RobotPAMCalculationOptimization_Foot.m'));
sec.foot = grabMuscles(Muscles); clear Muscles     % {tib_post, tib_ant, per_brev, per_long, per_tert}

clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts Sraw v1 surfPath mapPath host
divisions = divSet;
run(fullfile(JOS, 'RobotPAMCalculationOptimization_Toe.m'));
sec.toe = grabMuscles(Muscles); clear Muscles      % {FDL, FHL, EDL, EHL}

% Uni_Hip: the script ERRORS at Muscle4 (PLocation undefined; map documents
% 'presumed Location4'). try/catch, then build Psoas from the script's own
% workspace values (Location4 etc. as left at the error point). v1 verbatim.
clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts Sraw v1 surfPath mapPath host
divisions = divSet;
uniErr = '';
try
    run(fullfile(JOS, 'RobotPAMCalculationOptimization_Uni_Hip.m'));
catch e
    uniErr = e.message;
end
if ~isempty(uniErr)
    fprintf('Uni_Hip script stopped (expected PLocation slip): %s\n', uniErr);
end
assert(exist('Location4', 'var') == 1 && exist('T', 'var') == 1, ...
    'Uni_Hip workspace incomplete after run');
Muscle4p = PamData_20261022('Psoas X', Location4, PCrossPoints, PMIF, T, Axis1);
sec.unihip = grabMuscles({Muscle1, Muscle2, Muscle3, Muscle4p});
clear Muscles Muscle1 Muscle2 Muscle3 Muscle4p Location4 PCrossPoints PMIF Axis1 T

fprintf('\n=== PHASE A: section muscles (Dia/Rest/Long from Connor''s scripts, divisions=%d) ===\n', divSet);
fn = fieldnames(sec);
for f = 1:numel(fn)
    S = sec.(fn{f});
    fprintf('%-8s:', fn{f});
    for m = 1:numel(S)
        fprintf('  %s[%.0fmm r=%.4g L=%.4g]', S(m).Muscle, S(m).Diameter, S(m).RestingL, S(m).LongestL);
    end
    fprintf('\n');
end

%% 3. Constants + builders (v1 verbatim) --------------------------------------
d2r = pi/180;
C.p_b = [-0.1007; 0.0815; 0];                    % connor_Back.m:27
C.p_h = [-0.0707; -0.0661; 0.0835];              % connor_Bi_Hip.m:27
C.p_a = [0; -0.43; 0];                           % connor_Calves.m:30
C.p_s = [-0.04877; -0.04195; 0.00792];           % connor_Calves.m:36
C.p_m = [0.1788; -0.002; 0.00108];               % connor_Toe.m:38
C.omg_a = [-0.10501355; -0.17402245; 0.97912632];% connor_Calves.m:27
C.omg_s = [0.7871796; 0.60474746; -0.12094949];  % connor_Calves.m:33
C.omg_m = [-0.5809544; 0; 0.81393611];           % connor_Toe.m:35
% knee + patella fits (connor_Bi_Hip.m:39-44,108-113,152-153, VERBATIM data)
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x       = [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
C.KneeXFunc = fit(knee_angle_x, knee_x, 'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y       = [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
C.KneeYFunc = fit(knee_angle_y, knee_y, 'cubicspline');
rect_fem_x  = [0.0156367; 0.0179948; 0.0274274; 0.029683; 0.0306; 0.0366; 0.0422; 0.0451; 0.0484; 0.0533; 0.0617; 0.0634; 0.067; 0.0733];
rect_fem_xD = pi/180*[-120.118; -114.871; -90.068; -83.532; -80; -60; -40; -30; -20; -10; 0; 1.6; 5; 10];
C.fcn3 = fit(rect_fem_xD, rect_fem_x, 'smoothingspline');
rect_fem_y  = [0.0234; 0.0238; 0.0251; 0.0253; 0.025284; 0.0249; 0.0243; 0.0239; 0.0234; 0.0228; 0.0210; 0.0206; 0.0192; 0.0160];
rect_fem_yD = pi/180*[-120; -114.6; -90; -83.5; -80.01; -60; -40; -30; -20; -10; 0; 1.6; 5; 10 ];
C.fcn4 = fit(rect_fem_yD, rect_fem_y, 'smoothingspline');
C.fcn7 = fit(rect_fem_xD, rect_fem_x-(0.0156367-0.008285), 'smoothingspline');
C.fcn8 = fit(rect_fem_yD, rect_fem_y+(0.0256239-0.02346), 'smoothingspline');

%% 4. PHASE C: verify the T builders + 3-D Location rebuilds (v1 verbatim) ----
fprintf('\n=== PHASE C: verification vs the scripts'' own T / Location ===\n');
D = divSet;
ver = struct('name', {}, 'diff', {});
Tb = buildSzT('back', D, C); Tb(:,:,2,:) = Tb(:,:,1,:); Tb(:,:,3,:) = Tb(:,:,1,:);  % connor_Back.m:55-56 duplication
ver(end+1) = struct('name', 'back T',       'diff', maxAbs(Tb, sec.back(1).TransformationMat));
ver(end+1) = struct('name', 'bihip T',      'diff', maxAbs(buildSzT('bihip', D, C), sec.bihip(1).TransformationMat));
Tcs = buildSzT('calves_soleus', D, C);
ver(end+1) = struct('name', 'calves soleus T (slots 1-2)', 'diff', ...
    maxAbs(Tcs, sec.calves(3).TransformationMat(:, :, 1:2, :)));
ver(end+1) = struct('name', 'calves gastroc T', 'diff', maxAbs(buildSzT('calves_gastroc', D, C), sec.calves(1).TransformationMat));
Tft = buildSzT('foot', D, C);
ver(end+1) = struct('name', 'foot T (slots 1-2)', 'diff', ...
    maxAbs(Tft, sec.foot(1).TransformationMat(:, :, 1:2, :)));
ver(end+1) = struct('name', 'toe T',        'diff', maxAbs(buildSzT('toe', D, C), sec.toe(1).TransformationMat));
ver(end+1) = struct('name', 'unihip T',     'diff', maxAbs(buildSzT('unihip', D, C), sec.unihip(1).TransformationMat));
ver(end+1) = struct('name', 'unihip psoas T', 'diff', maxAbs(buildSzT('unihip_psoas', D, C), sec.unihip(4).TransformationMat));
[LocR, ~, ~] = rebuildRectFem(sec.bihip(5), D, C);
ver(end+1) = struct('name', 'rect_fem Location rebuild', 'diff', maxAbs(LocR, sec.bihip(5).Location));
[LocV, ~, ~] = rebuildVasInt(sec.bihip(8), D, C);
ver(end+1) = struct('name', 'vas_int Location rebuild', 'diff', maxAbs(LocV, sec.bihip(8).Location));
for v = 1:numel(ver)
    fprintf('  %-32s max|diff| = %.3g %s\n', ver(v).name, ver(v).diff, ...
        ternary(ver(v).diff <= 1e-12, 'EXACT (machine noise)', '*** DIFFERS ***'));
end
assert(all([ver.diff] <= 1e-12), 'builder verification failed');

%% 5. Per-actuator spec: 2-D graded T case + dof pair --------------------------
% Same section/index/row/col/s/Ben-comparison as v1's Phase B table, plus the
% (primary_dof, secondary_dof) pair each 2-D case sweeps (asserted against
% the surfaces CSV and the map below).
APR = load(fullfile(oldc, 'AttachPoints_Robot.mat'));
PP  = load(fullfile(oldc, 'PAM_Prop.mat'));
E = [];
E = addSpec(E, 1,'back',1,'back2d',1,3,+1,{'p7','ercspn_wr1','ercspn_wr2'}, 'ercspn',  'lumbar_extension','lumbar_bending');
E = addSpec(E, 2,'back',2,'back2d',1,3,+1,{'p2','b3'}, 'intobl',  'lumbar_extension','lumbar_bending');
E = addSpec(E, 3,'back',3,'back2d',1,3,+1,{'p4','b4'}, 'extobl',  'lumbar_extension','lumbar_bending');
E = addSpec(E, 4,'bihip',1,'bihip2d',2,3,+1,{'p5','twr2','twr3','t5'}, 'bifemlh', 'knee_angle_r','hip_flexion_r');
E = addSpec(E, 5,'bihip',2,'bihip2d',2,3,+1,{'p5','semimem_wr1','semimem_wr2','twr1','t1','t6'}, 'semimem', 'knee_angle_r','hip_flexion_r');
E = addSpec(E, 6,'bihip',3,'bihip2d',2,3,+1,{'p3','sar_wr1','sar_wr2','sar_wr3','twr1','t1','t6'}, 'sar', 'knee_angle_r','hip_flexion_r');
E = addSpec(E, 7,'bihip',4,'bihip2d',2,3,+1,{'p3','tfl_wr1','tfl_wr2','t2','t5'}, 'tfl', 'knee_angle_r','hip_flexion_r');
E = addSpec(E, 8,'bihip',5,'bihip2dR',2,3,+1,{'p10','rect_fem_wr1','rect_fem_wr2','rect_fem_wr3'}, 'rect_fem', 'knee_angle_r','hip_flexion_r');
E = addSpec(E, 9,'bihip',6,'bihip2d',2,3,+1,{'p8','grac_wr1','grac_wr2','twr1','t1','t6'}, 'grac', 'knee_angle_r','hip_flexion_r');
E = addSpec(E,10,'bihip',7,'bihip2d',2,3,+1,{'f2','twr2','twr3'}, 'bifemsh', 'knee_angle_r','hip_flexion_r');
E = addSpec(E,11,'bihip',8,'bihip2dV',2,3,+1,{'f1','vas_int_wr1','vas_int_wr2'}, 'vas_int', 'knee_angle_r','hip_flexion_r');
E = addSpec(E,12,'calves',3,'soleus2d',1,1,+1,{'t1','a1'}, 'soleus', 'ankle_angle_r','knee_angle_r');
E = addSpec(E,13,'calves',1,'gastroc2d',2,1,+1,{'f4','med_gas_wr1','a1'}, 'med_gas', 'ankle_angle_r','knee_angle_r');
E = addSpec(E,14,'calves',2,'gastroc2d',2,1,+1,{'f5','lat_gas_wr1','a1'}, 'lat_gas', 'ankle_angle_r','knee_angle_r');
E = addSpec(E,15,'foot',1,'foot2d',2,1,+1,{'t4','tib_post_wr1','tib_post_wr2','a4'}, 'tib_post', 'subtalar_angle_r','ankle_angle_r');
E = addSpec(E,16,'foot',2,'tibant2d',1,1,+1,{'t2','tib_ant_wr1','a3'}, 'tib_ant', 'ankle_angle_r','subtalar_angle_r');
E = addSpec(E,17,'foot',3,'foot2d',2,1,+1,{'t2','per_brev_wr1','per_brev_wr2','per_brev_wr3','a2'}, 'per_brev', 'subtalar_angle_r','ankle_angle_r');
E = addSpec(E,18,'foot',4,'foot2d',2,1,+1,{'t2','per_long_wr1','per_long_wr2','per_long_wr3','per_long_wr4','per_long_wr5','a3'}, 'per_long', 'subtalar_angle_r','ankle_angle_r');
E = addSpec(E,19,'foot',5,'foot2d',2,1,+1,{'t2','per_tert_wr1','a2'}, 'per_tert', 'subtalar_angle_r','ankle_angle_r');
E = addSpec(E,20,'toe',1,'toe2d',3,1,+1,{'t4','flex_dig_wr1','flex_dig_wr2','flex_dig_wr3','flex_dig_wr4','flex_dig_wr5','flex_dig_wr6','a7'}, 'flex_dig', 'mtp_angle_r','ankle_angle_r');
E = addSpec(E,21,'toe',2,'toe2d',3,1,+1,{'t4','flex_hal_wr1','flex_hal_wr2','flex_hal_wr3','flex_hal_wr4','flex_hal_wr5','a8'}, 'flex_hal', 'mtp_angle_r','ankle_angle_r');
E = addSpec(E,22,'toe',3,'toe2d',3,1,+1,{'t2','ext_dig_wr1','ext_dig_wr2','ext_dig_wr3','ext_dig_wr4','a5'}, 'ext_dig', 'mtp_angle_r','ankle_angle_r');
E = addSpec(E,23,'toe',4,'toe2d',3,1,+1,{'t2','ext_hal_wr1','ext_hal_wr2','ext_hal_wr3','ext_hal_wr4','ext_hal_wr5','a6'}, 'ext_hal', 'mtp_angle_r','ankle_angle_r');
E = addSpec(E,24,'unihip',1,'unihip2dZ',1,3,+1,{'p6','glut_max_wr1','glut_max_wr2','glut_max_wr3','f3'}, 'glut_max', 'hip_flexion_r','hip_adduction_r');
E = addSpec(E,25,'unihip',2,'unihip2dX',1,1,+1,{'p9','add_mag_wr1','add_mag_wr2','f4'}, 'add_mag', 'hip_adduction_r','hip_flexion_r');
E = addSpec(E,26,'unihip',3,'unihip2dZ',1,3,+1,{'p1','iliacus_wr1','f4'}, 'iliacus', 'hip_flexion_r','hip_adduction_r');
E = addSpec(E,27,'unihip',4,'psoas2d',2,3,+1,{'b1','psoas_wr1','f4'}, 'psoas', 'hip_flexion_r','hip_adduction_r');
assert(numel(E) == 27);

%% 6. Surfaces: per-actuator grid parse + verification -------------------------
fprintf('\n=== Surfaces verification (grids vs map rom_table_deg, dof pairs) ===\n');
surf = cell(27, 1);
gridOK = true;
idn = coerceNum(Sraw.actuator_id);
csvPriAll = coerceCell(Sraw.primary_dof);
csvSecAll = coerceCell(Sraw.secondary_dof);
v1muscle = coerceCell(v1.muscle);
v1source = coerceCell(v1.source);
v1sus    = coerceCell(v1.suspect);
v1note   = coerceCell(v1.note);
for k = 1:27
    a = acts(k); e = E(k);
    S = struct();
    S.p    = coerceNum(Sraw.p_node_deg(idn == a.id));
    S.s    = coerceNum(Sraw.s_node_deg(idn == a.id));
    S.tauH = coerceNum(Sraw.tau_h_Nm(idn == a.id));
    S.dom  = coerceNum(Sraw.in_graded_domain(idn == a.id)) == 1;
    S.flg  = (coerceNum(Sraw.flag_via(idn == a.id)) == 1) | ...
             (coerceNum(Sraw.flag_step(idn == a.id)) == 1) | ...
             (coerceNum(Sraw.flag_slack(idn == a.id)) == 1) | ...
             (coerceNum(Sraw.flag_short(idn == a.id)) == 1);
    assert(numel(S.p) == 25 && numel(S.tauH) == 25, 'act %d: need 25 surface rows', a.id);
    assert(all(isfinite(S.tauH)), 'act %d: non-finite tau_h in surfaces', a.id);
    % dof names: CSV vs spec table vs map (primary)
    csvPri = csvPriAll(idn == a.id);
    csvSec = csvSecAll(idn == a.id);
    assert(all(strcmp(csvPri, e.pri_dof)) && all(strcmp(csvSec, e.sec_dof)), ...
        'act %d: CSV dofs (%s,%s) != spec (%s,%s)', a.id, csvPri{1}, csvSec{1}, e.pri_dof, e.sec_dof);
    assert(strcmp(csvPri{1}, a.primary_dof), 'act %d: CSV primary != map primary', a.id);
    % node grids: 5 inclusive nodes per axis over the map rom rectangle
    [pNodes, ~, pI] = unique(round(S.p, 9));   % p outer, s inner (CSV order)
    [sNodes, ~, sI] = unique(round(S.s, 9));
    pI = pI(:)'; sI = sI(:)'; pNodes = pNodes(:)'; sNodes = sNodes(:)';   % force rows
    assert(numel(pNodes) == 5 && numel(sNodes) == 5, 'act %d: not a 5x5 grid', a.id);
    rp = mj.rom_table_deg.(e.pri_dof);
    rs = mj.rom_table_deg.(e.sec_dof);
    idealP = rp(1) + (rp(2) - rp(1)) * (0:4) / 4;
    idealS = rs(1) + (rs(2) - rs(1)) * (0:4) / 4;
    if max(abs(sort(pNodes) - idealP)) > 1e-4 || max(abs(sort(sNodes) - idealS)) > 1e-4
        gridOK = false;
        fprintf('GRID MISMATCH act %d: p [%s] vs [%s], s [%s] vs [%s]\n', a.id, ...
            mat2str(sort(pNodes'), 6), mat2str(idealP, 6), mat2str(sort(sNodes'), 6), mat2str(idealS, 6));
    end
    % canonical pose order: p outer (i), s inner (j) -> pose (i-1)*5+j
    S.pOrder = pNodes; S.sOrder = sNodes;
    S.ii = pI; S.jj = sI;
    S.pose = (pI - 1) * 5 + sI;                  % 1..25 canonical
    assert(isequal(sort(S.pose), 1:25), 'act %d: incomplete 5x5 cross product', a.id);
    [S.canon, orderIdx] = sort(S.pose);
    S.pC = S.p(orderIdx); S.sC = S.s(orderIdx); S.tauHC = S.tauH(orderIdx);
    S.domC = S.dom(orderIdx); S.flgC = S.flg(orderIdx);
    surf{k} = S;
    if ~all(S.dom), fprintf('  act %2d %-32s envelope-limited: %d/25 cells in graded domain (%d excluded)\n', ...
        a.id, a.master_name, sum(S.dom), sum(~S.dom)); end
end
fprintf('Grid/dof verification vs map rom_table_deg (both axes, all 27): %s\n', string(gridOK));
assert(gridOK, 'surfaces grid verification failed');
nDomAll = cellfun(@(S) sum(S.domC), surf);
fprintf('Domain cells: %d total (%d excluded by the hip coupling envelope, acts 24-27).\n', ...
    sum(nDomAll), 27*25 - sum(nDomAll));

%% 7. 2-D T anchor: buildT2D at secondary=0 == v1 1-D buildGradedT -------------
fprintf('\n=== 2-D anchor: buildT2D(thP, 0) == v1 buildGradedT(thP) ===\n');
anchorCases = {'back2d','back'; 'bihip2d','knee'; 'bihip2dR','knee3dR'; 'bihip2dV','knee3dV'; ...
    'soleus2d','ankle2'; 'gastroc2d','ankle3'; 'foot2d','subtal2'; 'tibant2d','ankle2'; ...
    'toe2d','mtp3'; 'unihip2dZ','hipflex1'; 'unihip2dX','hipadd1'; 'psoas2d','hipflex2'};
thTry = [-2.0, -0.35, 0.0, 0.47, 1.1];   % rad, within all swept ranges
anchorOK = true;
for c = 1:size(anchorCases, 1)
    d = maxAbs(buildT2D(anchorCases{c,1}, thTry, zeros(size(thTry)), C), ...
               buildGradedT(anchorCases{c,2}, thTry, C));
    ok = d <= 1e-12;
    anchorOK = anchorOK & ok;
    fprintf('  %-10s vs v1 %-9s max|diff| = %.3g %s\n', anchorCases{c,1}, anchorCases{c,2}, d, ...
        ternary(ok, 'EXACT (machine noise)', '*** DIFFERS ***'));
end
assert(anchorOK, '2-D T anchor against v1 failed');

%% 8. PHASE B: 2-D grading pass on the surfaces grids --------------------------
fprintf('\n=== PHASE B: 2-D grading on the surfaces 5x5 grids (620 kPa, domain-masked) ===\n');
res = [];
for k = 1:27
    a = acts(k); e = E(k); mo = sec.(e.section)(e.idx); S = surf{k};
    assert(isnumeric(mo.RestingL) && isnumeric(mo.Diameter), ...
        'act %d: non-numeric BPA size from Size()', a.id);
    g = gradeAct(e, mo, S, C);                 % canonical pose order (p outer, s inner)
    tauB = g.tauB; strain = g.strain; tauH = S.tauHC;
    dom = S.domC; flg = S.flgC;
    taut  = strain > 0 & isfinite(tauB);
    meet  = taut & isfinite(tauH) & (sign(tauB) == sign(tauH) | abs(tauH) <= 1e-9) ...
            & abs(tauB) >= abs(tauH) - 1e-12;
    activ = taut & abs(tauB) > 0;
    nDom = sum(dom); nEx = sum(~dom); nFlg = sum(dom & flg);
    assert(nDom > 0, 'act %d: empty graded domain', a.id);
    meetPct = 100 * sum(meet & dom) / nDom;
    actPct  = 100 * sum(activ & dom) / nDom;
    du = dom & ~flg;
    if any(du)
        meetPctU = 100 * sum(meet & du) / sum(du);
        actPctU  = 100 * sum(activ & du) / sum(du);
    else
        meetPctU = NaN; actPctU = NaN;
    end
    mSk = dom & taut & abs(tauH) > 1;
    if any(mSk), signAgree = 100 * mean(sign(tauB(mSk)) == sign(tauH(mSk))); else, signAgree = NaN; end
    r = struct('id', a.id, 'name', a.master_name, 'pri', e.pri_dof, 'sec', e.sec_dof, ...
        'dia', mo.Diameter, 'rest', mo.RestingL, 'long', mo.LongestL, 'F', g.F, ...
        'strainMin', min(strain), 'strainMax', max(strain), ...
        'tauMin', min(tauB), 'tauMax', max(tauB), ...
        'nDom', nDom, 'nEx', nEx, 'nFlg', nFlg, ...
        'meetPct', meetPct, 'actPct', actPct, 'meetPctU', meetPctU, 'actPctU', actPctU, ...
        'signAgree', signAgree, 'tauB', tauB, 'strain', strain, 'meet', meet, 'activ', activ);
    if k == 1, res = r; else, res(k) = r; end %#ok<AGROW>
    fprintf(['%2d %-32s [%s x %s] d%2.0fmm r=%.4g F=%7.1fN dom %2d/25 (excl %d, flg %d) ' ...
        'strain[%+.3f,%+.3f] tau[%+8.2f,%+8.2f] meet=%5.1f%% active=%5.1f%% ' ...
        '(unflg %5.1f%%/%5.1f%%) sign=%s\n'], ...
        r.id, r.name, e.pri_dof, e.sec_dof, r.dia, r.rest, r.F, nDom, nEx, nFlg, ...
        r.strainMin, r.strainMax, r.tauMin, r.tauMax, meetPct, actPct, ...
        meetPctU, actPctU, formatAgree(signAgree));
end

%% 9. REPLAY GATE (FIX_SPEC 6, 0.000 tolerance) --------------------------------
fprintf('\n=== REPLAY GATE: full grading pass recomputed from fresh objects ===\n');
replayOK = true;
for k = 1:27
    e = E(k); mo = sec.(e.section)(e.idx); S = surf{k};
    g2 = gradeAct(e, mo, S, C);
    d = max(maxAbs(g2.tauB, res(k).tauB), maxAbs(g2.strain, res(k).strain));
    if d > 0
        replayOK = false;
        fprintf('  act %d REPLAY DIFFERS: max|diff| = %.3g\n', k, d);
    end
end
fprintf('Replay gate (0.000 tolerance, all 27, fresh T/Location/PamData): %s\n', ...
    string(replayOK));
assert(replayOK, 'replay gate failed');

%% 10. PHASE D: Connor-vs-Ben value comparison + MIF sanity (v1 verbatim) ------
fprintf('\n=== PHASE D: Connor vs Ben (AttachPoints_Robot.mat) + MIF sanity ===\n');
cmp = cell(27, 1);
for k = 1:27
    e = E(k); mo = sec.(e.section)(e.idx);
    mm = struct('col', {}, 'ben', {}, 'connor', {}, 'ben_v', {});
    if isnumeric(mo.Location)
        nCols = size(mo.Location, 2);
        if ndims(mo.Location) == 2
            for c = 1:min(nCols, numel(e.ben))
                d = max(abs(mo.Location(:, c) - APR.(e.ben{c})));
                if d > 1e-12
                    mm(end+1) = struct('col', c, 'ben', e.ben{c}, ...
                        'connor', mo.Location(:, c), 'ben_v', APR.(e.ben{c})); %#ok<AGROW>
                end
            end
        else % 3-D: constant columns only (dynamic column covered by fit check)
            constCols = 1:(size(mo.Location, 2) - 1);
            for c = constCols
                d = max(abs(mo.Location(:, c, 1) - APR.(e.ben{c})));
                if d > 1e-12
                    mm(end+1) = struct('col', c, 'ben', e.ben{c}, ...
                        'connor', mo.Location(:, c, 1), 'ben_v', APR.(e.ben{c})); %#ok<AGROW>
                end
            end
        end
    end
    mifMap = acts(k).master_mif_n;
    mifBen = PP.([e.mifname '_mif']);
    cmp{k} = struct('mm', {mm}, 'mifMap', mifMap, 'mifBen', mifBen, 'mifConnor', mo.MIF);
    if isempty(mm)
        fprintf('%2d %-32s values MATCH Ben; MIF Connor %g / map %g / Ben %g\n', ...
            acts(k).id, acts(k).master_name, mo.MIF, mifMap, mifBen);
    else
        for j = 1:numel(mm)
            fprintf('%2d %-32s MISMATCH col %d vs Ben %s: Connor [%s] Ben [%s]\n', ...
                acts(k).id, acts(k).master_name, mm(j).col, mm(j).ben, ...
                mat2str(mm(j).connor', 4), mat2str(mm(j).ben_v', 4));
        end
        fprintf('%2d %-32s MIF Connor %g / map %g / Ben %g\n', ...
            acts(k).id, acts(k).master_name, mo.MIF, mifMap, mifBen);
    end
end
amLoc = sec.unihip(2).Location;
fprintf('add_mag Location columns 2 and 3 identical: %s\n', ...
    string(isequal(amLoc(:,2), amLoc(:,3))));
fprintf('ES: Ben originals carry 9 points (p7, ercspn_wr1..7, b2); Connor keeps 3.\n');
fprintf('ES: Ben hand-set ercspn_rest = 0.7428 (bipedal_kinematics_back.m:311); Connor flow uses Size output %.6g.\n', sec.back(1).RestingL);

% suspect set: fresh machine mismatches + v1's hand-coded cases (15, 25, 27),
% asserted IDENTICAL to the v1 suspect column (the nine flags carry over).
freshSuspect = false(1, 27);
for k = 1:27
    freshSuspect(k) = ~isempty(cmp{k}.mm) || ismember(k, [15 25 27]);
end
v1suspectName = v1muscle(strcmp(v1sus, 'yes'));
mySuspectName = {acts(freshSuspect).master_name};
fprintf('\nv1 suspects (%d): %s\n', numel(v1suspectName), strjoin(v1suspectName, '; '));
fprintf('fresh suspects (%d): %s\n', numel(mySuspectName), strjoin(mySuspectName, '; '));
assert(isequal(sort(v1suspectName), sort(mySuspectName)), ...
    'fresh Phase D suspect set differs from v1 -- investigate before trusting');

%% 11. Output CSV ---------------------------------------------------------------
if FAST
    outPath = fullfile(here, 'ht04_SMOKE_hand_tuned_2d_20261023.csv');
else
    outPath = fullfile(here, 'hand_tuned_metrics_20261023.csv');
end
fid = fopen(outPath, 'w');
fprintf(fid, '# hand-tuned metrics, 2-D re-grade 20261023 (FIX_SPEC_20261023 sections 1/1b/5)\n');
fprintf(fid, '# surfaces: human_surfaces_5x5_grading_20261023.csv (5x5 two-DOF grading grid; single source of truth for both rows)\n');
fprintf(fid, '# robot: old-code routing (Connor-corrected, Connor''s six section scripts verbatim) + old-code BPA sizes at 620 kPa\n');
fprintf(fid, '# meet/active per cell = v1/optimized-row booleans; area share over in_graded_domain cells only (hip-envelope cells excluded from denominator)\n');
fprintf(fid, '# unflagged = in-domain cells with no via/step/slack/short flag; suspect+note carried over VERBATIM from hand_tuned_metrics_20261022.csv\n');
fprintf(fid, 'muscle,source,primary_dof,secondary_dof,cells_domain,cells_excluded,cells_flagged_in_domain,meetPct,activePct,meetPct_unflagged,activePct_unflagged,suspect,note\n');
for k = 1:27
    r = res(k);
    ix = strcmp(v1muscle, r.name);
    assert(sum(ix) == 1, 'act %d: v1 row not found for %s', k, r.name);
    fprintf(fid, '"%s",%s,%s,%s,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%s,"%s"\n', ...
        strrep(r.name, '"', ''''), v1source{ix}, r.pri, r.sec, ...
        r.nDom, r.nEx, r.nFlg, r.meetPct, r.actPct, r.meetPctU, r.actPctU, ...
        v1sus{ix}, strrep(v1note{ix}, '"', ''''));
end
fclose(fid);
fprintf('\nWROTE %s\n', outPath);
nSus = sum(strcmp(v1sus, 'yes'));
fprintf('Coverage: 27/27 connor-corrected, 0 ben-original, 0 excluded (Connor''s six section scripts define all 27 master actuators).\n');
fprintf('Suspects flagged (v1, carried over unchanged): %d\n', nSus);
fprintf('Mean meetPct (domain-weighted cells): ');
wSum = 0; wMeet = 0; wAct = 0;
for k = 1:27
    wSum = wSum + res(k).nDom;
    wMeet = wMeet + res(k).nDom * res(k).meetPct / 100;
    wAct = wAct + res(k).nDom * res(k).actPct / 100;
end
fprintf('%.2f%% | mean activePct: %.2f%%  (over %d domain cells of %d)\n', ...
    100 * wMeet / wSum, 100 * wAct / wSum, wSum, 27*25);
fprintf('GRADE-2D DONE\n');

%% ---------------- local functions -------------------------------------------
function S = grabMuscles(M)
S = struct('Muscle', M{1}.Muscle, 'Location', M{1}.Location, ...
    'CrossPoints', M{1}.CrossPoints, 'MIF', M{1}.MIF, 'Axis', M{1}.Axis, ...
    'Diameter', M{1}.Diameter, 'RestingL', M{1}.RestingL, 'LongestL', M{1}.LongestL, ...
    'MuscleLength', M{1}.MuscleLength, 'MomentArm', M{1}.MomentArm, ...
    'Torque', M{1}.Torque, 'TransformationMat', M{1}.TransformationMat);
for k = 2:numel(M)
    S(k) = struct('Muscle', M{k}.Muscle, 'Location', M{k}.Location, ...
        'CrossPoints', M{k}.CrossPoints, 'MIF', M{k}.MIF, 'Axis', M{k}.Axis, ...
        'Diameter', M{k}.Diameter, 'RestingL', M{k}.RestingL, 'LongestL', M{k}.LongestL, ...
        'MuscleLength', M{k}.MuscleLength, 'MomentArm', M{k}.MomentArm, ...
        'Torque', M{k}.Torque, 'TransformationMat', M{k}.TransformationMat); %#ok<AGROW>
end
end

function x = coerceNum(v)
% readtable column -> double row vector (numeric or char cells).
if isnumeric(v), x = double(v(:))';
elseif iscell(v), x = str2double(v(:))';
elseif isstring(v), x = str2double(v(:))';
else, error('coerceNum: unsupported type');
end
end

function c = coerceCell(v)
% readtable text column -> cellstr row vector (cell, string, or char).
if iscell(v), c = v(:)';
elseif isstring(v), c = cellstr(v(:))';
elseif ischar(v), c = cellstr(v)';   % char matrix: one row per entry
else, error('coerceCell: unsupported type');
end
end

function R = Rax(ax, th)
% Rodrigues rotation (JointData.computeRotMat formula, unit axes).
R = [cos(th)+ax(1)^2*(1-cos(th)), ax(1)*ax(2)*(1-cos(th))-ax(3)*sin(th), ax(1)*ax(3)*(1-cos(th))+ax(2)*sin(th);
     ax(1)*ax(2)*(1-cos(th))+ax(3)*sin(th), cos(th)+ax(2)^2*(1-cos(th)), ax(2)*ax(3)*(1-cos(th))-ax(1)*sin(th);
     ax(1)*ax(3)*(1-cos(th))-ax(2)*sin(th), ax(2)*ax(3)*(1-cos(th))+ax(1)*sin(th), cos(th)+ax(3)^2*(1-cos(th))];
end

function T4 = tr4(R, p)
T4 = [R, p; 0 0 0 1];
end

function kh = kneeHomeOf(th, C)
kh = [C.KneeXFunc(th); C.KneeYFunc(th); 0];
end

function T = buildGradedT(caseName, th, C)
% v1's 1-D graded T (ht03 buildGradedT, VERBATIM): primary joint swept over
% th (rad), all other joints at default. Used only for the 2-D anchor gate.
n = numel(th);
switch caseName
    case 'back'
        T = zeros(4,4,1,n);
        for i = 1:n, T(:,:,1,i) = tr4(Rax([0;0;1], th(i)), C.p_b); end
    case {'knee', 'knee3dR', 'knee3dV'}
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(eye(3), C.p_h);
            T(:,:,2,i) = tr4(Rax([0;0;1], th(i)), kneeHomeOf(th(i), C));
        end
    case 'ankle2'
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(Rax(C.omg_a, th(i)), C.p_a);
            T(:,:,2,i) = tr4(eye(3), C.p_s);
        end
    case 'ankle3'
        T = zeros(4,4,3,n);
        kh0 = kneeHomeOf(0, C);
        for i = 1:n
            T(:,:,1,i) = tr4(eye(3), kh0);
            T(:,:,2,i) = tr4(Rax(C.omg_a, th(i)), C.p_a);
            T(:,:,3,i) = tr4(eye(3), C.p_s);
        end
    case 'subtal2'
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(eye(3), C.p_a);
            T(:,:,2,i) = tr4(Rax(C.omg_s, th(i)), C.p_s);
        end
    case 'mtp3'
        T = zeros(4,4,3,n);
        for i = 1:n
            T(:,:,1,i) = tr4(eye(3), C.p_a);
            T(:,:,2,i) = tr4(eye(3), C.p_s);
            T(:,:,3,i) = tr4(Rax(C.omg_m, th(i)), C.p_m);
        end
    case 'hipflex1'
        T = zeros(4,4,1,n);
        for i = 1:n, T(:,:,1,i) = tr4(Rax([0;0;1], th(i)), C.p_h); end
    case 'hipadd1'
        T = zeros(4,4,1,n);
        for i = 1:n, T(:,:,1,i) = tr4(Rax([1;0;0], th(i)), C.p_h); end
    case 'hipflex2'
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(eye(3), C.p_b);
            T(:,:,2,i) = tr4(Rax([0;0;1], th(i)), C.p_h);
        end
    otherwise
        error('bad graded case %s', caseName);
end
end

function T = buildT2D(caseName, thP, thS, C)
% 2-D graded T (NEW, 2026-10-23): primary node thP on the primary old-code
% hinge, secondary node thS on the crossed hinge; remaining modeled joints
% identity/default. Two-axis frames follow Connor's own composition order
% verbatim: back Rx*Rz (connor_Back.m:44-56), unihip hip Rx*Rz
% (connor_Uni_Hip.m:44-53). Anchor-gated against buildGradedT at thS=0.
n = numel(thP);
assert(numel(thS) == n, 'buildT2D: angle vectors must match');
switch caseName
    case 'back2d'      % primary lumbar extension (z), secondary bending (x); 1 slot
        T = zeros(4,4,1,n);
        for i = 1:n
            T(:,:,1,i) = tr4(Rax([1;0;0], thS(i)) * Rax([0;0;1], thP(i)), C.p_b);
        end
    case {'bihip2d', 'bihip2dR', 'bihip2dV'}  % primary knee (slot 2), secondary hip flexion (slot 1)
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(Rax([0;0;1], thS(i)), C.p_h);
            T(:,:,2,i) = tr4(Rax([0;0;1], thP(i)), kneeHomeOf(thP(i), C));
        end
    case 'soleus2d'    % primary ankle (slot 1); subtalar identity. The old
                       % code's soleus route has NO knee slot: the surface's
                       % secondary (knee) axis carries no coupling (as-is).
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(Rax(C.omg_a, thP(i)), C.p_a);
            T(:,:,2,i) = tr4(eye(3), C.p_s);
        end
    case 'gastroc2d'   % primary ankle (slot 2), secondary knee (slot 1), subtalar identity
        T = zeros(4,4,3,n);
        for i = 1:n
            T(:,:,1,i) = tr4(Rax([0;0;1], thS(i)), kneeHomeOf(thS(i), C));
            T(:,:,2,i) = tr4(Rax(C.omg_a, thP(i)), C.p_a);
            T(:,:,3,i) = tr4(eye(3), C.p_s);
        end
    case 'foot2d'      % primary subtalar (slot 2), secondary ankle (slot 1)
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(Rax(C.omg_a, thS(i)), C.p_a);
            T(:,:,2,i) = tr4(Rax(C.omg_s, thP(i)), C.p_s);
        end
    case 'tibant2d'    % primary ankle (slot 1), secondary subtalar (slot 2)
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(Rax(C.omg_a, thP(i)), C.p_a);
            T(:,:,2,i) = tr4(Rax(C.omg_s, thS(i)), C.p_s);
        end
    case 'toe2d'       % primary mtp (slot 3), secondary ankle (slot 1), subtalar identity (pinned 0)
        T = zeros(4,4,3,n);
        for i = 1:n
            T(:,:,1,i) = tr4(Rax(C.omg_a, thS(i)), C.p_a);
            T(:,:,2,i) = tr4(eye(3), C.p_s);
            T(:,:,3,i) = tr4(Rax(C.omg_m, thP(i)), C.p_m);
        end
    case 'unihip2dZ'   % single hip frame Rx(adduction)*Rz(flexion); PRIMARY flexion (z)
        T = zeros(4,4,1,n);
        for i = 1:n
            T(:,:,1,i) = tr4(Rax([1;0;0], thS(i)) * Rax([0;0;1], thP(i)), C.p_h);
        end
    case 'unihip2dX'   % same frame; PRIMARY adduction (x)
        T = zeros(4,4,1,n);
        for i = 1:n
            T(:,:,1,i) = tr4(Rax([1;0;0], thP(i)) * Rax([0;0;1], thS(i)), C.p_h);
        end
    case 'psoas2d'     % back slot identity (lumbar 0), hip composed; PRIMARY flexion (z)
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(eye(3), C.p_b);
            T(:,:,2,i) = tr4(Rax([1;0;0], thS(i)) * Rax([0;0;1], thP(i)), C.p_h);
        end
    otherwise
        error('bad 2-D case %s', caseName);
end
end

function g = gradeAct(e, mo, S, C)
% One actuator's 2-D grading pass: fresh T + Location + PamData object.
% Canonical pose order: p outer (i), s inner (j); knee-angle-dependent
% 3-D Locations rebuilt at the pose's knee angle (knee = primary for
% rect_fem/vas_int). tauB in OpenSim sign (e.s = +1, v1 convention).
d2r = pi/180;
thP = S.pC * d2r;
thS = S.sC * d2r;
T = buildT2D(e.tcase, thP, thS, C);
switch e.tcase
    case 'bihip2dR', Loc = rebuildRectFemG(mo, thP, C);
    case 'bihip2dV', Loc = rebuildVasIntG(mo, thP, C);
    otherwise,      Loc = mo.Location;
end
mG = PamData_20261022(mo.Muscle, Loc, mo.CrossPoints, mo.MIF, T, mo.Axis);
g.Lmt  = reshape(mG.MuscleLength, 1, []);
g.MA   = reshape(mG.MomentArm(e.row, :, e.col), 1, []);
g.F    = maxBPAforce(mo.RestingL, num2str(mo.Diameter), 620);
g.strain = (mo.LongestL - g.Lmt) / mo.RestingL;   % festo contraction strain (v1)
g.tauB = e.s * (g.MA * g.F);                       % OpenSim-sign torque
end

function T = buildSzT(caseName, D, C)
% Connor's own section grids (divisions = D), pose order ii = (i-1)*D + j.
% v1 verbatim.
switch caseName
    case 'back'
        thX = linspace(-15, 15, D) * pi/180;
        thZ = linspace(-45, 15, D) * pi/180;
        T = zeros(4,4,1,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax([1;0;0], thX(i)) * Rax([0;0;1], thZ(j)), C.p_b);
            end
        end
    case 'bihip'
        thH = linspace(-15, 110, D) * pi/180;
        thK = linspace(-2.0943951, 0.17453293, D);
        T = zeros(4,4,2,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax([0;0;1], thH(i)), C.p_h);
                T(:,:,2,ii) = tr4(Rax([0;0;1], thK(j)), kneeHomeOf(thK(j), C));
            end
        end
    case 'calves_soleus'
        thA = linspace(-50, 20, D) * pi/180;
        T = zeros(4,4,2,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax(C.omg_a, thA(i)), C.p_a);
                T(:,:,2,ii) = tr4(eye(3), C.p_s);
            end
        end
    case 'calves_gastroc'
        thA = linspace(-50, 20, D) * pi/180;
        thK = linspace(-2.0943951, 0.17453293, D);
        T = zeros(4,4,3,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax([0;0;1], thK(j)), kneeHomeOf(thK(j), C));
                T(:,:,2,ii) = tr4(Rax(C.omg_a, thA(i)), C.p_a);
                T(:,:,3,ii) = tr4(eye(3), C.p_s);
            end
        end
    case 'foot'
        thA = linspace(-50, 20, D) * pi/180;
        thS = linspace(-25, 35, D) * pi/180;
        T = zeros(4,4,2,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax(C.omg_a, thA(i)), C.p_a);
                T(:,:,2,ii) = tr4(Rax(C.omg_s, thS(j)), C.p_s);
            end
        end
    case 'toe'
        thA = linspace(-50, 20, D) * pi/180;
        thM = linspace(-30, 80, D) * pi/180;
        T = zeros(4,4,3,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax(C.omg_a, thA(i)), C.p_a);
                T(:,:,2,ii) = tr4(eye(3), C.p_s);
                T(:,:,3,ii) = tr4(Rax(C.omg_m, thM(j)), C.p_m);
            end
        end
    case 'unihip'
        thX = linspace(-30, 20.6, D) * pi/180;
        thZ = linspace(-15, 110, D) * pi/180;
        T = zeros(4,4,1,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax([1;0;0], thX(j)) * Rax([0;0;1], thZ(i)), C.p_h);
            end
        end
    case 'unihip_psoas'
        thX = linspace(-30, 20.6, D) * pi/180;
        thZ = linspace(-15, 110, D) * pi/180;
        T = zeros(4,4,2,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(eye(3), C.p_b);
                T(:,:,2,ii) = tr4(Rax([1;0;0], thX(j)) * Rax([0;0;1], thZ(i)), C.p_h);
            end
        end
    otherwise
        error('bad sizing case %s', caseName);
end
end

function [Loc, thH, thK] = rebuildRectFem(mo, D, C)
% rect_fem 3-D Location over the bihip grid (connor_Bi_Hip.m:115-124).
thH = linspace(-15, 110, D) * pi/180;
thK = linspace(-2.0943951, 0.17453293, D);
Loc = zeros(3, 5, D*D); ii = 0;
for i = 1:D
    for j = 1:D
        ii = ii + 1;
        Loc(:, :, ii) = [mo.Location(:, 1:4, 1), [C.fcn3(thK(j)); C.fcn4(thK(j)); mo.Location(3, 5, 1)]];
    end
end
end

function [Loc, thH, thK] = rebuildVasInt(mo, D, C)
% vas_int 3-D Location over the bihip grid (connor_Bi_Hip.m:154-163).
thH = linspace(-15, 110, D) * pi/180;
thK = linspace(-2.0943951, 0.17453293, D);
Loc = zeros(3, 4, D*D); ii = 0;
for i = 1:D
    for j = 1:D
        ii = ii + 1;
        Loc(:, :, ii) = [mo.Location(:, 1:3, 1), [C.fcn7(thK(j)); C.fcn8(thK(j)); mo.Location(3, 4, 1)]];
    end
end
end

function Loc = rebuildRectFemG(mo, thK, C)
n = numel(thK);
Loc = zeros(3, 5, n);
for i = 1:n
    Loc(:, :, i) = [mo.Location(:, 1:4, 1), [C.fcn3(thK(i)); C.fcn4(thK(i)); mo.Location(3, 5, 1)]];
end
end

function Loc = rebuildVasIntG(mo, thK, C)
n = numel(thK);
Loc = zeros(3, 4, n);
for i = 1:n
    Loc(:, :, i) = [mo.Location(:, 1:3, 1), [C.fcn7(thK(i)); C.fcn8(thK(i)); mo.Location(3, 4, 1)]];
end
end

function E = addSpec(E, id, section, idx, tcase, row, col, s, ben, mifname, pri_dof, sec_dof)
s0 = struct('id', id, 'section', section, 'idx', idx, 'tcase', tcase, ...
    'row', row, 'col', col, 's', s, 'ben', {ben}, 'mifname', mifname, ...
    'pri_dof', pri_dof, 'sec_dof', sec_dof);
if isempty(E)
    E = s0;
else
    E(end+1) = s0; %#ok<AGROW>
end
end

function d = maxAbs(A, B)
d = max(abs(A(:) - B(:)));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function s = formatAgree(x)
if isnan(x), s = 'n/a'; else, s = sprintf('%.0f%%', x); end
end
