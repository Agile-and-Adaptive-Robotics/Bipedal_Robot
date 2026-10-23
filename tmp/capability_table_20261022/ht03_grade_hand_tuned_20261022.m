% Hand-tuned BPA capability row (2026-10-22), graded on eb475ws4.
% Evaluates Ben's HAND-TUNED BPA designs (the old-code robot routing) against
% the human torque curves CSV, using the metric definitions of the optimized
% row (Opt_run_BiPulley_survey_20261022.m grade(): meet/active, 25-cell grid).
%
% SOURCES (read-only):
%   - Connor's corrected kinematics: repo Code/Matlab/Previous Optimization
%     Code/Joint_Optimization_Scripts/RobotPAMCalculationOptimization_*.m,
%     run VERBATIM here for the BPA sizing (their own 100x100 grids).
%   - Ben's originals: tmp/old_code_copy_20261022 (AttachPoints_Robot.mat,
%     MovePoints.mat, PAM_Prop.mat, Size.m) for the value comparison only.
%   - Human targets: tmp/capability_table_20261022/human_torque_curves_20261022.csv
%     (copied fresh from easteregg2 immediately before this run).
%   - Force model: Code/Matlab/Functions/maxBPAforce.m at 620 kPa.
%
% MECHANICS: the repo's PamData.m FAILS TO PARSE (line 271 '9hi' corruption,
% verified 2026-10-22). This driver uses the staged copy PamData_20261022.m
% (verbatim except class rename + that one token removed; Ben's original
% Size.m carries the clean line) and, to run Connor's section scripts
% verbatim, a staged shadow class named PamData.m in THIS folder (same
% content). The repo and old-code files are never modified.
%
% CONVENTIONS:
%   - Old-code joint signs vs OpenSim, fixed per DOF by an anatomical
%     cross-check (old-code MA signs vs each muscle's action):
%     knee/ankle/hip-flexion/lumbar-extension/mtp positive = extension/
%     dorsiflexion/flexion/extension/extension (s = +1, matches Ben's
%     comments AND the muscle-action check). For the SUBTALAR and HIP-X
%     axes Ben's robot-script comments ("inversion-/eversion+",
%     "adduction-/abduction+") are INVERTED vs the actual axis geometry:
%     the old-code tibialis-posterior MA is POSITIVE and the peronei MAs
%     are NEGATIVE (tib_post inverts, peronei evert), and add_mag's MA is
%     POSITIVE (it adducts) -- so the robot's +subtalar = INVERSION and
%     +hip-x = ADDUCTION, same as OpenSim (s = +1 for all DOFs). The
%     human CSV confirms: tib_post tau_h > 0, peronei tau_h < 0,
%     adductor tau_h > 0.
%   - Strain = (RestingL - (Lmt - tendon))/RestingL with tendon =
%     LongestL - RestingL (old festo() convention; same contraction-strain
%     structure as the optimized row's strain_p).
%   - meet = taut & (sign(tauB)==sign(tauH) | |tauH|<=1e-9) & |tauB|>=|tauH|-1e-12
%     with taut = strain>0 & isfinite(tauB);  active = taut & |tauB|>0.
%
% OUTPUT: hand_tuned_metrics_20261022.csv (muscle,source,meetPct,activePct,
% suspect,note) in this folder + a full stdout report.

FAST = ~isempty(getenv('HT_FAST'));   % HT_FAST=1: Phase A grids at divisions=25 (validation only).

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

%% 1. Map + human CSV ---------------------------------------------------------
mapPath = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Routing27_v2', 'actuator_map_27.json');
mj = jsondecode(fileread(mapPath));
acts = mj.actuators;
assert(numel(acts) == 27, 'map must have 27 actuators');

csvPath = fullfile(here, 'human_torque_curves_20261022.csv');
raw = readtable(csvPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char', 'CommentStyle', '#');
vars = raw.Properties.VariableNames;
iAng = find(startsWith(vars, 'theta_deg_'));
iTau = find(startsWith(vars, 'tau_h_Nm_'));
assert(numel(iAng) == 25 && numel(iTau) == 25, 'CSV must carry 25+25 grid columns');
angAll = table2array(raw(:, iAng));      % 27 x 25 (degrees, OpenSim)
tauAll = table2array(raw(:, iTau));      % 27 x 25 (N*m, OpenSim native sign)
aid = raw.actuator_id;
fprintf('Human CSV: %d rows, 25-pt grids.\n', height(raw));
% verify grid vs map rom per actuator (CSV prints 6 decimals -> 1e-4 tol)
gridOK = true;
for k = 1:27
    r = mj.rom_table_deg.(acts(k).primary_dof);
    a = angAll(aid == acts(k).id, :);
    ideal = r(1) + (r(2) - r(1)) * (0:24) / 24;
    if max(abs(a - ideal)) > 1e-4
        gridOK = false;
        fprintf('GRID MISMATCH act %d %s: csv [%g..%g] vs map [%g..%g] maxdev %.2g\n', ...
            acts(k).id, acts(k).primary_dof, a(1), a(25), r(1), r(2), max(abs(a - ideal)));
    end
    if ~strcmp(raw.primary_dof{aid == acts(k).id}, acts(k).primary_dof)
        gridOK = false;
        fprintf('DOF MISMATCH act %d\n', acts(k).id);
    end
end
fprintf('Grid/dof verification vs map rom_table_deg: %s\n', string(gridOK));

%% 2. PHASE A: run Connor's six section scripts verbatim ----------------------
% Sizing (Diameter/RestingL/LongestL) and the routing VALUES come straight
% from the scripts' own Muscles objects (zero transcription).
sec = struct();
divSet = 100; if FAST, divSet = 25; end

clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts raw vars iAng iTau angAll tauAll aid csvPath mapPath host
divisions = divSet;   % pre-set: scripts honor an existing 'divisions'
run(fullfile(JOS, 'RobotPAMCalculationOptimization_Back.m'));
sec.back = grabMuscles(Muscles); clear Muscles

clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts raw vars iAng iTau angAll tauAll aid csvPath mapPath host
divisions = divSet;
run(fullfile(JOS, 'RobotPAMCalculationOptimization_Bi_Hip.m'));
sec.bihip = grabMuscles(Muscles); clear Muscles

clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts raw vars iAng iTau angAll tauAll aid csvPath mapPath host
divisions = divSet;
run(fullfile(JOS, 'RobotPAMCalculationOptimization_Calves.m'));
sec.calves = grabMuscles(Muscles); clear Muscles   % {MG, LG, Soleus}

clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts raw vars iAng iTau angAll tauAll aid csvPath mapPath host
divisions = divSet;
run(fullfile(JOS, 'RobotPAMCalculationOptimization_Foot.m'));
sec.foot = grabMuscles(Muscles); clear Muscles     % {tib_post, tib_ant, per_brev, per_long, per_tert}

clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts raw vars iAng iTau angAll tauAll aid csvPath mapPath host
divisions = divSet;
run(fullfile(JOS, 'RobotPAMCalculationOptimization_Toe.m'));
sec.toe = grabMuscles(Muscles); clear Muscles      % {FDL, FHL, EDL, EHL}

% Uni_Hip: the script ERRORS at Muscle4 (PLocation undefined; map documents
% 'presumed Location4'). try/catch, then build Psoas from the script's own
% workspace values (Location4 etc. as left at the error point).
clearvars -except root here oldc repoPOC JOS sec FAST divSet mj acts raw vars iAng iTau angAll tauAll aid csvPath mapPath host
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

fprintf('\n=== PHASE A: section muscles (Dia/Rest/Long from Connor''s scripts) ===\n');
fn = fieldnames(sec);
for f = 1:numel(fn)
    S = sec.(fn{f});
    fprintf('%-8s:', fn{f});
    for m = 1:numel(S)
        fprintf('  %s[%.0fmm r=%.4g L=%.4g]', S(m).Muscle, S(m).Diameter, S(m).RestingL, S(m).LongestL);
    end
    fprintf('\n');
end

%% 3. Constants + builders ----------------------------------------------------
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

%% 4. PHASE C: verify my T builders + 3-D Location rebuilds -------------------
fprintf('\n=== PHASE C: verification vs the scripts'' own T / Location ===\n');
D = divSet;
ver = struct('name', {}, 'diff', {});
Tb = buildSzT('back', D, C); Tb(:,:,2,:) = Tb(:,:,1,:); Tb(:,:,3,:) = Tb(:,:,1,:);  % connor_Back.m:55-56 duplication
ver(end+1) = struct('name', 'back T',       'diff', maxAbs(Tb, sec.back(1).TransformationMat));
ver(end+1) = struct('name', 'bihip T',      'diff', maxAbs(buildSzT('bihip', D, C), sec.bihip(1).TransformationMat));
% connor_Calves allocates a 3-slot T and fills slots 1-2 for the soleus
% (slot 3 zeros, inert: soleus CrossPoints [2 2] never touches slot 3);
% compare the used slots only.
Tcs = buildSzT('calves_soleus', D, C);
ver(end+1) = struct('name', 'calves soleus T (slots 1-2)', 'diff', ...
    maxAbs(Tcs, sec.calves(3).TransformationMat(:, :, 1:2, :)));
ver(end+1) = struct('name', 'calves gastroc T', 'diff', maxAbs(buildSzT('calves_gastroc', D, C), sec.calves(1).TransformationMat));
% connor_Foot also allocates a 3-slot T, filling slots 1-2 (slot 3 zeros,
% inert for CrossPoints [3 3]/[4 4]); compare the used slots only.
Tft = buildSzT('foot', D, C);
ver(end+1) = struct('name', 'foot T (slots 1-2)', 'diff', ...
    maxAbs(Tft, sec.foot(1).TransformationMat(:, :, 1:2, :)));
ver(end+1) = struct('name', 'toe T',        'diff', maxAbs(buildSzT('toe', D, C), sec.toe(1).TransformationMat));
ver(end+1) = struct('name', 'unihip T',     'diff', maxAbs(buildSzT('unihip', D, C), sec.unihip(1).TransformationMat));
ver(end+1) = struct('name', 'unihip psoas T', 'diff', maxAbs(buildSzT('unihip_psoas', D, C), sec.unihip(4).TransformationMat));
% 3-D Location rebuilds at the bihip grid
[LocR, ~, ~] = rebuildRectFem(sec.bihip(5), D, C);
ver(end+1) = struct('name', 'rect_fem Location rebuild', 'diff', maxAbs(LocR, sec.bihip(5).Location));
[LocV, ~, ~] = rebuildVasInt(sec.bihip(8), D, C);
ver(end+1) = struct('name', 'vas_int Location rebuild', 'diff', maxAbs(LocV, sec.bihip(8).Location));
for v = 1:numel(ver)
    fprintf('  %-32s max|diff| = %.3g %s\n', ver(v).name, ver(v).diff, ...
        ternary(ver(v).diff <= 1e-12, 'EXACT (machine noise)', '*** DIFFERS ***'));
end
assert(all([ver.diff] <= 1e-12), 'builder verification failed');

%% 5. PHASE B: grading pass on the human CSV grids ----------------------------
% Per-actuator spec: section, muscle index in the section, grading T case,
% torque (row, axisCol), sign s (OpenSim = s * old-code), Ben comparison points.
APR = load(fullfile(oldc, 'AttachPoints_Robot.mat'));
PP  = load(fullfile(oldc, 'PAM_Prop.mat'));
E = [];
E = addSpec(E, 1,'back',1,'back',1,3,+1,{'p7','ercspn_wr1','ercspn_wr2'}, 'ercspn');
E = addSpec(E, 2,'back',2,'back',1,3,+1,{'p2','b3'}, 'intobl');
E = addSpec(E, 3,'back',3,'back',1,3,+1,{'p4','b4'}, 'extobl');
E = addSpec(E, 4,'bihip',1,'knee',2,3,+1,{'p5','twr2','twr3','t5'}, 'bifemlh');
E = addSpec(E, 5,'bihip',2,'knee',2,3,+1,{'p5','semimem_wr1','semimem_wr2','twr1','t1','t6'}, 'semimem');
E = addSpec(E, 6,'bihip',3,'knee',2,3,+1,{'p3','sar_wr1','sar_wr2','sar_wr3','twr1','t1','t6'}, 'sar');
E = addSpec(E, 7,'bihip',4,'knee',2,3,+1,{'p3','tfl_wr1','tfl_wr2','t2','t5'}, 'tfl');
E = addSpec(E, 8,'bihip',5,'knee3dR',2,3,+1,{'p10','rect_fem_wr1','rect_fem_wr2','rect_fem_wr3'}, 'rect_fem');
E = addSpec(E, 9,'bihip',6,'knee',2,3,+1,{'p8','grac_wr1','grac_wr2','twr1','t1','t6'}, 'grac');
E = addSpec(E,10,'bihip',7,'knee',2,3,+1,{'f2','twr2','twr3'}, 'bifemsh');
E = addSpec(E,11,'bihip',8,'knee3dV',2,3,+1,{'f1','vas_int_wr1','vas_int_wr2'}, 'vas_int');
E = addSpec(E,12,'calves',3,'ankle2',1,1,+1,{'t1','a1'}, 'soleus');
E = addSpec(E,13,'calves',1,'ankle3',2,1,+1,{'f4','med_gas_wr1','a1'}, 'med_gas');
E = addSpec(E,14,'calves',2,'ankle3',2,1,+1,{'f5','lat_gas_wr1','a1'}, 'lat_gas');
E = addSpec(E,15,'foot',1,'subtal2',2,1,+1,{'t4','tib_post_wr1','tib_post_wr2','a4'}, 'tib_post');
E = addSpec(E,16,'foot',2,'ankle2',1,1,+1,{'t2','tib_ant_wr1','a3'}, 'tib_ant');
E = addSpec(E,17,'foot',3,'subtal2',2,1,+1,{'t2','per_brev_wr1','per_brev_wr2','per_brev_wr3','a2'}, 'per_brev');
E = addSpec(E,18,'foot',4,'subtal2',2,1,+1,{'t2','per_long_wr1','per_long_wr2','per_long_wr3','per_long_wr4','per_long_wr5','a3'}, 'per_long');
E = addSpec(E,19,'foot',5,'subtal2',2,1,+1,{'t2','per_tert_wr1','a2'}, 'per_tert');
E = addSpec(E,20,'toe',1,'mtp3',3,1,+1,{'t4','flex_dig_wr1','flex_dig_wr2','flex_dig_wr3','flex_dig_wr4','flex_dig_wr5','flex_dig_wr6','a7'}, 'flex_dig');
E = addSpec(E,21,'toe',2,'mtp3',3,1,+1,{'t4','flex_hal_wr1','flex_hal_wr2','flex_hal_wr3','flex_hal_wr4','flex_hal_wr5','a8'}, 'flex_hal');
E = addSpec(E,22,'toe',3,'mtp3',3,1,+1,{'t2','ext_dig_wr1','ext_dig_wr2','ext_dig_wr3','ext_dig_wr4','a5'}, 'ext_dig');
E = addSpec(E,23,'toe',4,'mtp3',3,1,+1,{'t2','ext_hal_wr1','ext_hal_wr2','ext_hal_wr3','ext_hal_wr4','ext_hal_wr5','a6'}, 'ext_hal');
E = addSpec(E,24,'unihip',1,'hipflex1',1,3,+1,{'p6','glut_max_wr1','glut_max_wr2','glut_max_wr3','f3'}, 'glut_max');
E = addSpec(E,25,'unihip',2,'hipadd1',1,1,+1,{'p9','add_mag_wr1','add_mag_wr2','f4'}, 'add_mag');
E = addSpec(E,26,'unihip',3,'hipflex1',1,3,+1,{'p1','iliacus_wr1','f4'}, 'iliacus');
E = addSpec(E,27,'unihip',4,'hipflex2',2,3,+1,{'b1','psoas_wr1','f4'}, 'psoas');
assert(numel(E) == 27);

fprintf('\n=== PHASE B: grading on the human CSV 25-pt grids (620 kPa) ===\n');
res = [];
for k = 1:27
    a = acts(k); e = E(k);
    mo = sec.(e.section)(e.idx);   % muscle object values from Connor's script
    assert(isnumeric(mo.RestingL) && isnumeric(mo.Diameter), ...
        'act %d: non-numeric BPA size from Size()', a.id);
    q  = angAll(aid == a.id, :);            % 1x25 OpenSim degrees
    tauH = tauAll(aid == a.id, :);          % 1x25 N*m OpenSim
    thOld = e.s * q;                        % old-code joint angle (degrees)
    % grading T
    Tg = buildGradedT(e.tcase, thOld * d2r, C);
    % Location (object values as-is; 3-D rebuilds for rect_fem/vas_int)
    switch e.tcase
        case 'knee3dR', Loc = rebuildRectFemG(mo, thOld * d2r, C);
        case 'knee3dV', Loc = rebuildVasIntG(mo, thOld * d2r, C);
        otherwise,      Loc = mo.Location;
    end
    mG = PamData_20261022(mo.Muscle, Loc, mo.CrossPoints, mo.MIF, Tg, mo.Axis);
    Lmt = reshape(mG.MuscleLength, 1, []);
    MA  = reshape(mG.MomentArm(e.row, :, e.col), 1, []);
    F   = maxBPAforce(mo.RestingL, num2str(mo.Diameter), 620);
    strain = (mo.LongestL - Lmt) / mo.RestingL;      % festo contraction strain
    tauB = e.s * (MA * F);                            % OpenSim-sign torque
    taut  = strain > 0 & isfinite(tauB);
    meet  = taut & isfinite(tauH) & (sign(tauB) == sign(tauH) | abs(tauH) <= 1e-9) ...
            & abs(tauB) >= abs(tauH) - 1e-12;
    activ = taut & isfinite(strain) & strain > 0 & abs(tauB) > 0;
    meetPct = 100 * sum(meet) / 25;  actPct = 100 * sum(activ) / 25;
    mSk = taut & abs(tauH) > 1;
    if any(mSk), signAgree = 100 * mean(sign(tauB(mSk)) == sign(tauH(mSk))); else, signAgree = NaN; end
    r = struct('id', a.id, 'name', a.master_name, 'dof', a.primary_dof, ...
        'dia', mo.Diameter, 'rest', mo.RestingL, 'long', mo.LongestL, 'F', F, ...
        'strainMin', min(strain), 'strainMax', max(strain), ...
        'tauMin', min(tauB), 'tauMax', max(tauB), 'maMin', min(MA), 'maMax', max(MA), ...
        'meetPct', meetPct, 'actPct', actPct, 'signAgree', signAgree);
    if k == 1, res = r; else, res(k) = r; end %#ok<AGROW>
    fprintf(['%2d %-32s %-17s d%2.0fmm r=%.4g F=%7.1fN strain[%+.3f,%+.3f] ' ...
        'tau[%+8.2f,%+8.2f] meet=%5.1f%% active=%5.1f%% signAgree=%s\n'], ...
        r.id, r.name, r.dof, r.dia, r.rest, r.F, r.strainMin, r.strainMax, ...
        r.tauMin, r.tauMax, meetPct, actPct, formatAgree(signAgree));
end

%% 6. PHASE D: Connor-vs-Ben value comparison + MIF sanity --------------------
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
% duplicate-column detection (add_mag case)
amLoc = sec.unihip(2).Location;
fprintf('add_mag Location columns 2 and 3 identical: %s\n', ...
    string(isequal(amLoc(:,2), amLoc(:,3))));
fprintf('ES: Ben originals carry 9 points (p7, ercspn_wr1..7, b2); Connor keeps 3.\n');
fprintf('ES: Ben hand-set ercspn_rest = 0.7428 (bipedal_kinematics_back.m:311); Connor flow uses Size output %.6g.\n', sec.back(1).RestingL);

%% 7. Output CSV + report -----------------------------------------------------
if FAST
    outPath = fullfile(here, 'hand_tuned_metrics_20261022_FAST_validation.csv');
else
    outPath = fullfile(here, 'hand_tuned_metrics_20261022.csv');
end
fid = fopen(outPath, 'w');
fprintf(fid, 'muscle,source,meetPct,activePct,suspect,note\n');
for k = 1:27
    r = res(k); mm = cmp{k}.mm;
    suspect = 'no'; note = '';
    switch r.id
        case 1
            note = 'Connor routing keeps 3 of Ben''s 9 points (p7 wr1 wr2; Ben full: p7 wr1-wr7 b2); retained values match Ben; Ben hand-set rest 0.7428 m not in Connor flow (Size output used)';
        case 5
            note = 'Connor applies knee transform to twr1 where Ben''s inline formula applies hip (Connor corrected kinematics); values match Ben';
        case 8
            note = 'static points match Ben; patella point via fcn3/fcn4 identical to Ben''s MovePoints.mat fcns';
        case 9
            note = 'origin differs from Ben p8; Connor y z equal rect_fem origin p10''s (likely transcription slip)';
        case 11
            note = 'static points match Ben; patella point via fcn7/fcn8 identical to Ben''s MovePoints.mat fcns';
        case 15
            suspect = 'yes';
            note = 'MIF 3549 = map master_mif_n but Ben PAM_Prop has tib_post_mif 1588 (map documents 3549 as soleus copy-paste artifact; same 40 mm diameter either way); routing values match Ben';
        case 16
            note = 'Connor names the muscle Tibialis Posterior (duplicate-name bug documented in map; MIF 905 and route identify tib_ant); values match Ben';
        case 19
            note = 'three components differ from Ben (values appear shuffled between t2/per_tert_wr1/a2)';
        case 21
            note = 'Connor passes Axis1 instead of Axis2 (both [3 1 3]; numerically identical); values match Ben';
        case 25
            suspect = 'yes';
            note = 'Location columns 2 and 3 both equal add_mag_wr2; Ben routing uses add_mag_wr1 then add_mag_wr2 (wr1 absent; likely transcription slip)';
        case 27
            suspect = 'yes';
            note = 'PLocation undefined-variable slip (values = Location4; map: presumed Location4); b1/wr1 frame roles swapped vs Ben''s inline formula; point values match Ben';
    end
    % append the machine comparison facts to every note
    if isempty(mm)
        if ~contains(note, 'match Ben'), note = [note '; values match Ben']; end
    else
        suspect = 'yes';
        parts = '';
        for j = 1:numel(mm)
            parts = [parts sprintf(' col%d vs Ben %s: Connor [%s] Ben [%s];', ...
                mm(j).col, mm(j).ben, mat2str(mm(j).connor', 4), mat2str(mm(j).ben_v', 4))]; %#ok<AGROW>
        end
        note = [note '; VALUE MISMATCH:' parts];
    end
    note = strtrim(note);
    while startsWith(note, ';'), note = strtrim(note(2:end)); end
    while endsWith(note, ';'), note = strtrim(note(1:end-1)); end
    fprintf(fid, '"%s",connor-corrected,%.4f,%.4f,%s,"%s"\n', ...
        strrep(r.name, '"', ''''), r.meetPct, r.actPct, suspect, strrep(note, '"', ''''));
    res(k).suspect = suspect; res(k).note = note; %#ok<AGROW>
end
fclose(fid);
fprintf('\nWROTE %s\n', outPath);
nSus = sum(strcmp({res.suspect}, 'yes'));
fprintf('Coverage: 27/27 connor-corrected, 0 ben-original, 0 excluded (Connor''s six section scripts define all 27 master actuators).\n');
fprintf('Suspects flagged: %d\n', nSus);
fprintf('Mean meetPct (covered): %.2f | mean activePct: %.2f\n', ...
    mean([res.meetPct]), mean([res.actPct]));
fprintf('GRADE DONE\n');

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
% Grading T: primary joint swept over th (rad, OLD-code sign); all other
% joints at their default (0 deg, identity rotation + home translation).
n = numel(th);
switch caseName
    case 'back'        % lumbar extension (z), 1 slot
        T = zeros(4,4,1,n);
        for i = 1:n, T(:,:,1,i) = tr4(Rax([0;0;1], th(i)), C.p_b); end
    case {'knee', 'knee3dR', 'knee3dV'}  % knee swept (2-D or 3-D Location), hip 0; 2 slots
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(eye(3), C.p_h);
            T(:,:,2,i) = tr4(Rax([0;0;1], th(i)), kneeHomeOf(th(i), C));
        end
    case 'ankle2'      % ankle swept (soleus / tib_ant); 2 slots
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(Rax(C.omg_a, th(i)), C.p_a);
            T(:,:,2,i) = tr4(eye(3), C.p_s);
        end
    case 'ankle3'      % ankle swept (gastrocs); knee 0, subtalar 0; 3 slots
        T = zeros(4,4,3,n);
        kh0 = kneeHomeOf(0, C);
        for i = 1:n
            T(:,:,1,i) = tr4(eye(3), kh0);
            T(:,:,2,i) = tr4(Rax(C.omg_a, th(i)), C.p_a);
            T(:,:,3,i) = tr4(eye(3), C.p_s);
        end
    case 'subtal2'     % subtalar swept (foot muscles); ankle 0; 2 slots
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(eye(3), C.p_a);
            T(:,:,2,i) = tr4(Rax(C.omg_s, th(i)), C.p_s);
        end
    case 'mtp3'        % mtp swept (toe muscles); ankle 0, subtalar 0; 3 slots
        T = zeros(4,4,3,n);
        for i = 1:n
            T(:,:,1,i) = tr4(eye(3), C.p_a);
            T(:,:,2,i) = tr4(eye(3), C.p_s);
            T(:,:,3,i) = tr4(Rax(C.omg_m, th(i)), C.p_m);
        end
    case 'hipflex1'    % hip flexion (z) swept; 1 slot (glut_max, iliacus)
        T = zeros(4,4,1,n);
        for i = 1:n, T(:,:,1,i) = tr4(Rax([0;0;1], th(i)), C.p_h); end
    case 'hipadd1'     % hip adduction axis (x, old + = abduction) swept; 1 slot
        T = zeros(4,4,1,n);
        for i = 1:n, T(:,:,1,i) = tr4(Rax([1;0;0], th(i)), C.p_h); end
    case 'hipflex2'    % psoas: back 0 (slot 1), hip flexion swept (slot 2)
        T = zeros(4,4,2,n);
        for i = 1:n
            T(:,:,1,i) = tr4(eye(3), C.p_b);
            T(:,:,2,i) = tr4(Rax([0;0;1], th(i)), C.p_h);
        end
    otherwise
        error('bad graded case %s', caseName);
end
end

function T = buildSzT(caseName, D, C)
% Connor's own section grids (divisions = D), pose order ii = (i-1)*D + j.
switch caseName
    case 'back'   % BackX(i) x BackZ(j) (connor_Back.m:44-56)
        thX = linspace(-15, 15, D) * pi/180;
        thZ = linspace(-45, 15, D) * pi/180;
        T = zeros(4,4,1,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax([1;0;0], thX(i)) * Rax([0;0;1], thZ(j)), C.p_b);
            end
        end
    case 'bihip'  % hip z(i) x knee(j) (connor_Bi_Hip.m:51-59)
        thH = linspace(-15, 110, D) * pi/180;
        thK = linspace(-2.0943951, 0.17453293, D);  % Connor's exact knee constants (connor_Bi_Hip.m:34-35)
        T = zeros(4,4,2,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax([0;0;1], thH(i)), C.p_h);
                T(:,:,2,ii) = tr4(Rax([0;0;1], thK(j)), kneeHomeOf(thK(j), C));
            end
        end
    case 'calves_soleus' % ankle(i), subtalar const (connor_Calves.m:67-75)
        thA = linspace(-50, 20, D) * pi/180;
        T = zeros(4,4,2,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax(C.omg_a, thA(i)), C.p_a);
                T(:,:,2,ii) = tr4(eye(3), C.p_s);
            end
        end
    case 'calves_gastroc' % knee(j), ankle(i), subtalar const (connor_Calves.m:91-98)
        thA = linspace(-50, 20, D) * pi/180;
        thK = linspace(-2.0943951, 0.17453293, D);  % Connor's exact knee constants (connor_Bi_Hip.m:34-35)
        T = zeros(4,4,3,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax([0;0;1], thK(j)), kneeHomeOf(thK(j), C));
                T(:,:,2,ii) = tr4(Rax(C.omg_a, thA(i)), C.p_a);
                T(:,:,3,ii) = tr4(eye(3), C.p_s);
            end
        end
    case 'foot'   % ankle(i) x subtalar(j) (connor_Foot/foot_example.m:43-51)
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
    case 'toe'    % ankle(i), subtalar const, mtp(j) (connor_Toe.m:50-59)
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
    case 'unihip' % HipZ(i) x HipX(j): R = Rx(thX_j)*Rz(thZ_i) (connor_Uni_Hip.m:44-53)
        thX = linspace(-30, 20.6, D) * pi/180;
        thZ = linspace(-15, 110, D) * pi/180;
        T = zeros(4,4,1,D*D); ii = 0;
        for i = 1:D
            for j = 1:D
                ii = ii + 1;
                T(:,:,1,ii) = tr4(Rax([1;0;0], thX(j)) * Rax([0;0;1], thZ(i)), C.p_h);
            end
        end
    case 'unihip_psoas' % back const (slot 1), hip (slot 2) (connor_Uni_Hip.m:89-99)
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
thK = linspace(-2.0943951, 0.17453293, D);  % Connor's exact knee constants (connor_Bi_Hip.m:34-35)
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
thK = linspace(-2.0943951, 0.17453293, D);  % Connor's exact knee constants (connor_Bi_Hip.m:34-35)
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

function E = addSpec(E, id, section, idx, tcase, row, col, s, ben, mifname)
s0 = struct('id', id, 'section', section, 'idx', idx, 'tcase', tcase, ...
    'row', row, 'col', col, 's', s, 'ben', {ben}, 'mifname', mifname);
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
