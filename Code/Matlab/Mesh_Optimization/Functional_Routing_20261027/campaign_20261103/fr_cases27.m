function CASES = fr_cases27(C, kin, overrides, opts)
%FR_CASES27 Build the 27 contract-driven campaign cases (run_0; v2 adds
% the optional per-group sizing overrides, see fr_applyCfg27).
%
% One case per target_contract.json group. Every case is contract-aligned:
% PASS 1 sweeps the contract's primaryDof over the actuator-map ROM table
% (knee passes use C.kneeRangeDeg [-120 0] deg - master item 6: no
% hyperextension - which overrides the stock-model table's +10 deg), so
% fr_targets and fr_objective score the contract's primary dof. A second
% pass adds the route's other crossed joint as constraint coverage (all
% hard constraints are evaluated over ALL passes; only pass 1 is scored).
%
% Row mechanics sets (documented per case in cs.note):
%   - ConditionalPathPoints are DROPPED from the mechanics row set
%     everywhere except the gas groups (med_gas/lat_gas keep their
%     condyle-wrap row FIXED, the gate-verified knee_bowden structure);
%     the declared conditional behavior stays implemented in fr_evalRoute
%     and is tracked per pose (stage-1 convention, gate check (c)).
%   - MovingPathPoints (rect_fem_r-P3, vas_lat_r-P5) are KEPT as
%     FIXED-TRACKED mechanics rows: seed location = fr_evalRoute's
%     declared location at the pass-1 MID pose, owning-body frame (the
%     stage-1 fixed-row convention with tracked behavior; README
%     limitation 2). They are never movable attachments.
%   - Sleeve = rows 1..cross(1)-1. Short routes (2 declared rows: trunk
%     x3, soleus, add_mag3, bifemsh, psoas-family hip cases with nProx=1)
%     carry the same zero-length-sleeve-in-grid convention as the
%     gate-passing wired trunk/hip cases (hardware at brackets; the bone
%     obstacle is vacuous over zero segments).
%
% Tackle policy (recorded per case):
%   wired-feasibility-exemplar overrides (gate-passed configs):
%     psoas/add_mag3/ercspn/med_gas -> G=2, family L.
%   ALL others: mechanical rule - G=2 family L iff the declared route's
%     pass-1 length range dML exceeds the single-BPA contraction band
%     kmaxF * restEst (restEst = boundary-taut seed estimate). This
%     REPLACES the wired sar demonstrator's no-tackle config, which was
%     penalty-dominated (f 10135 in the gate layering proof), with a
%     mechanically feasible one; documented here and in the manifest.
%
% nBPA = contract sizing n20_ceil_MIF_over_Fmax20; diameter '10' iff the
% contract's 10 mm screen passed, else '20'. The envelope/soft-overlap
% branch ring caps at 3 offsets (fr_bundleOffsets/fr_ringOffsets), so
% 4-branch groups (glut_max3, vas_int) under-represent one branch in the
% envelope - recorded machinery cap, later campaign item.
%
% Transform convention (gate-verified fr_chainT direction): T1 =
% fr_chainT(kin, q, childFrame, proxFrame) maps distal-row coords into
% the proximal (sleeve-side) frame; constant expr matrices re-express
% rows whose assigned frame differs from their owning body (built at the
% default pose, exact for the coordinates held constant in that pass).

tc = jsondecode(fileread(fullfile(C.targetsDir, 'target_contract.json')));
am = jsondecode(fileread(fullfile(C.repo, 'Code', 'Matlab', ...
    'Mesh_Optimization', 'Routing27_v2', 'actuator_map_27.json')));
% optional per-group sizing overrides (fr_sizingSearch winners): containers.Map
% group name -> struct('diaStr', '10'|'20', 'nBPA', n, 'family', 'L'|'M'|'S', ...
% 'gain', g). Absent/empty = the contract seeds (Ben 2026-11-04 item 4: the
% contract counts are search centers, never mandatory).
if nargin < 3
    overrides = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
if nargin < 4
    opts = struct();
end

caseCell = {};
for g = 1:numel(tc.groups)
    grp = tc.groups(g);
    cs = fr_buildOne27(C, kin, grp, am, overrides, opts);
    caseCell{end + 1} = cs;
end
CASES = caseCell{1};
for i = 2:numel(caseCell)
    CASES(i) = caseCell{i};
end
end

% ==========================================================================
function cs = fr_buildOne27(C, kin, grp, am, overrides, opts)
% Dispatch on the contract group id (explicit, deterministic). Each
% pattern builder derives its row/cross structure from the ACTUAL parsed
% route (never hardcoded row indices).

[~, im] = max(structfun(@(x) x, grp.sizing_data.per_muscle_Lmax_m));
seedMuscle = grp.stockMuscles{im};
% route source of record: the stock human model by default (optimization
% seeds); opts.routeFromRobotbody = true switches the seed geometry to
% gait2392_robotbody.osim's OWN path points -- the 2019-vintage robot
% routings, for the BASELINE comparison (Ben 2026-11-04). Connor-era
% frozen conditional artifacts in robotbody are handled by
% fr_activeRows27 (declared-range rows dropped from the mechanics row
% set, behavior tracked in fr_evalRoute), same as for stock routes.
if nargin < 6
    opts = struct();
end
routeSrc = C.osimStock;
if isfield(opts, 'routeFromRobotbody') && opts.routeFromRobotbody
    routeSrc = C.osimRobot;
end
route = fr_muscleroute(routeSrc, seedMuscle);
% case name of record = the contract group key (checkpoint filenames use it)
grp.name = grp.group;
diaStr = fr_diaStrOf(grp);
nBPA = grp.sizing_data.n20_ceil_MIF_over_Fmax20;

switch grp.map_id
    case {24, 25, 26, 27}   % glut_max3, add_mag3, iliacus, psoas
        cs = pbHip(C, kin, grp, route, seedMuscle, diaStr, nBPA, am);
    case {1, 2, 3}          % ercspn, intobl, extobl
        cs = pbTrunk(C, kin, grp, route, seedMuscle, diaStr, nBPA, am);
    case {4, 5, 6, 7, 8, 9} % bifemlh, semimem, sar, tfl, rect_fem, grac
        cs = pbBi(C, kin, grp, route, seedMuscle, diaStr, nBPA, am);
    case {10, 11}           % bifemsh, vas_int
        cs = pbKneeMono(C, kin, grp, route, seedMuscle, diaStr, nBPA, am);
    case {13, 14}           % med_gas, lat_gas
        cs = pbGas(C, kin, grp, route, seedMuscle, diaStr, nBPA, am);
    case {12, 16}           % soleus, tib_ant
        cs = pbAnkleMono(C, kin, grp, route, seedMuscle, diaStr, nBPA, am);
    case {15, 17, 18, 19}   % tib_post, per_brev, per_long, per_tert
        cs = pbSubtalar(C, kin, grp, route, seedMuscle, diaStr, nBPA, am);
    case {20, 21, 22, 23}   % flex_dig, flex_hal, ext_dig, ext_hal
        cs = pbFoot(C, kin, grp, route, seedMuscle, diaStr, nBPA, am);
    otherwise
        error('fr_cases27:Map', 'unmapped contract group id %d', ...
            grp.map_id);
end
cs = fr_applyCfg27(cs, overrides);
cs = fr_finishCase27(C, kin, cs);
end

% ==========================================================================
function cs = fr_applyCfg27(cs, overrides)
% v2 sizing hook: apply a per-group config override AFTER the pattern
% builders and BEFORE fr_finishCase27 (the seed derivations in
% fr_finishCase27 use the FINAL txConfig gain and diameter). Unknown
% groups are ignored; absent fields keep the case's own values.
if isempty(overrides) || ~isa(overrides, 'containers.Map') ...
        || ~overrides.isKey(cs.name)
    return
end
ov = overrides(cs.name);
if isfield(ov, 'diaStr') && ~isempty(ov.diaStr)
    cs.diaStr = char(ov.diaStr);
    cs.diaNum = str2double(cs.diaStr);
end
if isfield(ov, 'nBPA') && ~isempty(ov.nBPA)
    cs.seed.nBPA = double(ov.nBPA);
end
fam = cs.txConfig.tackleFamily;
if strcmpi(fam, 'none'); fam = 'L'; end
if isfield(ov, 'family') && ~isempty(ov.family); fam = char(ov.family); end
g = NaN;
if isfield(cs, 'txConfig') && ~isempty(cs.txConfig.gain)
    g = max(cs.txConfig.gain);
end
if isfield(ov, 'gain') && ~isempty(ov.gain); g = double(ov.gain); end
if ~cs.hasTackle
    cs.txConfig = fr_txFor(false, cs.seed.nBPA, 'none');
else
    cs.txConfig = fr_txFor(true, cs.seed.nBPA, fam, g);
end
% optional DENSE verification grid (fr_verifyDense): resample every
% pass's sweep to ov.gridN cells BEFORE fr_fillQ27 runs in
% fr_finishCase27; targets are then resampled onto the dense grid by
% fr_targets' documented verification fallback
if isfield(ov, 'gridN') && ~isempty(ov.gridN)
    gn = round(double(ov.gridN));
    for p = 1:numel(cs.passes)
        th1o = cs.passes{p}.th1(:);
        th1n = linspace(min(th1o), max(th1o), gn)';
        cs.passes{p}.th1 = th1n';
        if isfield(cs.passes{p}, 'th2') && numel(cs.passes{p}.th2) == numel(th1o)
            cs.passes{p}.th2 = interp1(th1o, cs.passes{p}.th2(:), th1n, ...
                'linear', 'extrap')';
        end
    end
    cs.note = sprintf('%s | dense verification grid %d cells', cs.note, gn);
end
cs.note = sprintf('%s | sizing override applied: dia %s nBPA %d fam %s gain %g', ...
    cs.note, cs.diaStr, cs.seed.nBPA, fam, g);
end

% ============================== HIP ======================================
function cs = pbHip(C, kin, grp, route, seedMuscle, diaStr, nBPA, am)
% pelvis rows proximal, femur rows distal, the HIP rides T1 about the hip
% pivot (fr_chainT composition, validated exact by fr_validateKin).
rows = fr_activeRows27(route);
nRows = numel(rows);
nProx = sum(cellfun(@(p) ~strcmp(p.body, 'femur_r'), rows));
cs = fr_newCase27(grp.name, grp.group, grp.map_id, grp.stockMuscles, ...
    diaStr, nBPA, rows, grp.region);
cs.cross = [nProx + 1, nProx + 1];
cs.hasTackle = fr_tackleDecision(C, grp, route, kin, am, diaStr, ...
    fr_fixedTackle(grp.map_id));
cs.txConfig = fr_txFor(cs.hasTackle, nBPA, 'L');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];
dof = grp.primaryDof;
rg = romDeg(am, dof);
cs.passes{1} = fr_newPass27('hip_sweep', dof, '', ...
    degArr(rg, 9), zeros(1, 9), ...
    @(k, q) fr_chainT(k, q, 'femur_r', 'pelvis'), @(k, q) eye(4));
cs.movable(1) = fr_movable27(1, 'Pelvis_R', 30);       % origin
cs.movable(2) = fr_movable27(nRows, 'Femur', 30);      % insertion
if nRows >= 4
    cs.movable(3) = fr_movable27(2, 'Pelvis_R', 25);   % pelvis via (P2)
end
cs.spanKind = 'none';
cs.obstacleBones = {'Pelvis_R', 'Femur'};
cs.note = sprintf(['hip pattern; conditional rows dropped from the ' ...
    'mechanics row set (declared behavior in fr_evalRoute); tackle: %s'], ...
    fr_tackleNote(cs.hasTackle, fr_fixedTackle(grp.map_id)));
end

% ============================== TRUNK ====================================
function cs = pbTrunk(C, kin, grp, route, seedMuscle, diaStr, nBPA, am)
% STOCK two-point pelvis->torso route (geometry_20261103 seed decision:
% the stock route, NOT the robotbody legacy loop); the BACK joint rides
% T1 = torso->pelvis.
rows = route.points;
nRows = numel(rows);
cs = fr_newCase27(grp.name, grp.group, grp.map_id, grp.stockMuscles, ...
    diaStr, nBPA, rows, grp.region);
cs.cross = [nRows, nRows];
cs.hasTackle = fr_tackleDecision(C, grp, route, kin, am, diaStr, ...
    fr_fixedTackle(grp.map_id));
cs.txConfig = fr_txFor(cs.hasTackle, nBPA, 'L');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];   % lumbar extension axis z (osim back joint)
rg = romDeg(am, 'lumbar_extension');
cs.passes{1} = fr_newPass27('lumbar_sweep', 'lumbar_extension', '', ...
    degArr(rg, 9), zeros(1, 9), ...
    @(k, q) fr_chainT(k, q, 'torso', 'pelvis'), @(k, q) eye(4));
cs.movable(1) = fr_movable27(1, 'Sacrum', 25);
cs.movable(2) = fr_movable27(nRows, 'Spine', 30);
cs.spanKind = 'none';
cs.obstacleBones = {'Sacrum', 'Spine'};
cs.note = sprintf(['trunk pattern; STOCK route seed (not the robotbody ' ...
    'legacy loop); tackle: %s'], ...
    fr_tackleNote(cs.hasTackle, fr_fixedTackle(grp.map_id)));
end

% ====================== THIGH BIARTICULAR =================================
function cs = pbBi(C, kin, grp, route, seedMuscle, diaStr, nBPA, am)
% pelvis -> (femur) -> tibia biarticular; contract primary knee:
% pass 1 = knee sweep (both crossings at the knee, proximal rows held at
% the default hip), pass 2 = hip sweep (sar wired pass-A structure: femur
% rows stay femur-frame, tibia rows re-expressed in FEMUR at the default
% knee and ride T2*T1).
rows = fr_activeRows27(route);
rows = fr_fixMovingRows27(rows, route, kin, 'knee_angle_r');
nRows = numel(rows);
qDef = kin.qDefault;
XpelFem = fr_chainT(kin, qDef, 'femur_r', 'pelvis');  % femur -> pelvis
XfemTib = fr_chainT(kin, qDef, 'tibia_r', 'femur_r'); % tibia -> femur
cs = fr_newCase27(grp.name, grp.group, grp.map_id, grp.stockMuscles, ...
    diaStr, nBPA, rows, grp.region);
firstTib = find(cellfun(@(p) strcmp(p.body, 'tibia_r'), rows), 1);
firstNonPelvis = find(cellfun(@(p) ~strcmp(p.body, 'pelvis'), rows), 1);
% pass 1 (knee, contract primary): prox rows re-expressed in PELVIS at
% the default hip (they do not move with the knee); cross1 = cross2 =
% first tibia row (used-row indexing).
cs.rowUse = {true(1, nRows), true(1, nRows)};
exprA = repmat({[]}, 1, nRows);
for r = 1:nRows
    if strcmp(rows{r}.body, 'femur_r')
        exprA{r} = XpelFem;
    end
end
cs.rowExpr = {exprA, repmat({[]}, 1, nRows)};
rowFrameA = cell(1, nRows);
rowFrameB = cell(1, nRows);
for r = 1:nRows
    b = rows{r}.body;
    rowFrameA{r} = b;
    rowFrameB{r} = b;
    if strcmp(b, 'femur_r')
        rowFrameA{r} = 'pelvis';   % osim body name is pelvis (no _r)
    end
end
cs.rowFrame = {rowFrameA, rowFrameB};
cs.cross = [firstTib, firstTib];
cs.crossFor = {[firstTib, firstTib], [firstNonPelvis, nRows]};
cs.hasTackle = fr_tackleDecision(C, grp, route, kin, am, diaStr, []);
cs.txConfig = fr_txFor(cs.hasTackle, nBPA, 'L');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];
krg = C.kneeRangeDeg;
hrg = romDeg(am, 'hip_flexion_r');
cs.passes{1} = fr_newPass27('knee_sweep', 'knee_angle_r', '', ...
    degArr(krg, 9), zeros(1, 9), ...
    @(k, q) fr_chainT(k, q, 'tibia_r', 'pelvis'), @(k, q) eye(4));
cs.passes{2} = fr_newPass27('hip_sweep', 'hip_flexion_r', '', ...
    degArr(hrg, 7), zeros(1, 7), ...
    @(k, q) fr_chainT(k, q, 'femur_r', 'pelvis'), ...
    fr_constT(XfemTib));
% movables: origin, first-distal guide (femur via or distal anchor),
% insertion
cs.movable(1) = fr_movable27(1, 'Pelvis_R', 25);          % origin
if numel(rows) >= 2 && strcmp(rows{2}.body, 'femur_r')
    cs.movable(2) = fr_movable27(2, 'Femur', 25);         % femoral via
elseif numel(rows) >= 2
    cs.movable(2) = fr_movable27(2, 'Tibia', 25);         % distal anchor
end
cs.movable(end + 1) = fr_movable27(nRows, 'Tibia', 25);   % insertion
cs.spanKind = 'knee';
cs.obstacleBones = {'Pelvis_R', 'Femur'};
cs.envelopePassIdx = 1;
cs.note = sprintf(['thigh bi pattern; contract primary knee (pass 1), ' ...
    'hip coverage pass 2; conditional rows dropped, moving rows ' ...
    'fixed-tracked at the pass-1 mid pose; tackle: %s'], ...
    fr_tackleNote(cs.hasTackle, []));
end

% ======================== KNEE MONO (femur->tibia) ========================
function cs = pbKneeMono(C, kin, grp, route, seedMuscle, diaStr, nBPA, am)
% vas_int / bifemsh: femur sleeve rows, tibia insertion; knee rides T1.
% The vasti's MOVING patellar-tendon insertion row is fixed-tracked.
rows = fr_activeRows27(route);
rows = fr_fixMovingRows27(rows, route, kin, 'knee_angle_r');
nRows = numel(rows);
nProx = sum(cellfun(@(p) ~strcmp(p.body, 'tibia_r'), rows));
cs = fr_newCase27(grp.name, grp.group, grp.map_id, grp.stockMuscles, ...
    diaStr, nBPA, rows, grp.region);
cs.cross = [nProx + 1, nProx + 1];
cs.hasTackle = fr_tackleDecision(C, grp, route, kin, am, diaStr, []);
cs.txConfig = fr_txFor(cs.hasTackle, nBPA, 'L');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];
krg = C.kneeRangeDeg;
cs.passes{1} = fr_newPass27('knee_sweep', 'knee_angle_r', '', ...
    degArr(krg, 9), zeros(1, 9), ...
    @(k, q) fr_chainT(k, q, 'tibia_r', 'femur_r'), @(k, q) eye(4));
cs.movable(1) = fr_movable27(1, 'Femur', 30);           % origin
cs.movable(end + 1) = fr_movable27(nRows, 'Tibia', 25); % insertion
if nRows >= 3
    % guide row 2: cloud matches the row's owning body (femur rows ->
    % Femur cloud; a distal tibia row -> Tibia cloud)
    if strcmp(rows{2}.body, 'femur_r')
        cs.movable(end + 1) = fr_movable27(2, 'Femur', 25);
    else
        cs.movable(end + 1) = fr_movable27(2, 'Tibia', 25);
    end
end
cs.spanKind = 'knee';
cs.obstacleBones = {'Femur', 'Tibia'};
cs.note = sprintf(['knee-mono pattern; conditional rows dropped, ' ...
    'moving rows fixed-tracked at the pass-1 mid pose; tackle: %s'], ...
    fr_tackleNote(cs.hasTackle, []));
end

% ==================== GAS (femur->calcn, knee+ankle) ======================
function cs = pbGas(C, kin, grp, route, seedMuscle, diaStr, nBPA, am)
% med_gas/lat_gas: ALL declared rows kept (the condyle-wrap conditional
% row stays as a FIXED sleeve row, the gate-verified knee_bowden
% structure). Contract primary ankle: pass 1 ankle sweep (femur sleeve
% stationary, calcn insertion in TALUS frame rides T1 through the
% default knee), pass 2 = the wired knee case verbatim (insertion in
% TIBIA frame, knee rides T1).
rows = route.points;
nRows = numel(rows);
qDef = kin.qDefault;
XtalCalcn = fr_chainT(kin, qDef, 'calcn_r', 'talus_r'); % calcn -> talus
XtibCalcn = fr_chainT(kin, qDef, 'calcn_r', 'tibia_r'); % calcn -> tibia
cs = fr_newCase27(grp.name, grp.group, grp.map_id, grp.stockMuscles, ...
    diaStr, nBPA, rows, grp.region);
cs.cross = [3, 3];
cs.rowUse = {true(1, nRows), true(1, nRows)};
exprA = repmat({[]}, 1, nRows);
exprB = repmat({[]}, 1, nRows);
for r = 1:nRows
    if strcmp(rows{r}.body, 'calcn_r')
        exprA{r} = XtalCalcn;
        exprB{r} = XtibCalcn;
    end
end
cs.rowExpr = {exprA, exprB};
rowFrameA = cell(1, nRows);
rowFrameB = cell(1, nRows);
for r = 1:nRows
    b = rows{r}.body;
    rowFrameA{r} = b;
    rowFrameB{r} = b;
    if strcmp(b, 'calcn_r')
        rowFrameA{r} = 'talus_r';
        rowFrameB{r} = 'tibia_r';
    end
end
cs.rowFrame = {rowFrameA, rowFrameB};
cs.crossFor = {[nRows, nRows], [nRows, nRows]};
cs.hasTackle = true;   % wired knee_bowden exemplar (G=2, family L)
cs.txConfig = fr_txFor(true, nBPA, 'L');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];
arg = romDeg(am, 'ankle_angle_r');
krg = C.kneeRangeDeg;
cs.passes{1} = fr_newPass27('ankle_sweep', 'ankle_angle_r', '', ...
    degArr(arg, 9), zeros(1, 9), ...
    @(k, q) fr_chainT(k, q, 'talus_r', 'femur_r'), @(k, q) eye(4));
cs.passes{2} = fr_newPass27('knee_sweep', 'knee_angle_r', '', ...
    degArr(krg, 9), zeros(1, 9), ...
    @(k, q) fr_chainT(k, q, 'tibia_r', 'femur_r'), @(k, q) eye(4));
cs.movable(1) = fr_movable27(1, 'Femur', 30);          % origin
cs.movable(2) = fr_movable27(nRows, 'Calcaneus', 20);  % insertion
cs.spanKind = 'knee';
cs.obstacleBones = {'Femur', 'Tibia'};
cs.note = sprintf(['gas pattern; ALL rows kept (condyle conditional ' ...
    'row FIXED, tracked behavior; gate knee_bowden structure); pass 1 ' ...
    'contract primary ankle, pass 2 knee; tackle G=2 L (wired)']);
end

% ==================== ANKLE MONO (tibia->calcn) ===========================
function cs = pbAnkleMono(C, kin, grp, route, seedMuscle, diaStr, nBPA, am)
% soleus / tib_ant: tibia sleeve rows stationary during the ankle sweep;
% calcn insertion in TALUS frame rides T1 (subtalar at default).
rows = route.points;
nRows = numel(rows);
qDef = kin.qDefault;
XtalCalcn = fr_chainT(kin, qDef, 'calcn_r', 'talus_r');
cs = fr_newCase27(grp.name, grp.group, grp.map_id, grp.stockMuscles, ...
    diaStr, nBPA, rows, grp.region);
nProx = sum(cellfun(@(p) ~any(strcmp(p.body, ...
    {'calcn_r', 'talus_r', 'toes_r'})), rows));
cs.cross = [nProx + 1, nProx + 1];
cs.rowUse = {true(1, nRows)};
exprA = repmat({[]}, 1, nRows);
for r = 1:nRows
    if strcmp(rows{r}.body, 'calcn_r')
        exprA{r} = XtalCalcn;
    end
end
cs.rowExpr = {exprA};
rowFrameA = cell(1, nRows);
for r = 1:nRows
    b = rows{r}.body;
    rowFrameA{r} = b;
    if strcmp(b, 'calcn_r')
        rowFrameA{r} = 'talus_r';
    end
end
cs.rowFrame = {rowFrameA};
cs.hasTackle = fr_tackleDecision(C, grp, route, kin, am, diaStr, []);
cs.txConfig = fr_txFor(cs.hasTackle, nBPA, 'L');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];
arg = romDeg(am, 'ankle_angle_r');
cs.passes{1} = fr_newPass27('ankle_sweep', 'ankle_angle_r', '', ...
    degArr(arg, 9), zeros(1, 9), ...
    @(k, q) fr_chainT(k, q, 'talus_r', 'tibia_r'), @(k, q) eye(4));
cs.movable(1) = fr_movable27(1, 'Tibia', 25);            % origin
cs.movable(end + 1) = fr_movable27(nRows, 'Calcaneus', 20); % insertion
cs.spanKind = 'foot';
cs.obstacleBones = {'Tibia'};
cs.footTipRule = true;
cs.note = sprintf(['ankle-mono pattern; calcn rows in talus frame ' ...
    '(subtalar at default); tackle: %s'], fr_tackleNote(cs.hasTackle, []));
end

% ==================== SUBTALAR (tibia->calcn) =============================
function cs = pbSubtalar(C, kin, grp, route, seedMuscle, diaStr, nBPA, am)
% tib_post / per_* : tibia sleeve rows; calcn insertion rows; contract
% primary subtalar (calcn rows ride T1 = tibia<-calcn); pass 2 ankle
% coverage (calcn rows re-expressed in TALUS at the default subtalar).
rows = route.points;
nRows = numel(rows);
qDef = kin.qDefault;
XtalCalcn = fr_chainT(kin, qDef, 'calcn_r', 'talus_r');
cs = fr_newCase27(grp.name, grp.group, grp.map_id, grp.stockMuscles, ...
    diaStr, nBPA, rows, grp.region);
nProx = sum(cellfun(@(p) ~any(strcmp(p.body, ...
    {'calcn_r', 'talus_r', 'toes_r'})), rows));
cs.cross = [nProx + 1, nProx + 1];
cs.rowUse = {true(1, nRows), true(1, nRows)};
exprA = repmat({[]}, 1, nRows);
exprB = repmat({[]}, 1, nRows);
for r = 1:nRows
    if strcmp(rows{r}.body, 'calcn_r')
        exprB{r} = XtalCalcn;
    end
end
cs.rowExpr = {exprA, exprB};
rowFrameA = cell(1, nRows);
rowFrameB = cell(1, nRows);
for r = 1:nRows
    b = rows{r}.body;
    rowFrameA{r} = b;
    rowFrameB{r} = b;
    if strcmp(b, 'calcn_r')
        rowFrameB{r} = 'talus_r';
    end
end
cs.rowFrame = {rowFrameA, rowFrameB};
cs.crossFor = {[nProx + 1, nProx + 1], [nProx + 1, nProx + 1]};
cs.hasTackle = fr_tackleDecision(C, grp, route, kin, am, diaStr, []);
cs.txConfig = fr_txFor(cs.hasTackle, nBPA, 'L');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];
srg = romDeg(am, 'subtalar_angle_r');
arg = romDeg(am, 'ankle_angle_r');
cs.passes{1} = fr_newPass27('subtalar_sweep', 'subtalar_angle_r', '', ...
    degArr(srg, 9), zeros(1, 9), ...
    @(k, q) fr_chainT(k, q, 'calcn_r', 'tibia_r'), @(k, q) eye(4));
cs.passes{2} = fr_newPass27('ankle_sweep', 'ankle_angle_r', '', ...
    degArr(arg, 7), zeros(1, 7), ...
    @(k, q) fr_chainT(k, q, 'talus_r', 'tibia_r'), @(k, q) eye(4));
cs.movable(1) = fr_movable27(1, 'Tibia', 25);              % origin
cs.movable(end + 1) = fr_movable27(nRows, 'Calcaneus', 20); % insertion
cs.spanKind = 'foot';
cs.obstacleBones = {'Tibia'};
cs.footTipRule = true;
cs.note = sprintf(['subtalar pattern; pass 1 contract primary ' ...
    'subtalar, pass 2 ankle coverage; tackle: %s'], ...
    fr_tackleNote(cs.hasTackle, []));
end

% ==================== FOOT MTP (tibia->toes) ==============================
function cs = pbFoot(C, kin, grp, route, seedMuscle, diaStr, nBPA, am)
% ext_hal family: contract primary MTP: pass 1 mtp sweep (tibia sleeve +
% toes insertion; calcn guide rows DROPPED from the pass-1 mechanics row
% set - tendon guides, not sleeve hardware, wired ext_hal precedent);
% pass 2 ankle sweep (all rows, calcn+toes in talus frame).
rows = route.points;
nRows = numel(rows);
qDef = kin.qDefault;
XtalCalcn = fr_chainT(kin, qDef, 'calcn_r', 'talus_r');
XtalToes = fr_chainT(kin, qDef, 'toes_r', 'talus_r');
cs = fr_newCase27(grp.name, grp.group, grp.map_id, grp.stockMuscles, ...
    diaStr, nBPA, rows, grp.region);
nTib = sum(cellfun(@(p) strcmp(p.body, 'tibia_r'), rows));
maskA = false(1, nRows);
maskB = true(1, nRows);
for r = 1:nRows
    b = rows{r}.body;
    maskA(r) = strcmp(b, 'tibia_r') || strcmp(b, 'toes_r');
end
cs.rowUse = {maskA, maskB};
exprA = repmat({[]}, 1, nRows);
exprB = repmat({[]}, 1, nRows);
for r = 1:nRows
    b = rows{r}.body;
    if strcmp(b, 'calcn_r')
        exprB{r} = XtalCalcn;
    elseif strcmp(b, 'toes_r')
        exprB{r} = XtalToes;
    end
end
cs.rowExpr = {exprA, exprB};
rowFrameA = cell(1, nRows);
rowFrameB = cell(1, nRows);
for r = 1:nRows
    b = rows{r}.body;
    rowFrameA{r} = b;
    rowFrameB{r} = b;
    if strcmp(b, 'calcn_r') || strcmp(b, 'toes_r')
        rowFrameB{r} = 'talus_r';
    end
end
cs.rowFrame = {rowFrameA, rowFrameB};
nToes = sum(cellfun(@(p) strcmp(p.body, 'toes_r'), rows));
cs.cross = [nTib + 1, nTib + 1];
cs.crossFor = {[nTib + 1, nTib + 1], [nTib + 1, nTib + 1]};
cs.hasTackle = fr_tackleDecision(C, grp, route, kin, am, diaStr, []);
cs.txConfig = struct();
cfg = fr_txFor(cs.hasTackle, nBPA, 'none');
cfg.routingMode = {'bowden', 'moving_via'};   % wired foot structure
cs.txConfig = cfg;
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];
mrg = romDeg(am, 'mtp_angle_r');
arg = romDeg(am, 'ankle_angle_r');
cs.passes{1} = fr_newPass27('mtp_sweep', 'mtp_angle_r', '', ...
    degArr(mrg, 9), zeros(1, 9), ...
    @(k, q) fr_chainT(k, q, 'toes_r', 'tibia_r'), @(k, q) eye(4));
cs.passes{2} = fr_newPass27('ankle_sweep', 'ankle_angle_r', '', ...
    degArr(arg, 7), zeros(1, 7), ...
    @(k, q) fr_chainT(k, q, 'talus_r', 'tibia_r'), @(k, q) eye(4));
cs.movable(1) = fr_movable27(1, 'Tibia', 25);            % origin
if nTib >= 2
    cs.movable(2) = fr_movable27(2, 'Tibia', 25);        % distal sleeve anchor
end
cs.movable(end + 1) = fr_movable27(nRows, 'Toes', 15);   % insertion
cs.spanKind = 'foot';
cs.obstacleBones = {'Tibia'};
cs.footTipRule = true;
cs.envelopePassIdx = 1;
cs.note = sprintf(['foot mtp pattern; pass 1 contract primary mtp ' ...
    '(calcn guide rows dropped: tendon guides, not sleeve hardware); ' ...
    'pass 2 ankle coverage, all rows in talus frame; tackle: %s'], ...
    fr_tackleNote(cs.hasTackle, []));
end

% ==========================================================================
% finish (the gate-verified fr_finishCase pipeline, contract-aligned)
function cs = fr_finishCase27(C, kin, cs)
nRows = cs.nRows;
if ~isfield(cs, 'rowUse') || isempty(cs.rowUse)
    cs.rowUse = {true(1, nRows)};
    cs.rowExpr = {repmat({[]}, 1, nRows)};
end
if ~isfield(cs, 'rowFrame') || isempty(cs.rowFrame)
    cs.rowFrame = {};
    for p = 1:numel(cs.rowUse)
        rf = cell(1, nRows);
        for r = 1:nRows
            rf{r} = cs.rows{r}.body;
        end
        cs.rowFrame{p} = rf;
    end
end
if ~isfield(cs, 'crossFor') || isempty(cs.crossFor)
    cs.crossFor = {cs.cross};
end
if ~isfield(cs, 'secondaryCrossing')
    cs.secondaryCrossing = 0;
end
% bundleSpread (Ben 2026-11-04: "build 2 or 3 different vastus muscles
% that can each take a load of the combined vastus force and reduce
% deformation"): multi-BPA groups spread their branch ring into distinct
% corridors instead of one tight bundle -- each branch behaves as its own
% vastus unit; soft overlap still governs collisions.
if cs.groupMapId == 11 || cs.groupMapId == 24 || cs.groupMapId == 25
    cs.bundleSpread = 2.5;
end
for p = 1:numel(cs.passes)
    cs.passes{p} = fr_fillQ27(kin, cs.passes{p});
end
cs.rowOwning = cell(1, nRows);
cs.seedLoc = zeros(nRows, 3);
for r = 1:nRows
    cs.rowOwning{r} = cs.rows{r}.body;
    cs.seedLoc(r, :) = cs.rows{r}.loc;
end
% centerline-safe seed (Ben 2026-11-04: "add_mag3 should not cross
% centerline"): the add_mag3 femur-side seed crosses z=0 at high hip
% adduction. Search small offsets in the owning-body frame and keep the
% first whose worst-pose z_pelvis margin is >= 2 mm; the patch below then
% builds from the corrected seed.
if cs.groupMapId == 25
    % centerline-safe seed AT THE HOME POSITION (Ben 2026-11-04
    % correction: the seed is judged in the model's zero/neutral
    % configuration, not across the articulation sweep; the every-pose
    % z=0 rule remains the optimizer's hard constraint). If the femur-side
    % seed has z_pelvis >= 2 mm at home, the original seed stands
    % unchanged; otherwise recenter to the best home-z offset.
    % femur-side movable row (the distal attachment whose home-z is checked)
    fzRow = NaN;
    for m = 1:numel(cs.movable)
        if ~strcmpi(cs.rowOwning{cs.movable(m).row}, 'pelvis') ...
                && ~strcmpi(cs.rowOwning{cs.movable(m).row}, 'sacrum') ...
                && ~strcmpi(cs.rowOwning{cs.movable(m).row}, 'spine')
            fzRow = cs.movable(m).row;
        end
    end
    qHz = cs.passes{1}.qStructs{1};
    hzF = fieldnames(qHz);
    for i2 = 1:numel(hzF)
        qHz.(hzF{i2}) = 0;
    end
    Thz = fr_chainT(kin, qHz, 'femur_r', 'pelvis');
    zh0 = Thz(1:3, 1:3) * cs.seedLoc(fzRow, :)' + Thz(1:3, 4);
    zh0 = zh0(3);
    fprintf('FR_CASES27: add_mag3 femur seed home z_pelvis = %.1f mm', zh0 * 1000);
    fprintf(char(10));
    if zh0 * 1000 < 2
        bestZ = -inf;
        bestP = cs.seedLoc(fzRow, :);
        for ax = 1:3
            for sgn = [-1, 1]
                for dmm = [5, 10, 15, 20]
                    cand = cs.seedLoc(fzRow, :);
                    cand(ax) = cand(ax) + sgn * dmm / 1000;
                    zh = Thz(1:3, 1:3) * cand(:) + Thz(1:3, 4);
                    zh = zh(3);
                    if zh > bestZ
                        bestZ = zh;
                        bestP = cand;
                    end
                end
            end
        end
        if bestZ * 1000 >= 2
            cs.rows{fzRow}.loc = bestP;
            cs.seedLoc(fzRow, :) = bestP;
            fprintf('FR_CASES27: add_mag3 seed recentered at home (z %.1f -> %.1f mm)', zh0 * 1000, bestZ * 1000); fprintf(char(10));

        else
            fprintf('FR_CASES27: add_mag3 home z %.1f mm not fixable within 20 mm offsets - left for the optimizer', zh0 * 1000); fprintf(char(10));

        end
    end
end
for m = 1:numel(cs.movable)
    mv = cs.movable(m);
    r = mv.row;
    cs.movable(m).patch = fr_patch(C, mv.cloud, cs.rowOwning{r}, ...
        cs.seedLoc(r, :), mv.radius_mm);
end
qg = struct('passIdx', 1, 'qList', {cs.passes{1}.qStructs});
[ptsSeed, rowBSeed] = fr_locForPass(C, cs, fr_seedDesign27(cs), kin, 1);
if ~strcmp(cs.spanKind, 'none')
    cs.seedSpan = fr_spans(C, kin, ptsSeed, rowBSeed, qg.qList, ...
        cs.spanKind);
else
    midQ = cs.passes{1}.qStructs{ceil(numel(cs.passes{1}.qStructs) / 2)};
    spanMid = fr_exitSpanSeed27(cs, kin, midQ);
    cs.seedSpan = struct('spanMax_m', spanMid, ...
        'requiredSpan_m', spanMid, 'kind', 'info', ...
        'note', 'informational record; no span mandate for this route');
end
mLgrid = fr_mLGrid27(cs, kin);
mLmax = max(mLgrid);
dML = mLmax - min(mLgrid);
kmaxF = fr_kmaxFraction27(C, cs.diaStr);
G = max(max(cs.txConfig.gain), 1);
cs.fittingLength = min(C.fittingLength, max(0.008, 0.07 * mLmax));
cs.tendonSeed = 2 * cs.fittingLength + 0.005;
rest0 = max(0.03, mLmax - cs.tendonSeed - 2 * cs.fittingLength + 0.008);
dProbe = struct('L0', rest0, 'tendon', cs.tendonSeed, ...
    'seriesSplit', 0.5, 'uv', repmat([0, 0, 5], numel(cs.movable), 1));
qgProbe = struct('passIdx', 1, 'qList', {cs.passes{1}.qStructs});
[mechProbe, ~] = fr_mechBuild(C, cs, dProbe, kin, 'rigid', qgProbe);
mLclassMax = max(mechProbe.cell.MuscleLength(1:mechProbe.nG, 1));
cs.L0Seed = max(0.03, mLclassMax - cs.tendonSeed - ...
    2 * cs.fittingLength + 0.008);
if dML / (G * kmaxF) > cs.L0Seed
    fprintf(['FR_CASES27 note: %s needs rest >= %.4f m for ' ...
        'dML/(G*KMAX) but the taut seed gives %.4f m\n'], cs.name, ...
        dML / (G * kmaxF), cs.L0Seed);
end
cs.tendonHardware_m = cs.seedSpan.requiredSpan_m;
nMov = numel(cs.movable);
cs.nX = 3 * nMov + 3;
cs.xSeed = zeros(cs.nX, 1);
for m = 1:nMov
    cs.xSeed(3 * m - 2:3 * m) = [0; 0; 5];
end
cs.xSeed(3 * nMov + 1) = cs.L0Seed * 10;
cs.xSeed(3 * nMov + 2) = cs.tendonSeed * 10;
cs.xSeed(end) = 0.5;
lb = cs.xSeed; ub = cs.xSeed;
for m = 1:nMov
    lb(3 * m - 2:3 * m) = [-1; -1; C.standoff_mm(1)];
    ub(3 * m - 2:3 * m) = [1; 1; C.standoff_mm(2)];
end
lb(3 * nMov + 1) = 0.4 * cs.L0Seed * 10;
ub(3 * nMov + 1) = 2.5 * cs.L0Seed * 10;
lb(3 * nMov + 2) = 0.3 * cs.tendonSeed * 10;
ub(3 * nMov + 2) = 2.5 * cs.tendonSeed * 10;
lb(end) = 0.35; ub(end) = 0.65;
cs.lb = lb; cs.ub = ub;
cs.target = fr_targets(C, kin, cs);
end

% ==========================================================================
% helpers (self-contained; the fr_cases.m subfunctions are private)
function rows = fr_activeRows27(route)
% drop ConditionalPathPoints from the MECHANICS row set (declared
% behavior stays implemented in fr_evalRoute, tracked per pose)
rowList = {};
for i = 1:numel(route.points)
    P = route.points{i};
    if isfield(P, 'range') && ~isempty(P.range)
        continue
    end
    rowList{end + 1} = P;
end
rows = rowList;
end

function rows = fr_fixMovingRows27(rows, route, kin, qMidCoord)
% MovingPathPoints: keep the row FIXED-TRACKED - seed location = the
% declared fr_evalRoute location at the DEFAULT pose (mid-range for the
% knee primary, whose default is 0 = the pass-1 grid midpoint), owning-
% body frame (stage-1 fixed-row convention, tracked behavior). rows
% carries the already-dropped conditional set.
qMid = [];
out = rows;
for i = 1:numel(route.points)
    P = route.points{i};
    if ~isfield(P, 'fnXYZ')
        continue
    end
    % match this declared point into the kept row set by name
    nm = cellfun(@(p) p.name, rows, 'UniformOutput', false);
    idx = find(strcmp(nm, P.name), 1);
    if isempty(idx)
        continue
    end
    if isempty(qMid)
        qMid = kin.qDefault;
    end
    e = fr_evalRouteOne27(P, qMid);
    sRow = out{idx};
    sRow.loc = e;
    sRow.fixedTracked = sprintf( ...
        ['%s fixed-tracked at the default %s (declared moving ' ...
        'behavior in fr_evalRoute)'], P.name, qMidCoord);
    out{idx} = sRow;
end
end

function loc = fr_evalRouteOne27(P, q)
% declared location of ONE moving point at q (fr_evalRoute's rule)
lx = fr_evalFn27(P.fnXYZ.x, q.(P.coordXYZ{1}));
ly = fr_evalFn27(P.fnXYZ.y, q.(P.coordXYZ{2}));
lz = fr_evalFn27(P.fnXYZ.z, q.(P.coordXYZ{3}));
loc = [lx, ly, lz];
end

function y = fr_evalFn27(fn, x)
% SimmSpline natural cubic (same rule as fr_splinemm) / constant / linear
switch fn.type
    case 'spline'
        y = fr_splineEval27(fn.x, fn.y, x);
    case 'constant'
        y = fn.value;
    case 'linear'
        y = fn.coeffs(1) * x + fn.coeffs(2);
    otherwise
        error('fr_cases27:Fn', 'unsupported moving-point function');
end
end

function yi = fr_splineEval27(x, y, xi)
% natural cubic spline (tridiagonal second derivatives), one point
n = numel(x);
h = diff(x);
if n < 2
    yi = y(1);
    return
end
if n == 2
    yi = y(1) + (y(2) - y(1)) * (xi - x(1)) / h(1);
    return
end
% natural BC: second derivative 0 at both ends
A = zeros(n);
rhs = zeros(n, 1);
A(1, 1) = 1; rhs(1) = 0;
A(n, n) = 1; rhs(n) = 0;
for i = 2:n - 1
    A(i, i - 1) = h(i - 1);
    A(i, i) = 2 * (h(i - 1) + h(i));
    A(i, i + 1) = h(i);
    rhs(i) = 3 * ((y(i + 1) - y(i)) / h(i) - (y(i) - y(i - 1)) / ...
        h(i - 1));
end
m = A \ rhs;
k = find(xi <= x(2), 1);
if isempty(k)
    k = n - 1;
else
    k = max(1, min(k - 1, n - 1));
end
t = (xi - x(k)) / h(k);
a = y(k);
b = (y(k + 1) - y(k)) / h(k) - h(k) * (2 * m(k) + m(k + 1)) / 3;
c = m(k);
d = (m(k + 1) - m(k)) / (3 * h(k));
yi = a + b * (xi - x(k)) + c * (xi - x(k))^2 + d * (xi - x(k))^3;
end

function cs = fr_newCase27(name, group, mapId, stockMuscles, diaStr, ...
    nBPA, rows, region)
cs = struct();
cs.name = name;
cs.group = group;
cs.groupMapId = mapId;
cs.stockMuscles = stockMuscles;
cs.region = region;
cs.diaStr = diaStr;
cs.diaNum = str2double(diaStr);
cs.seed = struct('nBPA', nBPA);
cs.rows = rows;
cs.nRows = numel(rows);
cs.cross = [2, 2];
cs.rowUse = {};
cs.rowExpr = {};
cs.rowFrame = {};
cs.crossFor = {};
cs.passes = {};
cs.movable = struct('row', {}, 'cloud', {}, 'radius_mm', {}, ...
    'patch', {});
cs.hasTackle = false;
cs.txConfig = struct();
cs.spanKind = 'none';
cs.obstacleBones = {};
cs.footTipRule = false;
cs.bundleSpread = 1.0;
cs.envelopePassIdx = 1;
cs.secondaryCrossing = 0;
% Xi3 per-pose inside the equilibrium (Ben 2026-11-04 item 3): the
% derived class path (BiPam_pulley2_xi3pose) for the COMPLIANT branch;
% rigid is unaffected. fr_xi3_check PASSED 83/0 in a fresh process
% (tmp/xi3_check_final.log). 'legacyMaxEnvelope' remains available as
% the explicit A/B fallback in fr_mechBuild.
cs.xi3Mode = 'perPose';
cs.note = '';
end

function cfg = fr_txFor(hasTackle, nBPA, tackleFamily, gain)
% field-by-field (the struct-with-cell-value broadcast trap). gain: the
% tackle's travel/force ratio on the primary crossing (default 2; the
% class enforces gain <= tackle line parts structurally)
if nargin < 4
    gain = 2;
end
cfg = struct();
if hasTackle
    cfg.nPulleyBPA = [nBPA, 1];
    cfg.tackleLineParts = [max(1, double(gain)), 1];
    cfg.gain = [max(1, double(gain)), 1];
    cfg.routingMode = {'moving_via', 'moving_via'};
    cfg.tackleFamily = tackleFamily;
else
    cfg.nPulleyBPA = [nBPA, 1];
    cfg.tackleLineParts = [1, 1];
    cfg.gain = [1, 1];
    cfg.routingMode = {'moving_via', 'moving_via'};
    cfg.tackleFamily = 'none';
end
end

function mv = fr_movable27(row, cloud, radius_mm)
mv = struct('row', row, 'cloud', cloud, 'radius_mm', radius_mm, ...
    'patch', []);
end

function pass = fr_newPass27(name, coord1, coord2, th1, th2, T1fn, T2fn)
pass = struct();
pass.name = name;
pass.sweptCoord = coord1;
pass.sweptCoord2 = coord2;
pass.th1 = reshape(th1, 1, []);
pass.th2 = reshape(th2, 1, []);
pass.T1 = T1fn;
pass.T2 = T2fn;
pass.qStructs = {};
end

function pass = fr_fillQ27(kin, pass)
n = numel(pass.th1);
pass.qStructs = cell(1, n);
for i = 1:n
    q = kin.qDefault;
    q.(pass.sweptCoord) = pass.th1(i);
    if ~isempty(pass.sweptCoord2)
        q.(pass.sweptCoord2) = pass.th2(i);
    end
    pass.qStructs{i} = q;
end
end

function arr = degArr(degRange, n)
arr = linspace(degRange(1), degRange(2), n) * pi / 180;
arr = reshape(arr, 1, []);
end

function rg = romDeg(am, coord)
rg = am.rom_table_deg.(coord);
end

function T = fr_constT(Xm)
T = @(k, q) Xm; %#ok<NASGU>
end

function dia = fr_diaStrOf(grp)
if grp.sizing_data.ten_mm_screen_pass
    dia = '10';
else
    dia = '20';
end
end

function fixed = fr_fixedTackle(mapId)
% wired-feasibility-exemplar overrides (gate-passed configs); [] = rule
if ismember(mapId, [27, 25, 1, 13])
    fixed = true;
else
    fixed = [];
end
end

function tf = fr_tackleDecision(C, grp, route, kin, am, diaStr, fixed)
% G=2 family L iff the declared route's pass-1 length range exceeds the
% single-BPA contraction band (kmaxF * restEst); wired exemplar override
% when fixed is true/false.
if ~isempty(fixed)
    tf = fixed;
    return
end
if ~isempty(grp.primaryDof) && isfield(am.rom_table_deg, grp.primaryDof)
    rg = romDeg(am, grp.primaryDof);
    if strcmp(grp.primaryDof, 'knee_angle_r')
        rg = C.kneeRangeDeg;
    end
else
    rg = [0 0];
end
qs = cell(1, 9);
thv = linspace(rg(1), rg(2), 9) * pi / 180;
for i = 1:9
    q = kin.qDefault;
    q.(grp.primaryDof) = thv(i);
    qs{i} = q;
end
mL = zeros(1, 9);
for i = 1:9
    mL(i) = fr_evalRoute(route, kin, qs{i}).lengthM;
end
mLmax = max(mL);
dML = mLmax - min(mL);
fit = min(C.fittingLength, max(0.008, 0.07 * mLmax));
tendonEst = 2 * fit + 0.005;
restEst = max(0.03, mLmax - tendonEst - 2 * fit);
if strcmp(diaStr, '10')
    kmaxF = C.maxContractionPct.d10 / 100;
else
    kmaxF = C.maxContractionPct.d20 / 100;
end
tf = dML > kmaxF * restEst;
end

function s = fr_tackleNote(hasTackle, fixed)
if hasTackle
    src = 'rule';
    if ~isempty(fixed)
        src = 'wired-exemplar override';
    end
    s = sprintf('G=2 family L (%s)', src);
else
    if isempty(fixed)
        s = 'none (rule: dML within the contraction band)';
    else
        s = 'none (wired-exemplar override)';
    end
end
end

function d = fr_seedDesign27(cs)
nM = numel(cs.movable);
uv = repmat([0, 0, 5], nM, 1);
d = struct('L0', NaN, 'tendon', NaN, 'seriesSplit', 0.5, 'uv', uv);
end

function s = fr_exitSpanSeed27(cs, kin, q)
% exit-row-to-insertion span at q (pass 1) - informational record
T1 = cs.passes{1}.T1(kin, q);
T2 = cs.passes{1}.T2(kin, q);
crossV = cs.crossFor{1};
used = find(cs.rowUse{1});
exU = crossV(1) - 1;
inU = numel(used);
pE = cs.seedLoc(used(exU), :);
ex = cs.rowExpr{1}{used(exU)};
if ~isempty(ex)
    pE = (ex(1:3, 1:3) * pE(:) + ex(1:3, 4))';
end
pI = cs.seedLoc(used(inU), :);
exI = cs.rowExpr{1}{used(inU)};
if ~isempty(exI)
    pI = (exI(1:3, 1:3) * pI(:) + exI(1:3, 4))';
end
if inU < crossV(2)
    pI = fr_applyT27(T1, pI);
else
    pI = fr_applyT27(T2 * T1, pI);
end
s = norm(pI - pE);
end

function mL = fr_mLGrid27(cs, kin)
% patch-projected seed polyline length over the pass-1 pose list (the
% class's MuscleLength bookkeeping convention, same rule as fr_cases)
qL = cs.passes{1}.qStructs;
mL = zeros(numel(qL), 1);
for i = 1:numel(qL)
    mL(i) = fr_mLSeed27(cs, kin, qL{i});
end
end

function s = fr_mLSeed27(cs, kin, q) %#ok<INUSD>
% patch-projected seed polyline in LOCAL assigned-frame coords - exactly
% the wired fr_mLSeed bookkeeping convention (the class's MuscleLength
% grid calibrates L0Seed against the probe build, which removes any
% local-vs-grid measurement mismatch)
d = fr_seedDesign27(cs);
[P, ~] = fr_locForPass([], cs, d, kin, 1);
s = sum(sqrt(sum(diff(P).^2, 2)));
end

function p = fr_applyT27(Tm, prow)
p = (Tm(1:3, 1:3) * prow(:) + Tm(1:3, 4))';
end

function f = fr_kmaxFraction27(C, diaStr)
if strcmp(diaStr, '10')
    f = C.maxContractionPct.d10 / 100;
else
    f = C.maxContractionPct.d20 / 100;
end
end
