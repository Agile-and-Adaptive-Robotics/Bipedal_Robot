function CASES = fr_cases(C, kin)
%FR_CASES Build the stage-1 demonstrator cases (ask item 5).
%
% Wired per the ask:
%   hip_flex_ext    - psoas stock route, pure hip flexion/extension sweep
%   hip_flex_extrot - psoas stock route, flexion COMBINED with external
%                     rotation (paired coordinate profile)
%   hip_rot         - adductor-group stock route, internal/external hip
%                     rotation sweep
%   knee_bowden     - medial gastrocnemius stock route crossing the knee
%                     over 0..-120 deg, tackle + moving-via, minimum
%                     central knee span enforced
%   foot_span       - extensor hallucis longus (10 mm seed) honoring the
%                     minimum foot tendon span; BPA hardware stays on the
%                     tibia (no BPA body / fitting / tackle on the foot)
%   trunk_ercspn    - erector spinae seeded per geometry_20261103's
%                     erector-spinae decision: the STOCK two-point route
%                     (attention_items.json seed_selected), NOT the
%                     robotbody 9-point legacy loop
%   sar_layering    - sartorius stock route (hip+knee biarticular), the
%                     SUPERFICIAL demonstrator for the layering proof
%
% Seeds are the stock gait2392_simbody routes parsed from the XML (master
% item 1: original OpenSim coordinates). Conditional points are dropped
% from the MECHANICS row set when not active over the whole pass (their
% declared behavior stays implemented in fr_evalRoute and is tracked per
% pose); med_gas keeps its condyle-wrap row as a FIXED sleeve row with the
% conditional behavior tracked - documented, not silent.
%
% Frame conventions per case follow BiPam_pulley2's row/cross layout with
% per-pass cross indices (crossFor) and per-pass row masks (rowUse); the
% swept joint always rides T1, whose translation is the joint pivot (see
% fr_kinematics header). Re-expression of re-framed rows uses constant
% matrices built at the pass reference pose (exact for the held
% coordinates).

caseCell = {};

% =============================== PSOAS ==================================
pso = fr_muscleroute(C.osimStock, 'psoas_r');
caseCell = fr_addCase(caseCell, fr_buildHipCase(C, kin, pso, 'psoas', 27, ...
    {'psoas_r'}, '20', 1, 'hip_flex_ext', 'hip_flexion_r', '', ...
    deg2radArr([-30 90], 9), zeros(1, 9), 'trunk/deep hip'));
caseCell = fr_addCase(caseCell, fr_buildHipCase(C, kin, pso, 'psoas', 27, ...
    {'psoas_r'}, '20', 1, 'hip_flex_extrot', 'hip_flexion_r', ...
    'hip_rotation_r', deg2radArr([-30 60], 7), ...
    max(-25, min(25, 0.4 * deg2radArr([-30 60], 7) * 180 / pi)) ...
    * pi / 180, 'trunk/deep hip'));

% ============================ ADDUCTORS =================================
adm = fr_muscleroute(C.osimStock, 'add_mag3_r');
caseCell = fr_addCase(caseCell, fr_buildHipCase(C, kin, adm, 'add_mag3', 25, ...
    {'add_mag1_r', 'add_mag2_r', 'add_mag3_r', 'add_long_r', ...
    'add_brev_r'}, '20', 2, 'hip_rot', 'hip_rotation_r', '', ...
    deg2radArr([-30 30], 7), zeros(1, 7), 'thigh'));

% ========================== KNEE / BOWDEN ================================
mga = fr_muscleroute(C.osimStock, 'med_gas_r');
caseCell = fr_addCase(caseCell, fr_buildKneeCase(C, kin, mga));

% ============================== FOOT =====================================
ehl = fr_muscleroute(C.osimStock, 'ext_hal_r');
caseCell = fr_addCase(caseCell, fr_buildFootCase(C, kin, ehl));

% ============================== TRUNK ====================================
erc = fr_muscleroute(C.osimStock, 'ercspn_r');
caseCell = fr_addCase(caseCell, fr_buildTrunkCase(C, kin, erc));

% ========================= SARTORIUS (layering) ==========================
sar = fr_muscleroute(C.osimStock, 'sar_r');
caseCell = fr_addCase(caseCell, fr_buildSarCase(C, kin, sar));

for i = 1:numel(caseCell)
    caseCell{i} = fr_finishCase(C, kin, caseCell{i});
end
% struct array (uniform field set across cases: builders share fr_newCase
% + fr_finishCase)
CASES = caseCell{1};
for i = 2:numel(caseCell)
    CASES(i) = caseCell{i};
end
end

% ==========================================================================
function cs = fr_buildHipCase(~, ~, route, group, mapId, stockMuscles, ...
    diaStr, nBPA, caseName, coord1, coord2, th1, th2, ~)
% Hip cases: pelvis rows proximal, femur rows distal, the HIP rides T1
% about the hip pivot (osim hip_r axes: flexion z, rotation y). T2 const.
rows = fr_activeRows(route);
nRows = numel(rows);
nProx = sum(cellfun(@(p) ~strcmp(p.body, 'femur_r'), rows));
hipPivot = [-0.0707, -0.0661, 0.0835];   % hip_r pelvis_offset (osim)
cs = fr_newCase(caseName, group, mapId, stockMuscles, diaStr, nBPA, ...
    rows, 'thigh');
cs.cross = [nProx + 1, nProx + 1];
cs.hasTackle = true;
% G=2 tackle at the hip crossing: the route-length variation over the hip
% ROM exceeds a single BPA's contraction band, so the tackle divides the
% travel demand (PULLEY_SWEEP_DIRECTIVE_20261021 scheme). Family L: the
% bundle-line reaction (~kN) exceeds the M family's 667 N catalog capacity.
cs.txConfig = fr_txConfig([nBPA, 1], [2, 1], [2, 1], 'moving_via', ...
    'moving_via', 'L');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];
cs.passes{1} = fr_newPass('hip_sweep', coord1, coord2, th1, th2, ...
    @(k, q) fr_hipT(q, hipPivot), @(k, q) eye(4));
cs.movable(1) = fr_movable(1, 'Pelvis_R', 30);       % origin
cs.movable(2) = fr_movable(nRows, 'Femur', 30);      % insertion
if nRows >= 4
    cs.movable(3) = fr_movable(2, 'Pelvis_R', 25);   % pelvis via (P2)
end
cs.spanKind = 'none';
cs.obstacleBones = {'Pelvis_R', 'Femur'};
end

% ==========================================================================
function cs = fr_buildKneeCase(~, kin, route)
% med_gas: femur sleeve rows (the condyle-wrap row kept FIXED, conditional
% behavior tracked), calcn insertion expressed in TIBIA (knee child frame,
% origin = knee pivot). Knee rides T1. Tackle at crossing 1 (gain [2 1]).
rows = route.points;   % ALL declared rows incl. conditional P2 (documented)
nRows = numel(rows);
XtibCalcn = fr_chainT(kin, kin.qDefault, 'calcn_r', 'tibia_r');
cs = fr_newCase('knee_bowden', 'med_gas', 13, {'med_gas_r'}, '20', 2, ...
    rows, 'lower leg');
cs.cross = [3, 3];
cs.rowUse = {true(1, nRows)};
cs.rowExpr = {repmat({[]}, 1, nRows)};
cs.rowExpr{1}{nRows} = XtibCalcn;
cs.rowFrame = {{'femur_r', 'femur_r', 'tibia_r'}};
cs.hasTackle = true;
cs.txConfig = fr_txConfig([2, 1], [2, 1], [2, 1], 'moving_via', 'moving_via', 'L');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];
cs.passes{1} = fr_newPass('knee_sweep', 'knee_angle_r', '', ...
    deg2radArr([0 -120], 9), zeros(1, 9), ...
    @(k, q) fr_chainT(k, q, 'tibia_r', 'femur_r'), @(k, q) eye(4));
cs.movable(1) = fr_movable(1, 'Femur', 30);          % origin
cs.movable(2) = fr_movable(nRows, 'Calcaneus', 20);  % insertion (calcn)
cs.spanKind = 'knee';
cs.obstacleBones = {'Femur', 'Tibia'};
end

% ==========================================================================
function cs = fr_buildFootCase(~, kin, route)
% ext_hal: tibia sleeve rows; calcn guide rows; toes insertion.
% Pass A (ankle): calcn+toes rows in TALUS (ankle child frame).
% Pass B (mtp): calcn guide rows DROPPED (documented; they are tendon
% guides, not sleeve hardware), toes rows in TOES frame, MTP rides T1
% via tibia<-toes (translation = MTP pivot in tibia coords).
rows = route.points;
nRows = numel(rows);
qDef = kin.qDefault;
XtalCalcn = fr_chainT(kin, qDef, 'calcn_r', 'talus_r');
XtalToes = fr_chainT(kin, qDef, 'toes_r', 'talus_r');
cs = fr_newCase('foot_span', 'ext_hal', 23, {'ext_hal_r'}, '10', 1, ...
    rows, 'lower leg');
cs.cross = [3, 3];
nToes = sum(cellfun(@(p) strcmp(p.body, 'toes_r'), rows));

% row masks per pass
maskA = true(1, nRows);
maskB = false(1, nRows);
for r = 1:nRows
    maskB(r) = strcmp(rows{r}.body, 'tibia_r') || ...
        strcmp(rows{r}.body, 'toes_r');
end
cs.rowUse = {maskA, maskB};
% expressions per pass (constant matrices into the pass's distal frame)
exprA = repmat({[]}, 1, nRows);
exprB = repmat({[]}, 1, nRows);
for r = 1:nRows
    b = rows{r}.body;
    if strcmp(b, 'calcn_r')
        exprA{r} = XtalCalcn;
    elseif strcmp(b, 'toes_r')
        exprA{r} = XtalToes;
    end
end
cs.rowExpr = {exprA, exprB};
rowFrameA = cell(1, nRows);
rowFrameB = cell(1, nRows);
for r = 1:nRows
    b = rows{r}.body;
    rowFrameA{r} = b;
    if strcmp(b, 'calcn_r') || strcmp(b, 'toes_r')
        rowFrameA{r} = 'talus_r';
    end
    if strcmp(b, 'toes_r')
        rowFrameB{r} = 'toes_r';
    elseif strcmp(b, 'tibia_r')
        rowFrameB{r} = 'tibia_r';
    end
end
cs.rowFrame = {rowFrameA, rowFrameB};
% pass B used rows are tibia rows then toes rows
nTibF = sum(cellfun(@(p) strcmp(p.body, 'tibia_r'), rows));
cs.crossFor = {[3, 3], [nTibF + 1, nTibF + 1]};
cs.hasTackle = false;
cs.txConfig = fr_txConfig([1, 1], [1, 1], [1, 1], 'bowden', 'moving_via', 'none');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];
cs.passes{1} = fr_newPass('ankle_sweep', 'ankle_angle_r', '', ...
    deg2radArr([-30 25], 7), zeros(1, 7), ...
    @(k, q) fr_chainT(k, q, 'talus_r', 'tibia_r'), @(k, q) eye(4));
cs.passes{2} = fr_newPass('mtp_sweep', 'mtp_angle_r', '', ...
    deg2radArr([-20 60], 7), zeros(1, 7), ...
    @(k, q) fr_chainT(k, q, 'toes_r', 'tibia_r'), @(k, q) eye(4));
cs.movable(1) = fr_movable(1, 'Tibia', 25);      % origin
cs.movable(2) = fr_movable(2, 'Tibia', 25);      % distal sleeve anchor
cs.movable(3) = fr_movable(nRows, 'Toes', 15);   % insertion
cs.spanKind = 'foot';
cs.obstacleBones = {'Tibia'};
cs.footTipRule = true;
cs.envelopePassIdx = 1;
cs.note = sprintf(['pass B drops the %d calcaneus guide rows from the ' 'mechanics row set (tendon guides, not sleeve hardware); the full ' 'route geometry with guides is evaluated by fr_evalRoute.'], nRows - nToes - sum(cellfun(@(p) strcmp(p.body, 'tibia_r'), rows)));
end

% ==========================================================================
function cs = fr_buildTrunkCase(~, ~, route)
% erector spinae: STOCK two-point route (geometry_20261103 seed_selected).
rows = route.points;
nRows = numel(rows);
cs = fr_newCase('trunk_ercspn', 'ercspn', 1, {'ercspn_r'}, '20', 3, ...
    rows, 'trunk/deep hip');
cs.cross = [nRows, nRows];
cs.hasTackle = true;
% G=2 tackle (lumbar ROM length variation); family L. NOTE: the 3x20 mm
% bundle reaction at G=2 exceeds the L family's 3003 N catalog capacity at
% full bundle force - the rated-load hard constraint flags exactly the
% handoff's documented overload case; the discrete redesign
% (load-separated branches) is a later campaign item.
cs.txConfig = fr_txConfig([3, 1], [2, 1], [2, 1], 'moving_via', ...
    'moving_via', 'L');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];   % lumbar extension axis z (osim back joint)
cs.passes{1} = fr_newPass('lumbar_sweep', 'lumbar_extension', '', ...
    deg2radArr([-10 30], 7), zeros(1, 7), ...
    @(k, q) fr_chainT(k, q, 'torso', 'pelvis'), @(k, q) eye(4));
cs.movable(1) = fr_movable(1, 'Sacrum', 25);      % P1: 5.36 mm from sacrum
cs.movable(2) = fr_movable(nRows, 'Spine', 30);   % P2: 18.63 mm from spine
cs.spanKind = 'none';
cs.obstacleBones = {'Sacrum', 'Spine'};
end

% ==========================================================================
function cs = fr_buildSarCase(~, kin, route)
% sartorius: pelvis -> femur -> tibia biarticular; the SUPERFICIAL
% demonstrator for the layering proof (deep demonstrator = psoas).
rows = route.points;
nRows = numel(rows);
qDef = kin.qDefault;
XpelFem = fr_chainT(kin, qDef, 'femur_r', 'pelvis');
XfemTib = fr_chainT(kin, qDef, 'tibia_r', 'femur_r');

firstTib = find(cellfun(@(p) strcmp(p.body, 'tibia_r'), rows), 1);
cs = fr_newCase('sar_layering', 'sar', 6, {'sar_r'}, '10', 1, ...
    rows, 'thigh');
% pass A: genuine bi - pelvis rows proximal, femur rows middle,
% tibia rows distal (cross [2 nRows]); T1 = pelvis<-femur (hip),
% T2 = femur<-tibia CONST at default knee.
% pass B: femur rows re-expressed in PELVIS at default hip (they do not
% move with the knee), cross = [firstTib firstTib]; T1 = pelvis<-tibia.
cs.rowUse = {true(1, nRows), true(1, nRows)};
exprA = repmat({[]}, 1, nRows);
exprB = repmat({[]}, 1, nRows);
for r = 1:nRows
    if strcmp(rows{r}.body, 'femur_r')
        exprB{r} = XpelFem;
    end
end
cs.rowExpr = {exprA, exprB};
rowFrameA = cell(1, nRows);
rowFrameB = cell(1, nRows);
for r = 1:nRows
    b = rows{r}.body;
    rowFrameA{r} = b;
    rowFrameB{r} = b;
    if strcmp(b, 'femur_r')
        rowFrameB{r} = 'pelvis';   % osim body name is pelvis (no _r)
    end
end
cs.rowFrame = {rowFrameA, rowFrameB};
cs.cross = [2, nRows];
cs.crossFor = {[2, nRows], [firstTib, firstTib]};
cs.hasTackle = false;
cs.txConfig = fr_txConfig([1, 1], [1, 1], [1, 1], 'moving_via', 'moving_via', 'none');
cs.primaryCrossing = 1;
cs.primaryAxis = [0; 0; 1];
cs.passes{1} = fr_newPass('hip_sweep', 'hip_flexion_r', '', ...
    deg2radArr([-30 90], 7), zeros(1, 7), ...
    @(k, q) fr_chainT(k, q, 'femur_r', 'pelvis'), @(k, q) XfemTib);
cs.passes{2} = fr_newPass('knee_sweep', 'knee_angle_r', '', ...
    deg2radArr([0 -120], 7), zeros(1, 7), ...
    @(k, q) fr_chainT(k, q, 'tibia_r', 'pelvis'), @(k, q) eye(4));
cs.movable(1) = fr_movable(1, 'Pelvis_R', 25);   % ASIS origin
cs.movable(2) = fr_movable(2, 'Femur', 25);      % femoral-triangle via
cs.movable(3) = fr_movable(nRows, 'Tibia', 20);  % pes anserinus insertion
cs.spanKind = 'none';
cs.obstacleBones = {'Pelvis_R', 'Femur', 'Tibia'};
cs.envelopePassIdx = 1;
cs.secondaryCrossing = 2;
end

% ==========================================================================
function cs = fr_finishCase(C, kin, cs)
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
% fill pose structs per pass
for p = 1:numel(cs.passes)
    cs.passes{p} = fr_fillQ(kin, cs.passes{p});
end
% seeds (row owning bodies + seed locations)
cs.rowOwning = cell(1, nRows);
cs.seedLoc = zeros(nRows, 3);
for r = 1:nRows
    cs.rowOwning{r} = cs.rows{r}.body;
    cs.seedLoc(r, :) = cs.rows{r}.loc;
end
% surface patches FIRST (the seed-span and mL measurements below must see
% the SAME patch-projected geometry the mechanics will build)
for m = 1:numel(cs.movable)
    mv = cs.movable(m);
    r = mv.row;
    cs.movable(m).patch = fr_patch(C, mv.cloud, cs.rowOwning{r}, ...
        cs.seedLoc(r, :), mv.radius_mm);
end
% seed span of record, pass 1 (resolves the contract's span placeholder
% for knee/foot routes; 'none' routes carry an informational record only)
qg = struct('passIdx', 1, 'qList', {cs.passes{1}.qStructs});
[ptsSeed, rowBSeed] = fr_locForPass(C, cs, fr_seedDesign(cs), kin, 1);
if ~strcmp(cs.spanKind, 'none')
    cs.seedSpan = fr_spans(C, kin, ptsSeed, rowBSeed, qg.qList, ...
        cs.spanKind);
else
    % informational record for routes without a span mandate: the
    % exit-row-to-insertion span at the mid pose (feeds the tendon-length
    % seed only)
    midQ = cs.passes{1}.qStructs{ceil(numel(cs.passes{1}.qStructs) / 2)};
    spanMid = fr_exitSpanSeed(cs, kin, midQ);
    cs.seedSpan = struct('spanMax_m', spanMid, 'requiredSpan_m', spanMid, ...
        'kind', 'info', ...
        'note', 'informational record; no span mandate for this route');
end
% Tendon bookkeeping seeds from the route's own length range (the class's
% TendonL is the TRANSMISSION-BOOKKEEPING parameter: contraction(theta) =
% (rest - (mL(theta) - tendon))/rest must stay inside [~0, KMAX] over the
% grid, which pins rest = dML/(c_hi - c_lo) and tendon = mL_max -
% rest*(1 - c_lo) with c_lo = 0.1*KMAX at the longest route and
% c_hi = 0.9*KMAX at the shortest). The REQUIRED-SPAN hardware quantity
% (master item 6 / contract placeholder) is a SEPARATE procurement number:
% cs.tendonHardware_m = max chord + ferrule/stop allowance at the SEED
% geometry, and fr_objective enforces the OPTIMIZED routing never needs
% more than that procured tendon.
mLgrid = fr_mLGrid(cs, kin);
mLmax = max(mLgrid);
dML = mLmax - min(mLgrid);
kmaxF = fr_kmaxFraction(C, cs.diaStr);
G = max(max(cs.txConfig.gain), 1);
% case-scaled fitting stand-in (the 25 mm probe fixture fits a ~0.7 m
% route; master item 2's fittings of record come from CAD later)
cs.fittingLength = min(C.fittingLength, max(0.008, 0.07 * mLmax));
% SEEDS DERIVED FROM THE CLASS'S OWN CLOSURE (BiPam_pulley2_20261022
% fortzBi2): s(theta) = s0ref(theta) + deltaL(theta)/G with
% s0ref = rest - mL(theta) + tendon + 2*fit. With the theta0 (longest-
% route) cell as the taut reference, s0ref = rest - mLmax + tendon + 2*fit
% is CONSTANT over the sweep (mL and deltaL track exactly), so feasibility
% (0 <= s <= KMAX*rest at every pose) pins
%   m0 = rest - mLmax + tendon + 2*fit  in [0, KMAX*rest - dML/G]
% Choose the boundary-taut seed m0 = 0 (the class clamps sLo = 0; taut is
% not slack) with the smallest bookkeeping tendon:
%   tendon = 2*fit + slack; rest = mLmax - tendon - 2*fit,
% which maximizes the contraction headroom KMAX*rest for the tackle's
% dML/G demand. Numerical check of record: dML/(G*KMAX) <= rest.
% Initial estimates of record; the solver + the class's own feasibility
% flags refine them.
cs.tendonSeed = 2 * cs.fittingLength + 0.005;
% m0 = 8 mm positive taut margin (boundary-taut is unstable against the
% tendon-stretch term of the closure: s* = m0 + dL/G - F*cSpr*nBPA/G^2).
% CALIBRATED against the class's own MuscleLength grid (a probe build with
% the current rest estimate): this removes the last millimeters of
% polyline-measurement mismatch, so the taut margin is exact.
rest0 = max(0.03, mLmax - cs.tendonSeed - 2 * cs.fittingLength + 0.008);
dProbe = struct('L0', rest0, 'tendon', cs.tendonSeed, ...
    'seriesSplit', 0.5, 'uv', repmat([0, 0, 5], numel(cs.movable), 1));
qgProbe = struct('passIdx', 1, 'qList', {cs.passes{1}.qStructs});
[mechProbe, ~] = fr_mechBuild(C, cs, dProbe, kin, 'rigid', qgProbe);
mLclassMax = max(mechProbe.cell.MuscleLength(1:mechProbe.nG, 1));
cs.L0Seed = max(0.03, mLclassMax - cs.tendonSeed - ...
    2 * cs.fittingLength + 0.008);
if dML / (G * kmaxF) > cs.L0Seed
    % headroom shortfall: shrink the bookkeeping tendon is impossible
    % (floor), so record it - the feasibility flags will show it and the
    % solver adjusts within the rest bounds
    fprintf(['FR_CASES note: %s needs rest >= %.4f m for dML/(G*KMAX) ' 'but the taut seed gives %.4f m\n'], cs.name, dML / (G * kmaxF), cs.L0Seed);
end
cs.tendonHardware_m = cs.seedSpan.requiredSpan_m;
% decision vector [per movable: u v s_mm | L0 (dm) | tendon (dm) | split]
nMov = numel(cs.movable);
cs.nX = 3 * nMov + 3;
cs.xSeed = zeros(cs.nX, 1);
for m = 1:nMov
    cs.xSeed(3 * m - 2:3 * m) = [0; 0; 5];   % u v standoff_mm(seed 5)
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
% (patches were attached before the seed measurements above)
% human target torque curve (stock model, master item 8)
cs.target = fr_targets(C, kin, cs);
end

function d = fr_seedDesign(cs)
% seed design for span/location purposes: uv at the patch centers with the
% SAME 5 mm seed standoff the decision vector carries (so the measured
% mL/span match what fr_mechBuild builds); the length fields are filled
% later in fr_finishCase
nM = numel(cs.movable);
uv = repmat([0, 0, 5], nM, 1);
d = struct('L0', NaN, 'tendon', NaN, 'seriesSplit', 0.5, 'uv', uv);
end

function s = fr_exitSpanSeed(cs, kin, q)
% exit-row-to-insertion span at q (pass 1): the tendon-side physical span
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
    pI = fr_applyTT(T1, pI);
else
    pI = fr_applyTT(T2 * T1, pI);
end
s = norm(pI - pE);
end

function mL = fr_mLGrid(cs, kin)
% rigid polyline length (the class's MuscleLength bookkeeping) over the
% pass-1 pose list
qL = cs.passes{1}.qStructs;
mL = zeros(numel(qL), 1);
for i = 1:numel(qL)
    mL(i) = fr_mLSeed(cs, kin, qL{i});
end
end

function s = fr_mLSeed(cs, kin, q)
% rigid polyline length at q over the pass-1 used rows (the class's
% MuscleLength bookkeeping), measured on the PATCH-PROJECTED seed
% geometry (fr_locForPass) so it matches what fr_mechBuild builds, with
% the class's per-row transform convention.
d = fr_seedDesign(cs);
[P, ~] = fr_locForPass([], cs, d, kin, 1);   % C unused inside
T1 = cs.passes{1}.T1(kin, q);
T2 = cs.passes{1}.T2(kin, q);
crossV = cs.crossFor{1};
nU = size(P, 1);
Pw = zeros(nU, 3);
for u = 1:nU
    if u < crossV(1)
        Pw(u, :) = P(u, :);
    elseif u < crossV(2)
        Pw(u, :) = fr_applyTT(T1, P(u, :));
    else
        Pw(u, :) = fr_applyTT(T2 * T1, P(u, :));
    end
end
s = sum(sqrt(sum(diff(P).^2, 2)));   % scalar polyline length
end

function p = fr_applyTT(Tm, prow)
p = (Tm(1:3, 1:3) * prow(:) + Tm(1:3, 4))';
end

function f = fr_kmaxFraction(C, diaStr)
if strcmp(diaStr, '10')
    f = C.maxContractionPct.d10 / 100;
else
    f = C.maxContractionPct.d20 / 100;
end
end

% ========================================================================== helpers
function cs = fr_newCase(name, group, mapId, stockMuscles, diaStr, nBPA, ...
    rows, region)
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
cs.movable = struct('row', {}, 'cloud', {}, 'radius_mm', {}, 'patch', {});  % 0x0 struct array placeholder; replaced by fr_movable on first add
cs.hasTackle = false;
cs.txConfig = struct();
cs.spanKind = 'none';
cs.obstacleBones = {};
cs.footTipRule = false;
cs.envelopePassIdx = 1;
cs.secondaryCrossing = 0;
cs.note = '';
end


function cfg = fr_txConfig(nPulleyBPA, tackleLineParts, gain, rm1, rm2, ...
    tackleFamily)
% field-by-field (a struct(...) call with the routingMode CELL value would
% broadcast txConfig into a 1x2 struct array - the struct-with-cell trap)
cfg = struct();
cfg.nPulleyBPA = nPulleyBPA;
cfg.tackleLineParts = tackleLineParts;
cfg.gain = gain;
cfg.routingMode = {rm1, rm2};
cfg.tackleFamily = tackleFamily;
end

function mv = fr_movable(row, cloud, radius_mm)
mv = struct('row', row, 'cloud', cloud, 'radius_mm', radius_mm, ...
    'patch', []);
end

function pass = fr_newPass(name, coord1, coord2, th1, th2, T1fn, T2fn)
% NB: build field by field - struct('qStructs', {}) would create a 0x0
% struct array (the struct-with-cell-value trap)
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

function pass = fr_fillQ(kin, pass)
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

function rowList = fr_activeRows(route)
% drop ConditionalPathPoints from the MECHANICS row set (stage-1: rows
% active over the whole pass; the declared conditional behavior stays
% implemented in fr_evalRoute and is tracked per pose)
rowList = {};
for i = 1:numel(route.points)
    P = route.points{i};
    if isfield(P, 'range') && ~isempty(P.range)
        continue
    end
    rowList{end+1} = P; %#ok<AGROW>
end
end

function arr = deg2radArr(degRangeOrList, n)
if numel(degRangeOrList) == 2
    arr = linspace(degRangeOrList(1), degRangeOrList(2), n) * pi / 180;
else
    arr = degRangeOrList * pi / 180;
end
arr = reshape(arr, 1, []);
end

function T = fr_hipT(q, pivot)
% hip transform about the pivot: flexion (z) then rotation (y), extrinsic,
% exactly the osim SpatialTransform axes order (rotation1 z flexion,
% rotation3 y rotation) with adduction 0.
R = fr_rotAxis([0; 0; 1], q.hip_flexion_r) * ...
    fr_rotAxis([0; 1; 0], q.hip_rotation_r);
T = eye(4);
T(1:3, 1:3) = R;
T(1:3, 4) = pivot(:);
end

function cellList = fr_addCase(cellList, cs)
% accumulate in a cell; converted to a struct array by fr_cases
cellList{end + 1} = cs;
end
