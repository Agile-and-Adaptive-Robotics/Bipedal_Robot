%% Figures for the Table 4.12 geometry/configuration review
% READ-ONLY rendering from Components/folded_routes_20261031.mat (recorded
% values extracted by extract_table412_geometry_20261031.m) plus the three
% reference MATs re-read for frame composition (rect_fem_r: hip+knee
% transforms; soleus_r: ankle+subtalar transforms). No class construction,
% no solving, no torque computation.
%
% Composition disclosure: the whole-leg overview folds every route into the
% PELVIS frame using each site's own recorded T at its stored pose nearest
% the neutral pose, plus the recorded rect_fem_r transforms (hip, knee at
% their stored neutral poses) for femur/tibia-rooted sites, and the recorded
% soleus_r transforms (ankle, subtalar) for the foot. psoas_r (route 27) is
% NOT overlaid: its chain root (torso) has no saved transform registering it
% to the pelvis frame; it is drawn in its own closeup panel instead.
% Skeleton links connect recorded pivots only (the rolling-knee center is
% carried inside the recorded transforms; the drawn knee vertex uses the
% recorded neutral fixed-pivot fallback, stated in the caption).

clear; clc; close all;
R = 'D:/GitHub/Bipedal_Robot';
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));
outDir = fullfile(R, 'Documentation', 'Reports and Papers', 'Dissertation', ...
    'Results', 'KompareXXVII');
figDir = fullfile(outDir, 'Figures');
cmpDir = fullfile(outDir, 'Components');
if ~isfolder(figDir), mkdir(figDir); end

L = load(fullfile(cmpDir, 'folded_routes_20261031.mat'));
f = L.folded;
fns = fieldnames(f);

selCells = cell(27, 1);
for i = 1:numel(fns)
    s = f.(fns{i});
    if s.selected
        selCells{s.actuator_id} = s;
    end
end
assert(all(~cellfun(@isempty, selCells)), 'missing selected sites');
sel = [selCells{:}];
ids = arrayfun(@(x) x.actuator_id, sel);
[ids, order] = sort(ids);
sel = sel(order);
nS = numel(ids);
assert(nS == 27, 'expected 27 selected sites, got %d', nS);

% ---- reference transforms for composition (recorded, re-read) ----------
refHip = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'BiPulley_surveys', 'BiPulley_survey2_rect_fem_r_nB3_G2_20261023_051930.mat'));
rf = refHip.spec.kin;
N1 = size(refHip.T, 3); N2 = size(refHip.T, 4);
g1 = linspace(rf.thetaRangesDeg(1, 1), rf.thetaRangesDeg(1, 2), N1);
g2 = linspace(rf.thetaRangesDeg(2, 1), rf.thetaRangesDeg(2, 2), N2);
[~, iiH] = min(abs(g1)); [~, iiiK] = min(abs(g2));
Thip = refHip.T(:, :, iiH, 1);    % femur -> pelvis at stored hip ~ 0
Tknee = refHip.T(:, :, iiiK, 2);  % tibia -> femur at stored knee ~ 0
Ttib = Thip * Tknee;              % tibia -> pelvis
fprintf(['composition pose (rect_fem_r grid): hip %.2f deg, knee %.2f deg\n'], ...
    g1(iiH), g2(iiiK));

% soleus for the foot chain (its folded pivots are already in tibia frame;
% the ankle/subtalar angles used per site are recorded in the folded struct)
refFoot = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'BiPulley_surveys', 'BiPulley_survey2_soleus_r_nB2_G1_20261023_052012.mat'));
sf = refFoot.spec.kin;
N1f = size(refFoot.T, 3); N2f = size(refFoot.T, 4);
gf1 = linspace(sf.thetaRangesDeg(1, 1), sf.thetaRangesDeg(1, 2), N1f);
gf2 = linspace(sf.thetaRangesDeg(2, 1), sf.thetaRangesDeg(2, 2), N2f);
[~, iiA] = min(abs(gf1)); [~, iiiS] = min(abs(gf2));
fprintf('foot pose (soleus_r grid): ankle %.2f deg, subtalar %.2f deg\n', ...
    gf1(iiA), gf2(iiiS));

% ---- groups (by recorded frame chain) ----------------------------------
grp = strings(nS, 1);
for k = 1:nS
    fr = strjoin(sel(k).frames, '>');
    if contains(fr, 'torsoA')
        grp(k) = "torso";
    elseif strcmp(fr, "pelvis>femur>tibia") || contains(fr, 'femurA')
        grp(k) = "hip";
    elseif strcmp(fr, "femur>tibia>tibia_rigid")
        grp(k) = "knee";
    elseif contains(fr, 'femur>tibia>calcn')
        grp(k) = "kneeankle";
    elseif contains(fr, 'torso>')
        grp(k) = "psoas";
    else
        grp(k) = "ankle";
    end
end
% group style: [color, origin marker, line style]
G = { ...
    "torso",     '#9D02D7', 's', '-';   ... % purple
    "hip",       '#EA5F94', 'o', '-';   ... % pink
    "knee",      '#FA8775', '^', '-';   ... % coral
    "kneeankle", '#FFB14E', 'v', '-';   ... % orange
    "ankle",     '#0000FF', 'd', '-';   ... % indigo
    "psoas",     '#CD34B5', 'p', '--'};    % magenta
Gname = string(G(:, 1));

% ---- compose every selected route into the pelvis frame ----------------
Pel = struct();  % sites 1..26 overlay; 27 kept separate
for k = 1:nS
    s = sel(k);
    P = s.Pw; PV = s.pivW;
    switch grp(k)
        case {"torso", "hip"}
            % already pelvis-frame
        case {"knee", "kneeankle"}
            P = P * Thip(1:3, 1:3).' + Thip(1:3, 4).';
            PV = PV * Thip(1:3, 1:3).' + Thip(1:3, 4).';
        case "ankle"
            P = P * Ttib(1:3, 1:3).' + Ttib(1:3, 4).';
            PV = PV * Ttib(1:3, 1:3).' + Ttib(1:3, 4).';
        case "psoas"
            % not composed (no saved torso->pelvis transform)
    end
    Pel(k).P = P; Pel(k).PV = PV;
end

% foot vertices in pelvis frame: ankle + subtalar pivots from soleus
kSol = find(ids == 12, 1);
ankP = sel(kSol).pivW(1, :);                          % ankle pivot (tibia frame)
subP = sel(kSol).pivW(2, :);                          % subtalar pivot (tibia frame)
toeP = max(sel(find(ids == 20, 1)).Pw, [], 1);       % toes extent (recorded points)

% =================== FIGURE 1: whole-leg overview =======================
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 17 22], ...
    'Name', 'table412 whole-leg overview');
tl = tiledlayout(fig1, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
cx = [1 3];   % (A) sagittal: horizontal = x; (B) frontal: horizontal = z
ttl = {'\bf (A) Sagittal view (x anterior, y superior)', ...
       '\bf (B) Frontal view (z right, y superior)'};
for vp = 1:2
    ax = nexttile(tl, vp);
    hold(ax, 'on');
    cH = cx(vp);   % horizontal component
    % skeleton links (recorded pivots only)
    pelO = [0 0 0];
    links = {pelO, Pel(find(ids == 8, 1)).PV(1, :); ...   % pelvis->hip
             Pel(find(ids == 8, 1)).PV(1, :), Pel(find(ids == 8, 1)).PV(2, :); ... % hip->knee
             Pel(find(ids == 8, 1)).PV(2, :), ankP; ...    % knee->ankle
             ankP, subP; ...                               % ankle->subtalar
             subP, toeP};                                  % subtalar->toes
    for e = 1:size(links, 1)
        plot(ax, [links{e, 1}(cH) links{e, 2}(cH)], ...
            [links{e, 1}(2) links{e, 2}(2)], '-', ...
            'Color', [0.55 0.55 0.55], 'LineWidth', 5, ...
            'HandleVisibility', 'off');
    end
    plot(ax, ankP(cH), ankP(2), 'k.', 'MarkerSize', 18, ...
        'HandleVisibility', 'off');
    % routes: line + separate role markers (no line+marker objects)
    % one legend handle per group present (stable G order)
    hGrp = gobjects(0); gLab = strings(0, 1);
    for gi = 1:size(G, 1)
        if strcmp(G{gi, 1}, "psoas"), continue; end
        if ~any(grp == G{gi, 1}), continue; end
        hGrp(end + 1) = plot(ax, nan, nan, G{gi, 4}, 'Color', G{gi, 2}, ...
            'LineWidth', 1.6); %#ok<AGROW>
        gLab(end + 1) = sprintf("%s (%d sites)", G{gi, 1}, ... %#ok<AGROW>
            sum(grp == G{gi, 1}));
    end
    for k = 1:nS
        if strcmp(grp(k), "psoas"), continue; end
        gi = find(Gname == grp(k));
        col = G{gi, 2};
        P = Pel(k).P;
        plot(ax, P(:, cH), P(:, 2), G{gi, 4}, 'Color', col, ...
            'LineWidth', 1.4, 'HandleVisibility', 'off');
        % index label at the origin
        text(ax, P(1, cH), P(1, 2), sprintf(' %d', ids(k)), ...
            'FontSize', 9, 'FontName', 'Arial', 'Color', col, ...
            'Interpreter', 'none');
    end
    % role markers legend (one set)
    hR = gobjects(3, 1);
    kDemo = find(ids == 8, 1);
    Pd = Pel(kDemo).P;
    hR(1) = scatter(ax, Pd(1, cH), Pd(1, 2), 70, 's', 'filled', ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', ...
        'DisplayName', 'Origin');
    hR(2) = scatter(ax, Pd(2, cH), Pd(2, 2), 55, 'o', 'filled', ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', ...
        'DisplayName', 'Via / guide');
    hR(3) = scatter(ax, Pd(end, cH), Pd(end, 2), 70, 'd', ...
        'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', ...
        'DisplayName', 'Insertion');
    axis(ax, 'equal');
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1.5, ...
        'TickLength', [0.015 0.03]);
    if vp == 1
        xlabel(ax, 'x, m');
    else
        xlabel(ax, 'z, m');
    end
    ylabel(ax, 'y, m');
    title(ax, ttl{vp}, 'Interpreter', 'tex');
    hGrp = hGrp(:); gLab = gLab(:);
    lab = [{'Origin'}; {'Via / guide'}; {'Insertion'}; cellstr(gLab)];
    lg = legend(ax, [hR; hGrp], lab, ...
        'Location', 'northeast', 'FontSize', 8, 'Box', 'on');
    set(lg, 'FontName', 'Arial');
end
annotation(fig1, 'textbox', [0.02 0.004 0.96 0.025], 'String', ...
    'Recorded pivots and route points only, no bone geometry; route 27 (psoas) appears only in its closeup.', ...
    'FontName', 'Arial', 'FontSize', 9, 'EdgeColor', 'none', ...
    'Interpreter', 'tex');
exportgraphics(fig1, fullfile(figDir, 'wholeleg_saved_routes.pdf'), 'ContentType', 'vector');
exportgraphics(fig1, fullfile(figDir, 'wholeleg_saved_routes.png'), 'Resolution', 150);
savefig(fig1, fullfile(cmpDir, 'wholeleg_saved_routes.fig'));
fprintf('WROTE Figures/wholeleg_saved_routes.(pdf|png)\n');

% =================== FIGURE 2: group closeups ===========================
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 24 15], ...
    'Name', 'table412 group closeups');
tl2 = tiledlayout(fig2, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
panels = { ...
    "torso",     '\bf (A) Torso, sites 1--3 (pelvis frame)'; ...
    "hip",       '\bf (B) Hip-crossing, sites 4--9, 24--26 (pelvis)'; ...
    "knee",      '\bf (C) Knee, sites 10--11 (femur frame)'; ...
    "kneeankle", '\bf (D) Knee+ankle, sites 13--14 (femur frame)'; ...
    "ankle",     '\bf (E) Ankle/foot, sites 12, 15--23 (tibia frame)'; ...
    "psoas",     '\bf (F) Psoas, site 27 (torso frame)'};
for pp = 1:6
    ax = nexttile(tl2, pp);
    hold(ax, 'on');
    gsel = find(grp == panels{pp, 1});
    for kk = gsel'
        s = sel(kk);
        gi = find(Gname == panels{pp, 1});
        col = G{gi, 2}; mkr = G{gi, 3};
        P = s.Pw;
        plot(ax, P(:, 1), P(:, 2), '-', 'Color', col, 'LineWidth', 1.4, ...
            'HandleVisibility', 'off');
        text(ax, P(1, 1), P(1, 2), sprintf(' %d', ids(kk)), 'FontSize', 8, ...
            'FontName', 'Arial', 'Color', col, 'Interpreter', 'none');
        scatter(ax, P(1, 1), P(1, 2), 55, mkr, 'filled', ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');
        if size(P, 1) > 2
            scatter(ax, P(2:end - 1, 1), P(2:end - 1, 2), 40, 'o', 'filled', ...
                'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
        scatter(ax, P(end, 1), P(end, 2), 55, 'd', 'filled', ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');
        if s.nBPA > 1 && ~isempty(s.bundleOffsets)
            quiver(ax, P(1, 1), P(1, 2), s.bundleOffsets(1, 1), ...
                s.bundleOffsets(1, 2), 0, 'Color', col, 'LineWidth', 1, ...
                'MaxHeadSize', 0.1, 'HandleVisibility', 'off');
        end
    end
    axis(ax, 'equal'); grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1.5, ...
        'TickLength', [0.015 0.03]);
    xlabel(ax, 'x, m'); ylabel(ax, 'y, m');
    title(ax, panels{pp, 2}, 'Interpreter', 'tex');
end
% one legend for the whole figure (role markers), in the last panel
axL = nexttile(tl2, 6);
hold(axL, 'on');
hL = gobjects(0); sL = {};
hL(end + 1) = scatter(axL, 0, 0, 60, 's', 'filled', ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Origin (fixed masters waypoint)'); %#ok<AGROW>
sL{end + 1} = hL(end).DisplayName; %#ok<AGROW>
hL(end + 1) = scatter(axL, 0, 0, 45, 'o', 'filled', ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Via / guide'); %#ok<AGROW>
sL{end + 1} = hL(end).DisplayName; %#ok<AGROW>
hL(end + 1) = scatter(axL, 0, 0, 60, 'd', 'filled', ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Insertion'); %#ok<AGROW>
sL{end + 1} = hL(end).DisplayName; %#ok<AGROW>
gl = {'torso', 'Torso group'; 'hip', 'Hip group (incl. hip-coupled)'; ...
    'knee', 'Knee group'; 'kneeankle', 'Knee+ankle group'; ...
    'ankle', 'Ankle/foot group'; 'psoas', 'Psoas (own frame)'};
for gi = 1:size(gl, 1)
    hL(end + 1) = plot(axL, nan, nan, '-', 'Color', G{gi, 2}, 'LineWidth', 2); %#ok<AGROW>
    sL{end + 1} = gl{gi, 2}; %#ok<AGROW>
end
axis(axL, 'off');
lg2 = legend(axL, hL, sL, 'Location', 'northwest', 'FontSize', 9);
set(lg2, 'FontName', 'Arial');
exportgraphics(fig2, fullfile(figDir, 'route_group_closeups.pdf'), 'ContentType', 'vector');
exportgraphics(fig2, fullfile(figDir, 'route_group_closeups.png'), 'Resolution', 150);
savefig(fig2, fullfile(cmpDir, 'route_group_closeups.fig'));
fprintf('WROTE Figures/route_group_closeups.(pdf|png)\n');
