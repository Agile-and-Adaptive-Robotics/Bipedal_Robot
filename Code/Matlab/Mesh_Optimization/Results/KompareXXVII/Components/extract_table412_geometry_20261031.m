%% Extract the saved geometry + configurations behind dissertation Table 4.12
% READ-ONLY provenance extraction. Loads each of the 54 saved survey MATs
% named by tmp/table412_geom_20261031/_fetched/
% goal_table_torque_active_corrected_20261028.csv, reads recorded fields,
% and writes route_configurations.csv + route_points.csv + folded geometry
% for the figures. NO class construction, NO solving, NO refitting, NO
% torque computation: the only arithmetic is folding the recorded route
% points through the RECORDED T grid to one stated pose for drawing
% (same composition the class's mapVectorToFrame performs).
%
% Frame convention (Ben, 2026-10-27): x anterior, y superior, z right;
% knee flexion negative. Points are body-local per spec.frames.
%   T(:,:,ii,1): middle frame -> proximal frame at joint-1 grid angle ii
%   T(:,:,iii,2): distal frame -> middle frame at joint-2 grid angle iii
% Pose of record for the folded figures: the stored grid pose nearest
% theta1 = 0 deg, theta2 = 0 deg, with the exact stored angles recorded.

clear; clc;
R = 'D:/GitHub/Bipedal_Robot';
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));  % class on path for MAT load only
addpath(fullfile(R, 'Code', 'Matlab', 'Functions', 'ModernRobotics'));  % RowVecTrans
csvIn = fullfile(R, 'tmp', 'table412_geom_20261031', '_fetched', ...
    'goal_table_torque_active_corrected_20261028.csv');
outDir = fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'table412_geometry_20261030');
cmpDir = fullfile(outDir, 'Components');
if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(cmpDir), mkdir(cmpDir); end

T = readtable(csvIn, 'TextType', 'string');
nS = height(T);
assert(nS == 27, 'roster must have 27 rows, got %d', nS);

cfgRows = {};   % configuration rows (54 = 27 sites x 2 sides)
ptsRows = {};   % route-point rows
folded = struct();  % per-site folded geometry for the figures

for i = 1:nS
    site = char(T.site(i));
    muscle = char(T.muscle(i));
    fid = T.actuator_id(i);
    for side = 1:2
        if side == 1
            matPath = char(T.original_mat(i)); tag = 'baseline';
        else
            matPath = char(T.updated_mat(i));  tag = 'selected';
        end
        matPath = strrep(matPath, 'D:\Github\', 'D:\GitHub\');
        matPath = strrep(matPath, 'D:/Github/', 'D:/GitHub/');
        assert(isfile(matPath), 'missing MAT: %s', matPath);
        S = load(matPath);                    % read-only
        o = S.g.obj;                          % saved class object (fields only)
        sp = S.spec; sm = S.surveyMeta;

        % ---- configuration fields (recorded values, verbatim) ----------
        cfg = struct();
        cfg.actuator_id = fid; cfg.site = site; cfg.muscle = muscle;
        cfg.side = tag;
        cfg.mat_file = strrep(matPath, '\\', '/');
        cfg.spec_name = gf(sp, 'name', '');
        cfg.spec_source = gf(sp, 'source', '');
        cfg.master_name = gf(sp, 'masterName', '');
        cfg.diameter_nominal_mm = gf(sm, 'bpa_size_mm', gf(sp, 'diameter', NaN));
        cfg.diameter_saved_mm = o.Diameter;            % actual saved object
        cfg.nBPA = gf(sm, 'nBPA', NaN);
        cfg.G_survey = gf(sm, 'G', NaN);               % travel multiplier of record
        cfg.pulley_nPulleyBPA_c1 = sp.pulley.nPulleyBPA(1);
        cfg.pulley_nPulleyBPA_c2 = sp.pulley.nPulleyBPA(min(2, numel(sp.pulley.nPulleyBPA)));
        cfg.tackleLineParts_c1 = sp.pulley.tackleLineParts(1);
        cfg.tackleLineParts_c2 = sp.pulley.tackleLineParts(min(2, numel(sp.pulley.tackleLineParts)));
        cfg.routingMode = strjoin(string(sp.pulley.routingMode), '/');
        cfg.resting_length_m = o.RestingL;
        cfg.tendon_length_m = o.TendonL;
        cfg.kmax = o.Kmax;
        cfg.fitting_length_m = o.FittingLength;
        cfg.pressure_kPa = o.Pressure;
        cfg.Xi0_m = o.Xi0; cfg.Xi1_N_per_m = o.Xi1; cfg.Xi2_N_per_m = o.Xi2;
        cfg.wraps = o.Wraps;
        cfg.fmax_N = gf(sp, 'fmax', NaN);
        cfg.fmax_source = gf(sp, 'fmaxSource', '');
        cfg.n_points = size(sp.points, 1);
        cfg.frames = strjoin(string(sp.frames), '>');
        cfg.cross = strtrim(sprintf('%d ', sp.cross));
        cfg.bundle_offsets = sprintf('[%s]', strjoin( ...
            string(num2str(sp.bundleOffsets, '%.6g')), '; '));
        cfg.n_bundle_offsets = size(sp.bundleOffsets, 1);
        cfg.theta1_range_deg = sprintf('[%g %g]', sp.kin.thetaRangesDeg(1, :));
        cfg.theta2_range_deg = sprintf('[%g %g]', sp.kin.thetaRangesDeg(2, :));
        cfg.primary_dof = gf(sm, 'primary_dof', '');
        cfg.secondary_dof = gf(sm, 'secondary_dof', '');
        cfg.grading_slot = gf(sm, 'slot', NaN);
        cfg.domain_cells = gf(sm, 'domainCells', NaN);
        cfg.excluded_cells = gf(sm, 'excludedCells', NaN);
        cfg.eval_budget = gf(sm, 'evals', NaN);
        cfg.spec_note = gf(sp, 'notes', '');
        cfg.surveymeta_note = gf(sm, 'note', '');
        cfg.kin_joint_names = strjoin(string(sp.kin.jointNames), ',');
        cfg.kin_types = strjoin(string(sp.kin.types), ',');

        % ---- route-point rows (recorded body-local coordinates) --------
        nP = size(sp.points, 1);
        nFr = numel(sp.frames);
        for r = 1:nP
            fi = rowFrameIdx(r, sp.cross);       % class rowFrame rule
            fi = min(fi, nFr);
            pr = struct();
            pr.actuator_id = fid; pr.site = site; pr.side = tag;
            pr.point_index = r;
            if r == 1
                pr.role = 'origin';
            elseif r == nP
                pr.role = 'insertion';
            else
                pr.role = 'via';
            end
            pr.frame = char(sp.frames(fi));
            pr.frame_index = fi;
            pr.x_m = sp.points(r, 1);
            pr.y_m = sp.points(r, 2);
            pr.z_m = sp.points(r, 3);
            pr.point_class = 'fixed (masters waypoint)';
            pr.source = sprintf('%s:spec.points(%d,:) + spec.frames{%d}', ...
                cfg.mat_file, r, fi);
            ptsRows{end + 1, 1} = pr; %#ok<AGROW>
        end

        % ---- side comparison (recorded arrays, isequal only) -----------
        if side == 1
            base.points = sp.points; base.frames = string(sp.frames);
            base.cross = sp.cross; base.pulley = sp.pulley;
            base.diameter = o.Diameter; base.kin = sp.kin;
            base.rest = o.RestingL; base.tendon = o.TendonL;
            cfg.points_identical_to_baseline = "n/a (baseline side)";
            cfg.diameter_identical_to_baseline = "n/a (baseline side)";
            cfg.pulley_identical_to_baseline = "n/a (baseline side)";
            cfg.rest_delta_m = 0; cfg.tendon_delta_m = 0;
        else
            cfg.points_identical_to_baseline = string( ...
                isequal(sp.points, base.points) && ...
                isequal(string(sp.frames), base.frames) && ...
                isequal(sp.cross, base.cross) && ...
                isequal(sp.kin.thetaRangesDeg, base.kin.thetaRangesDeg) && ...
                isequal(sp.kin.pivots, base.kin.pivots) && ...
                isequal(sp.kin.axes, base.kin.axes));
            cfg.diameter_identical_to_baseline = string(o.Diameter == base.diameter);
            cfg.pulley_identical_to_baseline = string(isequal(sp.pulley, base.pulley));
            cfg.rest_delta_m = o.RestingL - base.rest;
            cfg.tendon_delta_m = o.TendonL - base.tendon;
        end

        % ---- folded geometry for drawing (recorded T, stated pose) -----
        % Stored pose nearest theta1 = 0, theta2 = 0:
        th1 = sp.kin.thetaRangesDeg(1, :);
        th2 = sp.kin.thetaRangesDeg(2, :);
        N1 = size(S.T, 3); N2 = size(S.T, 4);
        g1 = linspace(th1(1), th1(2), N1);
        g2 = linspace(th2(1), th2(2), N2);
        [~, ii] = min(abs(g1)); [~, iii] = min(abs(g2));
        T1 = S.T(:, :, ii, 1);   % frame2 -> frame1
        T2 = S.T(:, :, iii, 2);  % frame3 -> frame2
        nP = size(sp.points, 1);
        Pw = nan(nP, 3);
        for r = 1:nP
            fi = rowFrameIdx(r, sp.cross);
            v = sp.points(r, :);
            switch fi
                case 1, Pw(r, :) = v;
                case 2, Pw(r, :) = RowVecTrans(T1, v);
                otherwise
                    v2 = RowVecTrans(T2, v);
                    Pw(r, :) = RowVecTrans(T1, v2);
            end
        end
        % pivots to proximal frame: pivot1 as-is; pivot2 via T1
        pivW = [sp.kin.pivots(1, :); RowVecTrans(T1, sp.kin.pivots(2, :))];
        f = struct();
        f.site = site; f.actuator_id = fid; f.muscle = muscle;
        f.selected = (side == 2);
        f.frames = string(sp.frames);
        f.Pw = Pw; f.pivW = pivW;
        f.ii = ii; f.iii = iii;
        f.theta1_used_deg = g1(ii); f.theta2_used_deg = g2(iii);
        f.theta_ranges = sp.kin.thetaRangesDeg;
        f.bundleOffsets = sp.bundleOffsets;
        f.diameter_mm = o.Diameter; f.nBPA = gf(sm, 'nBPA', NaN); f.G = gf(sm, 'G', NaN);
        f.rest = o.RestingL; f.tendon = o.TendonL;
        fn = sprintf('fold_%02d_%s_%s', fid, site, tag);
        folded.(fn) = f;

        cfgRows{end + 1, 1} = cfg; %#ok<AGROW>
        fprintf('%02d %-10s %-9s rest %.4f tendon %.4f nB%d G%d d%gmm pts%d\n', ...
            fid, site, tag, o.RestingL, o.TendonL, gf(sm, 'nBPA', -1), ...
            gf(sm, 'G', -1), o.Diameter, nP);
    end
    % side comparison sanity: every site must have inherited routing
    assert(isfield(cfg, 'points_identical_to_baseline'), 'no comparison for %s', site);
end

cfgT = struct2table([cfgRows{:}]);
writetable(cfgT, fullfile(outDir, 'route_configurations.csv'));
ptsT = struct2table([ptsRows{:}]);
writetable(ptsT, fullfile(outDir, 'route_points.csv'));
save(fullfile(cmpDir, 'folded_routes_20261031.mat'), 'folded', '-v7');
fprintf('\nWROTE route_configurations.csv (%d rows), route_points.csv (%d rows)\n', ...
    height(cfgT), height(ptsT));
fprintf('WROTE Components/folded_routes_20261031.mat (%d folded sites)\n', ...
    numel(fieldnames(folded)));

%% ---- local helper -----------------------------------------------------------
function val = gf(s, name, default)
% Guarded field read: recorded value when present, default + flag otherwise.
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end

function f = rowFrameIdx(r, cross)
% Verbatim BiPam_pulley2 rowFrame rule: rows before Cross(1) live in the
% proximal frame; duplicate crossings (Cross(1) == Cross(2)) have an empty
% middle frame and rows at or beyond Cross(1) live in the distal frame.
if r < cross(1)
    f = 1;
elseif cross(1) == cross(2)
    f = 3;
elseif r < cross(2)
    f = 2;
else
    f = 3;
end
end
