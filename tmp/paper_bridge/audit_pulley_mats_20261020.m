% audit_pulley_mats_20261020.m -- READ-ONLY pulley machinery audit (EDIT_SPEC
% 2026-10-20 section 4). Loads the 09-26 batch mats, the mono bifemsh mat and
% the checkpoint; prints the pulley-config fields verbatim. Writes NOTHING.
% (Scratch per spec rule: tmp/paper_bridge/ only. No save, no savepath.)

root = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');

fprintf('=== HOST: %s | %s ===\n', getenv('COMPUTERNAME'), version);

%% 0. whos -file inventory (no class path needed for this)
batchFiles = { ...
    'BiPulley_opensim_med_gas_r_20260926_1655.mat'; ...
    'BiPulley_opensim_soleus_r_20260926_1657.mat'; ...
    'BiPulley_opensim_tib_ant_r_20260926_1659.mat'; ...
    'BiPulley_opensim_rect_fem_r_20260926_1701.mat'; ...
    'BiPulley_opensim_bifemsh_r_20260926_1702.mat'};
fprintf('\n########## 0. whos -file inventory ##########\n');
for k = 1:numel(batchFiles)
    W = whos('-file', fullfile(resDir, batchFiles{k}));
    fprintf('--- %s ---\n', batchFiles{k});
    for j = 1:numel(W)
        fprintf('    %-14s %-10s %s\n', W(j).name, W(j).class, mat2str(W(j).size));
    end
end
W = whos('-file', fullfile(resDir, 'Bifemsh_20mm_Result_pulley_20260926_1652.mat'));
fprintf('--- Bifemsh_20mm_Result_pulley_20260926_1652.mat (first 40 vars of %d) ---\n', numel(W));
for j = 1:min(40, numel(W))
    fprintf('    %-18s %-12s %s\n', W(j).name, W(j).class, mat2str(W(j).size));
end

%% 1. The five batch mats: pulley config + outcome fields, verbatim
fprintf('\n########## 1. Batch mats: pulley fields verbatim ##########\n');
for k = 1:numel(batchFiles)
    S = load(fullfile(resDir, batchFiles{k}));   %#ok<NASGU>
    fprintf('===== %s =====\n', batchFiles{k});
    fo = S.finalObj;
    sp = S.spec;
    fprintf('spec: name=%s source=%s fmax=%g N (%s)\n', sp.name, sp.source, ...
        sp.fmax, sp.fmaxSource);
    fprintf('spec.pulley (builder input): nPulleyBPA=%s tackleLineParts=%s routingMode=%s\n', ...
        mat2str(sp.pulley.nPulleyBPA), mat2str(sp.pulley.tackleLineParts), ...
        mat2str(sp.pulley.routingMode));
    fprintf('route: frames={%s, %s, %s} cross=%s npoints=%d\n', sp.frames{1}, ...
        sp.frames{2}, sp.frames{3}, mat2str(sp.cross), size(sp.points, 1));
    fprintf('kin: types={%s, %s} kneeSlot=%d thetaRangesDeg=%s\n', sp.kin.types{1}, ...
        sp.kin.types{2}, sp.kin.kneeSlot, mat2str(sp.kin.thetaRangesDeg));
    fprintf('finalObj pulley fields: NPulleyBPA=%s TackleLineParts=%s PulleyGain=%s\n', ...
        mat2str(fo.NPulleyBPA), mat2str(fo.TackleLineParts), mat2str(fo.PulleyGain));
    fprintf('  RoutingMode={%s, %s} PulleyExitIndex=%s\n', fo.RoutingMode{1}, ...
        fo.RoutingMode{2}, mat2str(fo.PulleyExitIndex));
    fprintf('  BowdenBossDia=%g m BowdenRunClearance=%g m\n', fo.BowdenBossDia, ...
        fo.BowdenRunClearance);
    fprintf('  PulleyActive=%s ActiveCrossing=%d Pressure=%g kPa Diameter=%g mm\n', ...
        mat2str(fo.PulleyActive), fo.ActiveCrossing, fo.Pressure, fo.Diameter);
    fprintf('  RestingL=%.4f m TendonL=%.4f m Kmax=%.4f Wraps=%g Xi0=%g Xi1=%g Xi2=%g\n', ...
        fo.RestingL, fo.TendonL, fo.Kmax, fo.Wraps, fo.Xi0, fo.Xi1, fo.Xi2);
    fprintf('  grid Torque_ins size=%s ; PulleyInfeasible cells=%d ; PulleySlack cells=%d (of %d)\n', ...
        mat2str(size(fo.Torque_ins)), sum(fo.PulleyInfeasible(:)), ...
        sum(fo.PulleySlack(:)), numel(fo.PulleyInfeasible));
    aK = fo.ActiveCrossing;
    tau = squeeze(abs(fo.Torque_ins(:, :, aK)));
    fprintf('  |Torque_ins| @active crossing %d: min=%.4g max=%.4g N*m (NaNs: %d of %d)\n', ...
        aK, min(tau(:)), max(tau(:)), sum(isnan(tau(:))), numel(tau));
    otherK = 3 - aK;
    tauO = squeeze(abs(fo.Torque_ins(:, :, otherK)));
    fprintf('  |Torque_ins| @other crossing %d: min=%.4g max=%.4g N*m (NaNs: %d of %d)\n', ...
        otherK, min(tauO(:)), max(tauO(:)), sum(isnan(tauO(:))), numel(tauO));
    if isfield(S, 'xBest')
        fprintf('driver outcome: xBest=%s exitflagS=%g exitflagP=%g\n', ...
            mat2str(S.xBest(:).', 5), S.exitflagS, S.exitflagP);
    end
    if isfield(S, 'st')
        fprintf('driver outcome: margin=%+.6g worstC=%+.6g tauTarget=%g N*m\n', ...
            S.st.margin, S.st.worstC, sp.fmax * 0.05);
    end
    fprintf('notes: %s\n\n', sp.notes);
end

%% 2. Checkpoint
fprintf('########## 2. BiPulley_checkpoint.mat ##########\n');
C = load(fullfile(resDir, 'BiPulley_checkpoint.mat'));
fprintf('doneKeys = ');
for j = 1:numel(C.doneKeys)
    fprintf('%s  ', C.doneKeys{j});
end
fprintf('\n\n');

%% 3. Mono bifemsh mat (Opt_run_pulley product)
fprintf('########## 3. Bifemsh_20mm_Result_pulley_20260926_1652.mat ##########\n');
M = load(fullfile(resDir, 'Bifemsh_20mm_Result_pulley_20260926_1652.mat'));
c = M.ctx;
fprintf('ctx pulley knobs: pulleyBPACount=%g tackleLineParts=%g pulleyRoutingMode=%s\n', ...
    c.pulleyBPACount, c.tackleLineParts, c.pulleyRoutingMode);
fprintf('  pulleyExitIndex=%g bowdenBossDia=%g bowdenRunClearance=%g\n', ...
    c.pulleyExitIndex, c.bowdenBossDia, c.bowdenRunClearance);
fprintf('  optimizePulleyGain=%d pulleyGainBounds=%s pulleyConfigSweep=%d\n', ...
    c.optimizePulleyGain, mat2str(c.pulleyGainBounds), c.pulleyConfigSweep);
fprintf('ctx pressure fields: ');
fn = fieldnames(c);
for j = 1:numel(fn)
    if ~isempty(regexpi(fn{j}, 'pres', 'once'))
        v = c.(fn{j});
        if isnumeric(v) && numel(v) <= 12
            fprintf('%s=%s  ', fn{j}, mat2str(v));
        else
            fprintf('%s=<nonscalar>  ', fn{j});
        end
    end
end
fprintf('\n');
pb = M.predBest;
fprintf('predBest fields: ');
disp(fieldnames(pb).');
if isfield(pb, 'pulleyConfig')
    pc = pb.pulleyConfig;
    fprintf('predBest.pulleyConfig: ');
    pcn = fieldnames(pc);
    for j = 1:numel(pcn)
        v = pc.(pcn{j});
        if ischar(v) || (isstring(v))
            fprintf('%s=%s  ', pcn{j}, char(v));
        elseif isnumeric(v) && numel(v) <= 4
            fprintf('%s=%s  ', pcn{j}, mat2str(v));
        end
    end
    fprintf('\n');
end
if isfield(pb, 'pulleyGain')
    fprintf('predBest.pulleyGain = %s\n', mat2str(pb.pulleyGain));
end
% Torque curve inventory in the mono result
fn = fieldnames(pb);
for j = 1:numel(fn)
    v = pb.(fn{j});
    if isnumeric(v) && numel(v) > 4 && ndims(v) <= 3 && any(~isempty(regexpi(fn{j}, 'torq|tau|theta|angle|orient', 'once')))
        fprintf('predBest.%s: size=%s min=%.4g max=%.4g\n', fn{j}, mat2str(size(v)), ...
            min(v(:)), max(v(:)));
    end
end
% Opt_run_pulley context geometry: orientation grid size
if isfield(c, 'geo') && isfield(c.geo, 'nOrientations')
    fprintf('ctx.geo.nOrientations = %s\n', mat2str(c.geo.nOrientations));
end
% Objective/margin outcome if stored
if isfield(M, 'fBest')
    fprintf('fBest = %.6g\n', M.fBest);
end
if isfield(pb, 'worstMarginPct')
    fprintf('predBest.worstMarginPct = %.6g %%\n', pb.worstMarginPct);
end
% what pressure does the mono objective sweep?
if isfield(c, 'pressures')
    fprintf('ctx.pressures = %s kPa\n', mat2str(c.pressures));
end
fprintf('\nDone. (Read-only audit; nothing written.)\n');
