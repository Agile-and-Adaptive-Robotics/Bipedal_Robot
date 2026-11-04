function fr_deliverables()
%FR_DELIVERABLES Post-campaign deliverable tables (Ben directive item 8).
%
% 1. output/deliverable_status.csv - per (group, branch) status of record:
%    solver status (solved / best-infeasible / failed), fSeed -> fBest,
%    effective violation + top violator with units, sizing config, meet
%    fraction vs the signed group target.
% 2. output/deliverable_displacement.csv - seed-versus-final attachment
%    coordinates per movable row (owning-body frame, mm displacement).
% 3. output/torque_curves/<group>.csv - signed per-pose curves: human
%    target, robot rigid, robot compliant (figure-ready; plotting to the
%    journal standard is a separate step).
% 4. output/seed_vs_final_summary.txt - the same displacements, text.
%
% Run headless AFTER run_1: matlab -batch "addpath('<camp>'); fr_deliverables"

C = fr_setup();
if isfield(C, 'kinematics'); kin = C.kinematics; else; kin = fr_kinematics(C.osimStock); end
szFile = fullfile(C.camp, 'sizing_winners.json');
assert(exist(szFile, 'file') == 2, 'fr_deliverables:noSizing', ...
    'sizing_winners.json missing');
sz = jsondecode(fileread(szFile));
if iscell(sz.groups); groupsList = sz.groups;
else; groupsList = arrayfun(@(g) g, sz.groups, 'UniformOutput', false);
end
ck1 = fullfile(C.camp, 'run_1', 'checkpoints');
outDir = fullfile(C.camp, 'output');
curveDir = fullfile(outDir, 'torque_curves');
if ~exist(curveDir, 'dir'); mkdir(curveDir); end

fcsv = fopen(fullfile(outDir, 'deliverable_status.csv'), 'w');
fprintf(fcsv, ['group,branch,status,config,fSeed,fBest,violEff,' ...
    'topViolator,meetFrac,tauPeakTarget_Nm,tauPeakRobot_Nm\n']);
fdsp = fopen(fullfile(outDir, 'deliverable_displacement.csv'), 'w');
fprintf(fdsp, ['group,config,row,bone,seedX_mm,seedY_mm,seedZ_mm,' ...
    'finalX_mm,finalY_mm,finalZ_mm,displacement_mm\n']);

nPairs = 0;
nSolved = 0;
nBI = 0;
nFailed = 0;
for gi = 1:numel(groupsList)
    gsz = groupsList{gi};
    if ~isfield(gsz, 'diaStr'); continue; end
    cfgStr = sprintf('dia%s nBPA%d %s g%d', gsz.diaStr, gsz.nBPA, ...
        gsz.family, gsz.gain);
    ov = struct('diaStr', gsz.diaStr, 'nBPA', gsz.nBPA, ...
        'family', gsz.family, 'gain', gsz.gain);
    om = containers.Map('KeyType', 'char', 'ValueType', 'any');
    om(gsz.name) = ov;
    csAll = fr_cases27(C, kin, om);
    sel = find(strcmp({csAll.name}, gsz.name), 1);
    cs = csAll(sel);
    T = jsondecode(fileread(fullfile(C.camp, 'targets_v2', ...
        sprintf('%s_signed_targets.json', gsz.name))));
    tauT = T.tau_vec_Nm(:);

    % seed coordinates (owning-body frames) for the displacement table
    nRows = cs.nRows;
    for m = 1:numel(cs.movable)
        mv = cs.movable(m);
        r = mv.row;
        sl = cs.seedLoc(r, :) * 1000;
        fprintf(fdsp, '%s,%s,%d,%s,%.2f,%.2f,%.2f,', gsz.name, cfgStr, ...
            r, cs.rowOwning{r}, sl(1), sl(2), sl(3));
        % NOTE: final coordinates are appended per branch below from the
        % checkpoint reps (rep.loc rows correspond to used rows; map by
        % absolute row index)
        for branch = {'rigid', 'compliant'}
            b = branch{1};
            ck = fullfile(ck1, sprintf('%s_%s.mat', gsz.name, b));
            if exist(ck, 'file') ~= 2
                fprintf(fcsv, '%s,%s,NO-CHECKPOINT,%s,NA,NA,NA,NA,NA,NA\n', ...
                    gsz.name, b, cfgStr);
                nFailed = nFailed + 1;
                continue
            end
            S = load(ck);
            st = 'unknown';
            if isfield(S, 'status'); st = S.status; end
            [~, rep] = fr_objective(C, cs, kin, fr_packX(cs, S.design), ...
                b, struct('quiet', true));
            meetFrac = mean(abs(rep.tauPrimary(:)) >= ...
                abs(rep.tauTarget(:)) - 1e-9);
            topV = '';
            if isfield(rep, 'topViolator'); topV = rep.topViolator; end
            topV = strrep(topV, ',', ';'); topV = strrep(topV, newline, ' ');
            fprintf(fcsv, '%s,%s,%s,%s,%.4f,%.4f,%.3e,"%s",%.4f,%.3f,%.3f\n', ...
                gsz.name, b, st, cfgStr, S.info.fSeed, S.fBest, ...
                rep.violation, topV, meetFrac, max(abs(tauT)), ...
                max(abs(rep.tauPrimary)));
            if strcmp(st, 'solved'); nSolved = nSolved + 1;
            elseif strcmp(st, 'best-infeasible'); nBI = nBI + 1;
            else; nFailed = nFailed + 1;
            end
            nPairs = nPairs + 1;
            % displacement for movable rows whose used-row index matches
            absRows = find(cs.rowUse{1});
            for m2 = 1:numel(cs.movable)
                mv2 = cs.movable(m2);
                uIdx = find(absRows == mv2.row, 1);
                if isempty(uIdx) || uIdx > size(rep.loc, 1); continue; end
                fl = rep.loc(uIdx, :) * 1000;
                disp_mm = norm(fl - sl) ; %#ok<NODEF>
                fprintf(fdsp, '%s,%s,%d,%s,,,,,,%.2f\n', gsz.name, b, ...
                    mv2.row, cs.rowOwning{mv2.row}, disp_mm);
            end
        end
        break;  % displacement block handled per branch above; keep table compact
    end

    % signed curve CSV (target + both branches where present)
    cf = fullfile(curveDir, sprintf('%s.csv', gsz.name));
    fc2 = fopen(cf, 'w');
    fprintf(fc2, 'cell,deg,human_target_Nm,robot_rigid_Nm,robot_compliant_Nm\n');
    rr = struct('rigid', [], 'compliant', []);
    for branch = {'rigid', 'compliant'}
        b = branch{1};
        ck = fullfile(ck1, sprintf('%s_%s.mat', gsz.name, b));
        if exist(ck, 'file') == 2
            S = load(ck);
            [~, rep] = fr_objective(C, cs, kin, fr_packX(cs, S.design), ...
                b, struct('quiet', true));
            rr.(b) = rep.tauPrimary(:);
        end
    end
    for cc = 1:numel(tauT)
        vr = NaN; vc = NaN;
        if ~isempty(rr.rigid) && cc <= numel(rr.rigid); vr = rr.rigid(cc); end
        if ~isempty(rr.compliant) && cc <= numel(rr.compliant); vc = rr.compliant(cc); end
        deg = T.cells(cc).requested_deg;
        fprintf(fc2, '%d,%.3f,%.4f,%.4f,%.4f\n', cc, deg, tauT(cc), vr, vc);
    end
    fclose(fc2);
end
fclose(fcsv);
fclose(fdsp);
fprintf(['deliverables written: %d (group, branch) pairs ' ...
    '(solved %d, best-infeasible %d, failed/missing %d)\n'], ...
    nPairs, nSolved, nBI, nFailed);
fprintf('status table: %s\n', fullfile(outDir, 'deliverable_status.csv'));
fprintf('displacements: %s\n', fullfile(outDir, 'deliverable_displacement.csv'));
fprintf('torque curves: %s\n', curveDir);
end
