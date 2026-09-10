%Dig_resultsTable.m — master pivot-ready results table (long format, one row per
%solution x test case; all 3 GoF + baseline columns + ratios). FVU is the primary
%GoF per Ben; MaxResidual informational. Regenerate after new runs.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
R = {};   %long rows: {Category, RunID, Evaluator, Convention, Pbr, Xi0_m, Xi1, Xi2, Xi3, Test, RMSE, FVU, MaxResidual, RMSEbase, FVUbase, MaxRbase, RMSEratio, FVUratio, Notes}

%% A. Flexor pinned CV picks (4 driver mats, per-test + baseline)
cfgs = { ...
 'minimizeFlxPin10_results_20260908_2brkt_1trans_noT3.mat',   '1trans', 'noT3 (new bounds)';
 'minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat',   '2trans', 'noT3 (new bounds)';
 'minimizeFlxPin10_results_20260908_2brkt_1trans_noT3noT5.mat','1trans','noT3noT5 (old bounds)';
 'minimizeFlxPin10_results_20260908_2brkt_2trans_noT3noT5.mat','2trans','noT3noT5 (old bounds)'};
for c = 1:4
    S = load(cfgs{c,1}, 'k1','k2','k3','f','a0','labels');
    lab = string(S.labels);
    rid = ['minimizeFlxPin10mm_2brk ', char(cfgs{c,3})];
    for j = 1:5
        R(end+1,:) = {'FlexorPin CV pick=1', rid, 'minimizeFlxPin2brk', cfgs{c,2}, 'Pbri cantilever + corr47(exp)', ...
            S.k1, S.k2, S.k3, NaN, char(lab(j)), S.f(j,1), S.f(j,2), S.f(j,3), ...
            '', '', '', 'Xi1/Xi2 locked to this pick in the extensor runs'}; %#ok<AGROW>
        R(end+1,:) = {'FlexorPin baseline', rid, 'minimizeFlxPin2brk', cfgs{c,2}, 'Pbri cantilever + corr47(exp)', ...
            0, NaN, NaN, NaN, char(lab(j)), S.a0(j,1), S.a0(j,2), S.a0(j,3), ...
            '', '', '', 'rigid baseline'}; %#ok<AGROW>
    end
end

%% B. Extensor CV picks (2 front mats: rib midpoint + lower bolt hole), per-test + baseline
a0e = minimizeExtX3(0, Inf, Inf, 0);   %pinned ext baseline, 9 tests
lab9 = ["40cm","40cm-tendon","42cm","42cm-tendon","43cm","43cm-tendon","46cm","47cm","48cm"];
exts = { ...
 'minimizeExt10mmX3_results_20260910_noT3.mat', '2trans', 'rib midpoint [-3.84,-46.44,62.5]';
 'minimizeExt10mmX3_results_20260909_1trans_Lbolt.mat', '1trans', 'lower bolt hole [-6.26,-29.69,75.06]'};
for c = 1:2
    S = load(exts{c,1}, 'sol_actual', 'g', 'filtered_results', 'results_sort_actual');
    sa = S.sol_actual; g = S.g;
    rid = ['minimizeExt10mmX3 ', char(exts{c,2})];
    f_all = minimizeExtX3(sa(1), sa(2), sa(3), sa(4));   %all 9 tests
    for j = 1:9
        R(end+1,:) = {'ExtensorPin CV pick=1', rid, 'minimizeExtX3', '2trans', exts{c,3}, ...
            sa(1), sa(2), sa(3), sa(4), char(lab9(j)), f_all(j,1), f_all(j,2), f_all(j,3), ...
            a0e(j,1), a0e(j,2), a0e(j,3), f_all(j,1)/a0e(j,1), f_all(j,2)/a0e(j,2), 'Xi3>0.02 filter applies'}; %#ok<AGROW>
    end
    % bio-ext validation of the same pick
    b = minimizeExt(sa(1), sa(2), sa(3), sa(4), 1);
    R(end+1,:) = {'ExtensorPin CV pick=1', rid, 'minimizeExt', '2trans', exts{c,3}, ...
        sa(1), sa(2), sa(3), sa(4), 'bio-ext 52cm', b(1), b(2), b(3), ...
        3.032, 4.486, 5.587, b(1)/3.032, b(2)/4.486, 'further validation'}; %#ok<AGROW>
end

%% C. Bio-flexor candidate evaluations (329 unique x 3 cases, from dub_filt)
Dg = load('Dig_out/Dig_FlxBio_dubfilt_results.mat', 'dub_filt_results', 'met', 'score');
casesB = ["bio-flex 10mm", "bio-flex 620kPa", "bio-flex 325kPa"];
n = size(Dg.dub_filt_results,1);
for i = 1:n
    for c = 1:3
        vals = Dg.met(i, 3*(c-1)+(1:3));
        R(end+1,:) = {'BioFlexor candidate', sprintf('pooled front row %d', i), 'minimizeFlx', 'n/a', ...
            'n/a (minimizeFlx)', Dg.dub_filt_results(i,1), Dg.dub_filt_results(i,2), Dg.dub_filt_results(i,3), NaN, ...
            char(casesB(c)), vals(1), vals(2), vals(3), '', '', '', sprintf('score %.3f', Dg.score(i))}; %#ok<AGROW>
    end
end

%% D. Flipped screening (10 candidates, best Xi3 per candidate)
Sc = load('Dig_out/Dig_ExtPinX3_screen_results.mat', 'cand', 'out');
for i = 1:size(Sc.cand,1)
    name = char(Sc.cand{i,1});
    R(end+1,:) = {'ExtensorScreen (flipped signs)', name, 'minimizeExtX3 + minimizeExt', '2trans', ...
        'rib midpoint (at time of run)', Sc.cand{i,2}, Sc.cand{i,3}, Sc.cand{i,4}, Sc.out(i,4), ...
        'pinned pool (best Xi3)', NaN, NaN, NaN, '', '', '', ...
        sprintf('best Xi3 %.2f; pool RMSEratio %.3f FVUratio %.3f; bio RMSE %.3f FVU %.3f MaxR %.3f', ...
        Sc.out(i,4), Sc.out(i,5), Sc.out(i,6), Sc.out(i,7), Sc.out(i,8), Sc.out(i,9))}; %#ok<AGROW>
end

%% Emit
T = cell2table(R, 'VariableNames', {'Category','RunID','Evaluator','Convention','Pbr', ...
    'Xi0_m','Xi1','Xi2','Xi3','Test','RMSE','FVU','MaxResidual', ...
    'RMSEbase','FVUbase','MaxRbase','RMSEratio','FVUratio','Notes'});
writetable(T, 'Dig_results_20260910.xlsx');
writetable(T, 'Dig_results_20260910.csv');
fprintf('Wrote %d rows to Dig_results_20260910.xlsx / .csv\n', height(T));
