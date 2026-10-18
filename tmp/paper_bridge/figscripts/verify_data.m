%verify_data.m — run the four evaluators at the plan coefficients; print GoF + curve
%ranges for comparison against the chapter values of record and the exemplar .fig data.
%READ-ONLY: no files written.
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');
setenv('FLX_T5_YMM', '');   % ensure off
setenv('EXTX3_TRANS', '');  % ensure default (2trans)

fprintf('=== 1) minimizeFlxPin(0.0119, 1.02e5, 1.06e4) tests 1-2 ===\n');
[fF, bpaF] = minimizeFlxPin(0.0119, 1.02e5, 1.06e4, [1 2]);
for i = 1:2
    fprintf('test %d (rest=%.3f): GoF RMSE=%.4f FVU=%.4f Max=%.4f\n', i, bpaF(i).rest, fF(i,1), fF(i,2), fF(i,3));
    fprintf('   M range [%.4g..%.4g] | M_p range [%.4g..%.4g] | Mexp range [%.4g..%.4g] n=%d | M_h range [%.4g..%.4g]\n', ...
        min(bpaF(i).M), max(bpaF(i).M), min(bpaF(i).M_p(:,3)), max(bpaF(i).M_p(:,3)), ...
        min(bpaF(i).Mexp), max(bpaF(i).Mexp), numel(bpaF(i).Mexp), min(bpaF(i).M_h), max(bpaF(i).M_h));
end

% Mold (original old-surface prediction) via Dig_FlxPin_2brkt pattern
load('ForceStrainForFit.mat', 'z');
z(isnan(z)) = 0;
X = linspace(0,620,20); X = X(2:20);
Y = linspace(0,1,30);
for i = 1:2
    Yq = bpaF(i).strain./((bpaF(i).rest - bpaF(i).Kmax)/bpaF(i).rest);
    Xq = bpaF(i).P.*ones(size(Yq));
    Vq = interp2(X, Y', z, Xq, Yq, 'spline');
    Fold = Vq.*bpaF(i).unitD;
    Fq = hypot(Fold(:,1), Fold(:,2));
    Mold = -bpaF(i).mA.*Fq;
    Mold(Yq>1) = 0;
    Mold(bpaF(i).strain < 0) = NaN;
    fprintf('   Mold(test %d) range [%.4g..%.4g] (exemplar tile %d Original range)\n', i, min(Mold), max(Mold), i+2);
end

fprintf('\n=== 2) minimizeExtX3(-0.0109, 1.02e5, 1.06e4, 0.123) tests [9 7 3] ===\n');
[fE, bpaE] = minimizeExtX3(-0.0109, 1.02e5, 1.06e4, 0.123, [9 7 3]);
pmap = [9 7 3];
for k = 1:3
    i = pmap(k);
    fprintf('ke(%d) (rest=%.3f ten=%.3f): GoF RMSE=%.4f FVU=%.4f Max=%.4f | M [%.4g..%.4g] | M_p [%.4g..%.4g] | Mexp [%.4g..%.4g] n=%d | M_h [%.4g..%.4g]\n', ...
        i, bpaE(i).rest, bpaE(i).ten, fE(i,1), fE(i,2), fE(i,3), ...
        min(bpaE(i).M), max(bpaE(i).M), min(bpaE(i).M_p(:,3)), max(bpaE(i).M_p(:,3)), ...
        min(bpaE(i).Mexp), max(bpaE(i).Mexp), numel(bpaE(i).Mexp), min(bpaE(i).M_h), max(bpaE(i).M_h));
end

fprintf('\n=== 3) minimizeExt(-0.0109, 1.02e5, 1.06e4, 0.123, 1) — bio extensor 51.8 cm ===\n');
[fB, bpaB] = minimizeExt(-0.0109, 1.02e5, 1.06e4, 0.123, 1);
fprintf('GoF RMSE=%.4f FVU=%.4f Max=%.4f  (chapter initial: 0.8920 / 0.3396 / 1.4188)\n', fB(1,1), fB(1,2), fB(1,3));
fprintf('rest=%.4f | M [%.4g..%.4g] | M_p [%.4g..%.4g] | Mexp n=%d [%.4g..%.4g] | M_h [%.4g..%.4g] | Ak range [%.4g..%.4g]\n', ...
    bpaB(1).rest, min(bpaB(1).M), max(bpaB(1).M), min(bpaB(1).M_p(:,3)), max(bpaB(1).M_p(:,3)), ...
    numel(bpaB(1).Mexp), min(bpaB(1).Mexp), max(bpaB(1).Mexp), min(bpaB(1).M_h), max(bpaB(1).M_h), ...
    min(bpaB(1).Ak), max(bpaB(1).Ak));

fprintf('\n=== 4) minimizeFlx(0.0119, 1.02e5, 1.06e4) — bio flexor ===\n');
[f10, f20full, f20part, bpaL] = minimizeFlx(0.0119, 1.02e5, 1.06e4);
fprintf('10mm  GoF RMSE=%.4f FVU=%.4f Max=%.4f  (chapter: 0.3007 / 0.0102 / 0.4647)\n', f10(1), f10(2), f10(3));
fprintf('20mm 620 GoF RMSE=%.4f FVU=%.4f Max=%.4f  (chapter: 1.8024 / 0.0499 / 3.5678)\n', f20full(1), f20full(2), f20full(3));
fprintf('20mm part GoF RMSE=%.4f FVU=%.4f Max=%.4f\n', f20part(1), f20part(2), f20part(3));
for j = 1:3
    fprintf('bpaL(%d): rest=%.4f ten=%.3f P=%.1f dia=%.1f | M [%.4g..%.4g] | M_p [%.4g..%.4g] | Mexp n=%d [%.4g..%.4g] x[%.4g..%.4g]\n', ...
        j, bpaL(j).rest, bpaL(j).ten, bpaL(j).P, bpaL(j).dBPA*1000, ...
        min(bpaL(j).M), max(bpaL(j).M), min(bpaL(j).M_p(:,3)), max(bpaL(j).M_p(:,3)), ...
        numel(bpaL(j).Mexp), min(bpaL(j).Mexp), max(bpaL(j).Mexp), min(bpaL(j).Aexp), max(bpaL(j).Aexp));
end
fprintf('(exemplar: panel A Predicted [?..%.4g] Mexp n=10 x[-81..-2]; panel B Predicted620 [-41.89..-4.433] Meas620 n=11 y[-29.9..-5.646]; pink Meas n=4 x[-9,0,3,7.7])\n', 0);

fprintf('\n=== 5) Human line from FlxGrp.fig ===\n');
hG = openfig(fullfile(festodir, 'FlxGrp.fig'), 'reuse');
axList = findall(hG, 'Type', 'axes');
for k = 1:numel(axList)
    kids = axList(k).Children;
    if isempty(kids), continue; end
    for c = 1:numel(kids)
        o = kids(c);
        if isa(o, 'matlab.graphics.chart.primitive.Line') && strcmp(o.DisplayName, 'Human Model')
            hx = o.XData; hy = o.YData;
            fprintf('Human Model n=%d x[%.4g..%.4g] y[%.4g..%.4g]\n', numel(hx), min(hx), max(hx), min(hy), max(hy));
            humanX = hx; humanY = hy; %#ok<NASGU>
        end
    end
end
close(hG);
% annotation strings of FlxPin_group.fig (panel letters / row labels)
hP = openfig(fullfile(festodir, 'FlxPin_group.fig'), 'reuse');
tbs = findall(hP, 'Type', 'textbox');
fprintf('FlxPin_group.fig textboxes:\n');
for k = 1:numel(tbs)
    try
        fprintf('  str=[%s] pos=[%s] rot=%s\n', tbs(k).String, num2str(round(tbs(k).Position*1000)/1000), num2str(tbs(k).Rotation));
    catch, end
end
ttl = get(hP, 'Title'); %#ok<NASGU> % no fig-level title
close(hP);
fprintf('\nDONE_VERIFY\n');
