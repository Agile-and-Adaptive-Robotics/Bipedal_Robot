%Dig_ExtPin_passMatch.m — match step between extensor CV passes.
%Pass 1 (wide Xi1/Xi2): find front candidates near the previous solution
%(Xi0 ~ -10 mm, Xi3 ~ 0.62) and report the Xi1/Xi2 they carry -> pass-2 lock.
%Also: Xi3 solution-space view of the pass-1 front.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

S = load('minimizeExt10mmX3_front_pass1wide.mat', 'filtered_results', 'xCols', 'sol_actual');
F = S.filtered_results; xC = S.xCols; N = size(F,1);
sa = S.sol_actual;
fprintf('Pass 1 pick=1: Xi0=%.4f Xi1=%.3g Xi2=%.3g Xi3=%.3f\n', sa(1), sa(2), sa(3), sa(4));

%dedup (±1mm Xi0, ±5% Xi1/Xi2, ±0.02 Xi3)
dub = zeros(0,4); dubRow = zeros(0,1);
for i = 1:N
    c = F(i,xC);
    dup = false;
    for j = 1:size(dub,1)
        k = dub(j,:);
        if abs(c(1)-k(1)) < 1e-3 && abs(c(2)-k(2)) < 0.05*max(abs(k(2)),eps) && ...
           abs(c(3)-k(3)) < 0.05*max(abs(k(3)),eps) && abs(c(4)-k(4)) < 0.02
            dup = true; break;
        end
    end
    if dup == false
        dub(end+1,:) = c; dubRow(end+1) = i; %#ok<AGROW>
    end
end
nd = size(dub,1);
fprintf('Dedup: %d -> %d unique\n', N, nd);

%previous solution to match
prev = [-0.0101, 0.621];
fprintf('\n-- Candidates near previous pick (Xi0 ~ -10 mm, Xi3 ~ 0.62), best pinned val first --\n');
fprintf('%-6s %8s %9s %9s %7s %8s %8s\n', 'front','Xi0(m)','Xi1','Xi2','Xi3','valRMSE','dist');
dval = sqrt(((dub(:,1) - prev(1))/0.01).^2 + ((dub(:,4) - prev(2))/0.62).^2);
[~, ord] = sort(dval, 'ascend');
lockCands = [];
for q = 1:min(10, nd)
    i = ord(q);
    fprintf('%-6d %8.4f %9.3g %9.3g %7.3f %8.3f %8.3f\n', ...
        dubRow(i), dub(i,1), dub(i,2), dub(i,3), dub(i,4), ...
        F(dubRow(i), xC(end-2)), F(dubRow(i), xC(end))); %#ok<NODEF>
    lockCands(end+1,:) = [dub(i,2), dub(i,3)]; %#ok<AGROW>
end

%Xi3 solution-space view: Xi3 vs pinned validation RMSE across the whole front
fprintf('\n-- Xi3 profile (all unique candidates, sorted by Xi3) --\n');
[~, o3] = sort(dub(:,4), 'ascend');
fprintf('%7s %9s %9s %7s %8s\n', 'Xi3','Xi0(m)','Xi1','Xi2','valRMSE');
for q = 1:nd
    i = o3(q);
    fprintf('%7.3f %9.4f %9.3g %9.3g %8.3f\n', dub(i,4), dub(i,1), dub(i,2), dub(i,3), ...
        F(dubRow(i), xC(end-2)));
end
fprintf('\nSuggested pass-2 lock (best match): XI1=%.4g XI2=%.4g\n', lockCands(1,1), lockCands(1,2));
