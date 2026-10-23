% Probe (2026-10-22): compare the survey's med_gas graded-grid route
% length against the 2026-09-26 anchor mat (same builder + fold).
root = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));

A = load(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'BiPulley_opensim_med_gas_r_20260926_1655.mat'));
fprintf('ANCHOR med_gas: points (fold applied)\n');
disp(A.points);
if isfield(A, 'Lrig')
    fprintf('anchor Lrig range %.4f .. %.4f (grid %s)\n', min(A.Lrig(:)), ...
        max(A.Lrig(:)), mat2str(size(A.Lrig)));
end
fprintf('anchor xBest = [%.4f %.4f], margin %+.4f\n', A.xBest, A.st.margin);

spec0 = biPulleySpecsFromOpenSim({'med_gas_r'});
fprintf('fresh spec points:\n');
disp(spec0.points);
fprintf('cross = [%d %d], folds{2} = %s\n', spec0.cross, ...
    mat2str(spec0.kin.pointFolds{2}.'));

% My graded-grid Lrig: joint 1 (knee) pinned 0, joint 2 (ankle) swept
% over the 25 CSV angles.
raw = readtable(fullfile(root, 'tmp', 'paper_bridge', ...
    'capability_table_20261022', 'human_torque_curves_20261022.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', ...
    'TextType', 'char', 'CommentStyle', '#');
r13 = raw(raw.actuator_id == 13, :);
angCols = raw.Properties.VariableNames( ...
    startsWith(raw.Properties.VariableNames, 'theta_deg_'));
anglesDeg = cellfun(@(c) r13.(c), angCols);

kin = spec0.kin;
thetaP = anglesDeg(:)' * pi / 180;
nA = numel(thetaP);
N1 = nA; N2 = nA;
T1pin = RpToTrans(axisRotDbg(kin.axes(1, :), 0), kin.pivots(1, :)');
T = zeros(4, 4, N1, N2);
for i1 = 1:N1
    T(:, :, i1, 1) = T1pin;
end
for i2 = 1:N2
    T(:, :, i2, 2) = RpToTrans(axisRotDbg(kin.axes(2, :), thetaP(i2)), ...
        kin.pivots(2, :)');
end

points = spec0.points;
points(spec0.cross(2):end, :) = points(spec0.cross(2):end, :) + ...
    spec0.kin.pointFolds{2}';
cross = spec0.cross;
Lrig = zeros(N1, N2);
for ii = 1:N1
    for iii = 1:N2
        tot = 0;
        for rr = 1:size(points, 1) - 1
            pA = mapRowDbg(points, cross, rr, T(:, :, ii, 1), T(:, :, iii, 2));
            pB = mapRowDbg(points, cross, rr + 1, T(:, :, ii, 1), T(:, :, iii, 2));
            tot = tot + norm(pA - pB);
        end
        Lrig(ii, iii) = tot;
    end
end
fprintf('graded-grid Lrig range %.4f .. %.4f -> Lchar %.4f\n', ...
    min(Lrig(:)), max(Lrig(:)), max(Lrig(:)));
fprintf('bounds: rest [%.4f %.4f], tendon [%.4f %.4f]\n', ...
    0.30 * max(Lrig(:)), 0.90 * max(Lrig(:)), 0.02, max(0.30 * max(Lrig(:)), 0.06));

% What take-up demand does the anchor's own xBest face on THIS grid?
rest = A.xBest(1); tendon = A.xBest(2);
mL = mean(Lrig(:));
s0 = rest - mL + 0.005 + tendon + 2 * 0.025;
fprintf('anchor sizes on this grid: mean mL %.4f -> s0 %.4f vs KMAX*rest %.4f\n', ...
    mL, s0, 0.25 * rest);
fprintf('PROBE DONE\n');

function p = mapRowDbg(points, cross, r, T1, T2)
if r < cross(1)
    f = 1;
elseif cross(1) == cross(2)
    f = 3;
elseif r < cross(2)
    f = 2;
else
    f = 3;
end
p = points(r, :);
while f > 1
    if f == 3
        p = RowVecTrans(T2, p);
        f = 2;
    else
        p = RowVecTrans(T1, p);
        f = 1;
    end
end
end

function R = axisRotDbg(axis, theta)
a = axis(:) / norm(axis);
Km = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * Km + (1 - cos(theta)) * (Km * Km);
end
