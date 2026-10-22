% Debug (2026-10-22): reproduce the survey grade() for tib_ant_r at the
% initial sizes and print every size/NaN fact. Throwaway scratch.
scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot'), break; end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(fullfile(root, 'tmp', 'paper_bridge', 'opt_improve_20261022'));

% Load the human curve row for actuator 16
csvPath = fullfile(scriptDir, 'human_torque_curves_20261022.csv');
raw = readtable(csvPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char', 'CommentStyle', '#');
r16 = raw(raw.actuator_id == 16, :);
angCols = raw.Properties.VariableNames(startsWith(raw.Properties.VariableNames, 'theta_deg_'));
tauCols = raw.Properties.VariableNames(startsWith(raw.Properties.VariableNames, 'tau_h_Nm_'));
anglesDeg = cellfun(@(c) r16.(c), angCols);
tauH = cellfun(@(c) r16.(c), tauCols);
fprintf('anglesDeg %s, NaNs %d, tauH %s, NaNs %d\n', mat2str(size(anglesDeg)), ...
    sum(isnan(anglesDeg)), mat2str(size(tauH)), sum(isnan(tauH)));
fprintf('angles: %.2f .. %.2f | tauH: %.2f .. %.2f\n', min(anglesDeg), ...
    max(anglesDeg), min(tauH), max(tauH));

spec0 = biPulleySpecsFromOpenSim({'tib_ant_r'});
spec = spec0(1);
disp(spec.kin)
fprintf('points:\n'); disp(spec.points);
fprintf('cross = [%d %d], frames = %s | %s | %s\n', spec.cross, ...
    spec.frames{:});

% Transforms: joint 1 swept over the CSV grid, joint 2 pinned at 0.
kin = spec.kin;
th1 = anglesDeg(:)' * pi / 180;
th2 = 0;
N1 = numel(th1); N2 = max(numel(th2), 2);   % duplicated pinned slice
T = zeros(4, 4, N1, N2);
for i1 = 1:N1
    T(:, :, i1, 1) = RpToTrans(axisRot(kin.axes(1, :), th1(i1)), kin.pivots(1, :)');
end
for i2 = 1:numel(th2)
    T(:, :, i2, 2) = RpToTrans(eye(3), kin.pivots(2, :)');
end
T(:, :, 2, 2) = T(:, :, 1, 2);
fprintf('T %s, any NaN: %d\n', mat2str(size(T)), any(isnan(T(:))));

axisU = kin.axes(1, :) / norm(kin.axes(1, :));
fprintf('axisU = %s (numel %d)\n', mat2str(axisU, 4), numel(axisU));

points = spec.points;
Lrig = zeros(N1, N2);
for ii = 1:N1
    for iii = 1:N2
        tot = 0;
        for rr = 1:size(points, 1) - 1
            pA = mapRowDbg(points, spec.cross, rr, T(:, :, ii, 1), T(:, :, iii, 2));
            pB = mapRowDbg(points, spec.cross, rr + 1, T(:, :, ii, 1), T(:, :, iii, 2));
            tot = tot + norm(pA - pB);
        end
        Lrig(ii, iii) = tot;
    end
end
Lchar = max(Lrig(:));
x0 = [0.60 * Lchar; 0.05];
fprintf('Lchar %.4f, x0 = [%.4f %.4f]\n', Lchar, x0);

pulley = struct('nPulleyBPA', 1, 'tackleLineParts', [1, 1], 'routingMode', 'moving_via');
obj = BiPam_pulley2_20261022('tib_ant_r', points, spec.cross, 20, T, ...
    x0(1), 0.75 * x0(1), x0(2), 0.025, 620, 0.005, 5e5, 1e4, 6, pulley);
fprintf('Torque_ins %s | strain_p %s | PulleySlack %s\n', ...
    mat2str(size(obj.Torque_ins)), mat2str(size(obj.strain_p)), ...
    mat2str(size(obj.PulleySlack)));
T3 = obj.Torque_ins(:, :, :, 1);
fprintf('T3 %s numel %d; N1*N2 = %d\n', mat2str(size(T3)), numel(T3), N1 * N2);
proj = sum(T3 .* reshape(axisU, 1, 3, 1), 2);
fprintf('proj %s numel %d\n', mat2str(size(proj)), numel(proj));
tauB = reshape(proj, N1, N2);
tauHg = repmat(reshape(tauH, numel(tauH), 1), 1, N2);
fprintf('tauB range %.3f .. %.3f; tauH range %.3f .. %.3f\n', ...
    min(tauB(:)), max(tauB(:)), min(tauHg(:)), max(tauHg(:)));
sgn = sign(tauB) == sign(tauHg);
fprintf('sign agreement: %d of %d cells\n', sum(sgn(:)), numel(sgn));
fprintf('col1 == col2: %d\n', isequal(tauB(:, 1), tauB(:, 2)));
meet = ~obj.PulleySlack & ~obj.PulleyInfeasible & isfinite(tauB) & ...
    sgn & abs(tauB) >= abs(tauHg);
strain = reshape(obj.strain_p, N1, N2);
active = ~obj.PulleySlack & ~obj.PulleyInfeasible & strain > 0 & abs(tauB) > 0;
fprintf('strain range %.4f .. %.4f | slack %d | infeas %d\n', ...
    min(strain(:)), max(strain(:)), sum(obj.PulleySlack(:)), sum(obj.PulleyInfeasible(:)));
fprintf('at x0: meet %d/%d (%.1f%%), active %d/%d (%.1f%%)\n', ...
    sum(meet(:)), numel(meet), 100 * sum(meet(:)) / numel(meet), ...
    sum(active(:)), numel(active), 100 * sum(active(:)) / numel(active));
fprintf('angle  tauB(col1)   tauH  signB signH slack\n');
for i = 1:N1
    fprintf('%7.2f %11.4f %8.3f   %+.0f  %+.0f   %d\n', anglesDeg(i), ...
        tauB(i, 1), tauHg(i, 1), sign(tauB(i, 1)), sign(tauHg(i, 1)), ...
        obj.PulleySlack(i, 1));
end
fprintf('DEBUG DONE\n');

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

function R = axisRot(axis, theta)
a = axis(:) / norm(axis);
Km = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * Km + (1 - cos(theta)) * (Km * Km);
end
