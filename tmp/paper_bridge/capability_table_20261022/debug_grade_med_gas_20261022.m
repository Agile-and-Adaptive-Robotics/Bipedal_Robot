% Debug (2026-10-22): med_gas_r slot-2 graded grid -- sizes and T facts.
scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot'), break; end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));

csvPath = fullfile(scriptDir, 'human_torque_curves_20261022.csv');
raw = readtable(csvPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char', 'CommentStyle', '#');
r13 = raw(raw.actuator_id == 13, :);
angCols = raw.Properties.VariableNames(startsWith(raw.Properties.VariableNames, 'theta_deg_'));
anglesDeg = cellfun(@(c) r13.(c), angCols);

spec = biPulleySpecsFromOpenSim({'med_gas_r'});
spec = spec(1);
kin = spec.kin;
fprintf('types {%s,%s} kneeSlot %d cross [%d %d]\n', kin.types{1}, ...
    kin.types{2}, kin.kneeSlot, spec.cross);
fprintf('transTheta size %s, range %.1f..%.1f\n', mat2str(size(kin.transTheta)), ...
    min(kin.transTheta), max(kin.transTheta));

th1 = 0;                              % pinned knee (slot 2 swept)
th2 = anglesDeg(:)' * pi / 180;       % swept ankle
% LAYOUT: joint-2 transforms live at T(:,:,i2,2) (dim3 index i2), so
% dim3 must cover N2; the pinned joint-1 transform is duplicated across
% ALL dim3 slices (the class pairs T1(ii) x T2(iii) over the full
% cross product).
N2 = numel(th2);
N1 = max(numel(th1), N2);
T = zeros(4, 4, N1, N2);
T1pin = RpToTrans(axisRot(kin.axes(1, :), th1), kin.pivots(1, :)');
for i1 = 1:N1
    T(:, :, i1, 1) = T1pin;
end
for i2 = 1:numel(th2)
    T(:, :, i2, 2) = RpToTrans(axisRot(kin.axes(2, :), th2(i2)), kin.pivots(2, :)');
end
fprintf('T %s | T11 det %.6f | anyNaN %d\n', mat2str(size(T)), ...
    det(T(:, :, 1, 1)), any(isnan(T(:))));

axisU = kin.axes(2, :) / norm(kin.axes(2, :));
pulley = struct('nPulleyBPA', 1, 'tackleLineParts', [1, 1], 'routingMode', 'moving_via');
x0 = [0.28; 0.05];
tic;
for rep = 1:5
    obj = BiPam_pulley2_20261022('med_gas_r', spec.points, spec.cross, 20, T, ...
        x0(1), 0.75 * x0(1), x0(2), 0.025, 620, 0.005, 5e5, 1e4, 6, pulley);
end
tel = toc / 5;
fprintf('per-eval time (625-cell grid): %.2f s -> 300 evals ~ %.1f min\n', ...
    tel, 300 * tel / 60);
[N1c, N2c] = obj.gridSize();
fprintf('class grid %d x %d | Torque_ins %s | strain_p %s\n', N1c, N2c, ...
    mat2str(size(obj.Torque_ins)), mat2str(size(obj.strain_p)));
T3 = obj.Torque_ins(:, :, :, 2);
proj = sum(T3 .* reshape(axisU, 1, 3, 1), 2);
tauB = reshape(proj, N1, N2);
fprintf('tauB %s range %.3f .. %.3f\n', mat2str(size(tauB)), ...
    min(tauB(:)), max(tauB(:)));
% column-constancy: the graded 25-point curve is any row
fprintf('rows all equal: %d\n', isequal(tauB(1, :), tauB(25, :)));
strain = reshape(obj.strain_p, N1, N2);
fprintf('strain %.4f..%.4f | slack %d | infeas %d\n', min(strain(:)), ...
    max(strain(:)), sum(obj.PulleySlack(:)), sum(obj.PulleyInfeasible(:)));
fprintf('DEBUG DONE\n');

function R = axisRot(axis, theta)
a = axis(:) / norm(axis);
Km = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * Km + (1 - cos(theta)) * (Km * Km);
end
