function probe_chord_20261103()
%PROBE_CHORD_20261103 Localize the chord-FD residual for the bi fixture:
% compare at sample cells (1) the object's Torque_ins, (2) the analytic
% triple product from the object's OWN deformed arrays, (3) the chord FD,
% (4) the analytic triple product from rigid geometry.

repo = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(repo, 'Code', 'Matlab')));
addpath(fullfile(repo, 'Code', 'Matlab', 'Mesh_Optimization'), '-end');
addpath(fullfile(repo, 'Testing_Data', '2022_02_Festo'), '-end');

a1 = [0 0 1]; p1 = [0.01 0.01 0];
a2 = [0 1 0]; p2 = [0.005 0.02 0];
N = 41;
th1 = linspace(0, 0.35, N);
th2 = linspace(0, 0.35, N);
ptsBi = [0.05 0.04 0; 0.28 0.02 0; 0.30 0.10 0; 0.02 0.20 0; 0.01 0.21 0];
crossBi = [3 4];
Tbi = zeros(4, 4, N, N);
for i1 = 1:N
    Tbi(:, :, i1, 1) = RpToTrans(axisRotM(a1, th1(i1)), p1');
end
for i2 = 1:N
    Tbi(:, :, i2, 2) = RpToTrans(axisRotM(a2, th2(i2)), p2');
end
Lrig = rigidLen(ptsBi, crossBi, Tbi, N, N);
tendonBi = 0.05; fitn = 0.025; xi0 = 0.005;
restBi = (0.75 * max(Lrig(:)) + 0.25 * min(Lrig(:)) - xi0 - tendonBi ...
    - 2 * fitn) / (1 - 0.12);
cfgBi = struct('nPulleyBPA', 3, 'tackleLineParts', 2, 'gain', [1 2], ...
    'routingMode', 'moving_via');
obj = BiPam_pulley2_20261022('probe_bi', ptsBi, crossBi, 20, Tbi, ...
    restBi, 0.76 * restBi, tendonBi, fitn, 620, xi0, 5e5, 1e4, 6, cfgBi);

% Rigid chord grids.
npts = size(ptsBi, 1);
chord = zeros(N, N, 2);
for kk = 1:2
    for ii = 1:N
        T1 = Tbi(:, :, ii, 1);
        for jj = 1:N
            T2 = Tbi(:, :, jj, 2);
            pE = mapRowL(ptsBi, crossBi, obj.PulleyExitIndex(kk), T1, T2);
            pI = mapRowL(ptsBi, crossBi, npts, T1, T2);
            chord(ii, jj, kk) = norm(pI - pE);
        end
    end
end

tauY = reshape(obj.Torque_ins(:, 2, :, 2), N, N);
ut2 = squeeze(obj.u_t(:, :, :, 2));        % N x 3 x N
Tl2 = reshape(obj.segTension(:, 2), N, N);
cells = [21 21; 21 5; 21 35; 5 21; 35 21];
h2 = th2(2) - th2(1);
fprintf('cell (i1,i2) | tau_obj | tau_deformed-arrays | tau_rigid-arrays | -T*FDchord | chord scale\n')
for c = 1:size(cells, 1)
    ii = cells(c, 1); jj = cells(c, 2);
    % (2) analytic triple from the object's own deformed arrays:
    pIdef = reshape(obj.L_p(npts, :, ii, jj), 1, 3);   % slot-3 row: frame 3 coords as stored
    % L_p is in the ROW'S OWN frame: row 5 -> frame 3. map to frame-3 view:
    pI3 = obj.mapVectorToFrame(pIdef, 3, ii, jj, obj.rowFrame(npts));
    u3 = reshape(ut2(ii, :, jj), 1, 3);
    tauDef = Tl2(ii, jj) * cross(pI3, u3);
    % (4) rigid-geometry triple:
    pI3r = obj.mapVectorToFrame(ptsBi(npts, :), 3, ii, jj, 3);
    pE3r = obj.mapVectorToFrame(ptsBi(obj.PulleyExitIndex(2), :), 3, ii, jj, 2);
    u3r = (pE3r - pI3r) / norm(pE3r - pI3r);
    tauRig = Tl2(ii, jj) * cross(pI3r, u3r);
    % (3) chord FD:
    fdCh = -(chord(ii, jj + 1, 2) - chord(ii, jj - 1, 2)) / (2 * h2) * Tl2(ii, jj);
    fprintf('(%2d,%2d) | %9.4f | %9.4f | %9.4f | %9.4f | chord %.4f m\n', ...
        ii, jj, tauY(ii, jj), tauDef(2), tauRig(2), fdCh, chord(ii, jj, 2))
end

% Max |tau_y| over taut cells and the equivalent moment-arm error of the FD.
h2f = th2(2) - th2(1);
dC2 = (chord(:, 3:N, 2) - chord(:, 1:N - 2, 2)) / (2 * h2f);   % centers i2 = 2:N-1
exp2 = -Tl2(:, 2:N - 1) .* dC2;
tauRef2 = tauY(:, 2:N - 1);
mask = ~logical(obj.PulleySlack(:, :)); mask = mask(:, 2:N - 1);
eAbs = max(abs(tauRef2(mask) - exp2(mask)));
fprintf('\nmax |tau_y| = %.4f N*m; FD eAbs = %.4e N*m; equiv arm err = %.3e m; max Tline2 = %.1f N\n', ...
    max(abs(tauY(:))), eAbs, eAbs / max(Tl2(:)), max(Tl2(:)))

% ---- Crossing 1, joint 1: same 4-way comparison + argmax cell --------
tauZ = reshape(obj.Torque_ins(:, 3, :, 1), N, N);
Tl1 = reshape(obj.segTension(:, 1), N, N);
h1 = th1(2) - th1(1);
dC1 = (chord(3:N, :, 1) - chord(1:N - 2, :, 1)) / (2 * h1);    % centers i1 = 2:N-1
exp1 = -Tl1(2:N - 1, :) .* dC1;
tauRef1 = tauZ(2:N - 1, :);
mask1 = ~logical(obj.PulleySlack(:, :)); mask1 = mask1(2:N - 1, :);
res1 = abs(tauRef1(mask1) - exp1(mask1));
eAbs1 = max(res1);
[~, ix] = max(res1);
[irow, icol] = find(mask1);
ii1 = irow(ix) + 1; jj1 = icol(ix);
fprintf('crossing1: max |tau_z| = %.4f N*m; FD eAbs = %.4e N*m; equiv arm err = %.3e m; max Tline1 = %.1f N\n', ...
    max(abs(tauZ(:))), eAbs1, eAbs1 / max(Tl1(:)), max(Tl1(:)))
fprintf('crossing1 argmax cell (i1=%d, i2=%d):\n', ii1, jj1)
pI3 = obj.mapVectorToFrame(ptsBi(npts, :), 2, ii1, jj1, 3);
pE3 = obj.mapVectorToFrame(ptsBi(obj.PulleyExitIndex(1), :), 2, ii1, jj1, 1);
u3 = (pE3 - pI3) / norm(pE3 - pI3);
tauRig1 = Tl1(ii1, jj1) * cross(pI3, u3);
fdCh1 = -(chord(ii1 + 1, jj1, 1) - chord(ii1 - 1, jj1, 1)) / (2 * h1) * Tl1(ii1, jj1);
fprintf('  tau_obj = %.6f | tau_rigid-triple = %.6f | -T*FDchord = %.6f\n', ...
    tauZ(ii1, jj1), tauRig1(3), fdCh1)
fprintf('  5-point stencil check at that cell: ');
fd5 = (chord(ii1 - 2, jj1, 1) - 8 * chord(ii1 - 1, jj1, 1) + 8 * chord(ii1 + 1, jj1, 1) ...
    - chord(ii1 + 2, jj1, 1)) / (12 * h1) * Tl1(ii1, jj1);
fprintf('tau_obj = %.6f | -T*FD5chord = %.6f\n', tauZ(ii1, jj1), -fd5)
end

function L = rigidLen(points, cross, T, N1, N2)
nPts = size(points, 1);
L = zeros(N1, N2);
for ii = 1:N1
    T1 = T(:, :, ii, 1);
    for jj = 1:N2
        T2 = T(:, :, jj, 2);
        tot = 0;
        for r = 1:nPts - 1
            pA = mapRowL(points, cross, r, T1, T2);
            pB = mapRowL(points, cross, r + 1, T1, T2);
            tot = tot + norm(pA - pB);
        end
        L(ii, jj) = tot;
    end
end
end

function p = mapRowL(points, cross, r, T1, T2)
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

function R = axisRotM(axis, theta)
a = axis(:) / norm(axis);
K = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);
end
