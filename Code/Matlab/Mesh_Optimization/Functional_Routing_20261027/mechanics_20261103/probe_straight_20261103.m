function probe_straight_20261103()
repo = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(repo, 'Code', 'Matlab')));
addpath(fullfile(repo, 'Code', 'Matlab', 'Mesh_Optimization'), '-end');
addpath(fullfile(repo, 'Testing_Data', '2022_02_Festo'), '-end');
a1 = [0 0 1]; p1 = [0.01 0.01 0]; a2 = [0 1 0]; p2 = [0.005 0.02 0];
N = 41; th1 = linspace(0, 0.35, N); th2 = linspace(0, 0.35, N);
ptsBi = [0.05 0.04 0; 0.28 0.02 0; 0.30 0.10 0; 0.02 0.20 0; 0.01 0.21 0];
crossBi = [3 4];
Tbi = zeros(4, 4, N, N);
for i1 = 1:N, Tbi(:,:,i1,1) = RpToTrans(axisRotM(a1, th1(i1)), p1'); end
for i2 = 1:N, Tbi(:,:,i2,2) = RpToTrans(axisRotM(a2, th2(i2)), p2'); end
Lrig = rigidLen(ptsBi, crossBi, Tbi, N, N);
tendonBi = 0.05; fitn = 0.025; xi0 = 0.005;
restBi = (0.75*max(Lrig(:)) + 0.25*min(Lrig(:)) - xi0 - tendonBi - 2*fitn) / (1 - 0.12);
cfgBi = struct('nPulleyBPA', 3, 'tackleLineParts', 2, 'gain', [1 2], 'routingMode', 'moving_via');
ptsBiS = ptsBi; ptsBiS(4,:) = ptsBiS(3,:) + 0.6*(ptsBiS(5,:) - ptsBiS(3,:));
objS = BiPam_pulley2_20261022('probe_straight', ptsBiS, crossBi, 20, Tbi, ...
    restBi, 0.76*restBi, tendonBi, fitn, 620, xi0, 5e5, 1e4, 6, cfgBi);
npts = size(ptsBiS, 1);
chord = zeros(N, N);
for ii = 1:N
    T1 = Tbi(:,:,ii,1);
    for jj = 1:N
        T2 = Tbi(:,:,jj,2);
        pE = mapRowL(ptsBiS, crossBi, objS.PulleyExitIndex(2), T1, T2);
        pI = mapRowL(ptsBiS, crossBi, npts, T1, T2);
        chord(ii, jj) = norm(pI - pE);
    end
end
h2 = th2(2) - th2(1);
dC = (chord(:, 3:N) - chord(:, 1:N-2)) / (2*h2);
dL = objS.deltaL(:, :, 2);
dDL = (dL(:, 3:N) - dL(:, 1:N-2)) / (2*h2);
fprintf('max |dDL + dC| = %.3e m/rad (should be ~0 if polyline == chord)\n', max(max(abs(dDL + dC))))
fprintf('max |dC| = %.3e, max |dDL| = %.3e m/rad\n', max(abs(dC(:))), max(abs(dDL(:))))
tauY = reshape(objS.Torque_ins(:, 2, :, 2), N, N);
Tl2 = reshape(objS.segTension(:, 2), N, N);
res = tauY(:, 2:N-1) - Tl2(:, 2:N-1) .* dDL;
fprintf('polyline-form residual max = %.4e N*m; chord-form max = %.4e N*m\n', ...
    max(abs(res(:))), max(abs(tauY(:, 2:N-1) + Tl2(:, 2:N-1).*dC), [], 'all'))
fprintf('sample cell (21,21): tau = %.6f | T*dDL = %.6f | -T*dC = %.6f\n', ...
    tauY(21,21), Tl2(21,21)*dDL(20,20), -Tl2(21,21)*dC(20,20))
end
function L = rigidLen(points, cross, T, N1, N2)
nPts = size(points, 1); L = zeros(N1, N2);
for ii = 1:N1
    T1 = T(:,:,ii,1);
    for jj = 1:N2
        T2 = T(:,:,jj,2); tot = 0;
        for r = 1:nPts-1
            pA = mapRowL(points, cross, r, T1, T2);
            pB = mapRowL(points, cross, r+1, T1, T2);
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
    if f == 3, p = RowVecTrans(T2, p); f = 2;
    else, p = RowVecTrans(T1, p); f = 1; end
end
end
function R = axisRotM(axis, theta)
a = axis(:)/norm(axis);
K = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta)*K + (1 - cos(theta))*(K*K);
end
