function probe_span_20261103()
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
obj = BiPam_pulley2_20261022('probe_span', ptsBi, crossBi, 20, Tbi, ...
    restBi, 0.76*restBi, tendonBi, fitn, 620, xi0, 5e5, 1e4, 6, cfgBi);
npts = size(ptsBi, 1);
% The class's own crossing-2 span (verbatim loop from updateBiPulleyGeometry):
spanC = zeros(N, N);
for ii = 1:N
    for jj = 1:N
        tot = 0;
        for r = obj.PulleyExitIndex(2):npts-1
            pA = obj.pointInFrame(r, 3, ii, jj);
            pB = obj.pointInFrame(r + 1, 3, ii, jj);
            tot = tot + norm(pA - pB);
        end
        spanC(ii, jj) = tot;
    end
end
% My chord (exit row 3 -> row 5 distance, frame-invariant):
chord = zeros(N, N);
for ii = 1:N
    T1 = Tbi(:,:,ii,1);
    for jj = 1:N
        T2 = Tbi(:,:,jj,2);
        pE = mapRowL(ptsBi, crossBi, obj.PulleyExitIndex(2), T1, T2);
        pI = mapRowL(ptsBi, crossBi, npts, T1, T2);
        chord(ii, jj) = norm(pI - pE);
    end
end
fprintf('obj.PulleyExitIndex = [%d %d], Cross = [%d %d]\n', obj.PulleyExitIndex(1), obj.PulleyExitIndex(2), obj.Cross(1), obj.Cross(2))
fprintf('row 21, cols 19:23:\n  spanC : '); fprintf('%9.6f ', spanC(21,19:23)); fprintf('\n')
fprintf('  chord : '); fprintf('%9.6f ', chord(21,19:23)); fprintf('\n')
fprintf('  deltaL: '); fprintf('%9.6f ', obj.deltaL(21,19:23,2)); fprintf('\n')
fprintf('max |spanC - chord| over grid = %.3e m\n', max(abs(spanC(:) - chord(:)), [], 'all'))
fprintf('RefIndex = [%d %d]; ii0=%d iii0=%d\n', obj.RefIndex(1), obj.RefIndex(2), ...
    mod(obj.RefIndex(2)-1, N)+1, floor((obj.RefIndex(2)-1)/N)+1)
% ALSO: the row-frame of row 4 and whether span uses rows Cross(2)..end? print span rows 3..5 segments:
for jj = 21
    for r = obj.PulleyExitIndex(2):npts-1
        pA = obj.pointInFrame(r, 3, 21, jj); pB = obj.pointInFrame(r+1, 3, 21, jj);
        fprintf('  seg %d-%d len = %.6f (rowFrame %d->%d)\n', r, r+1, norm(pA-pB), obj.rowFrame(r), obj.rowFrame(r+1))
    end
end
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
