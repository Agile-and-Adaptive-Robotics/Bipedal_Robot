% probe_gate2_20261027.m -- per-pose diff of the study1 closure port vs
% MonoPam_pulley at the design of record (1,1,aware). Extracts the
% class's cb from PulleyTravel = sCon - s0ref - Fmag.*cb.
function probe_gate2_20261027()
R = 'D:\GitHub\Bipedal_Robot';
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));
C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), 'xFullBest', 'ctxCamBest', ...
    'XiUsed');
xRec = C.xFullBest(:).';
p1Rec = xRec(1:3); pEndRec = xRec(4:6);
restRec = xRec(7); tendonRec = xRec(8);
ctx = C.ctxCamBest;
Xi = C.XiUsed(:).';
Dia = 20; Pres = 620; KMAX = 0.255; fitting = 0.021; wraps = 3;

[Location, ~, ~] = buildDistalRingLocation20mm(p1Rec, pEndRec, ...
    tendonRec, ctx);
N = size(Location, 3);

% ---- class side
pulCfg = struct('nPulleyBPA', 1, 'tackleLineParts', 1, 'gain', 1, ...
    'routingMode', 'moving_via');
obj = MonoPam_pulley('gate', Location, ctx.CrossPoint, Dia, ctx.T_Pam, ...
    restRec, restRec * (1 - KMAX), tendonRec, fitting, Pres, Xi(1), ...
    Xi(2), Xi(3), wraps, pulCfg);
ftClass = double(obj.Ftendon(:));
sClass = double(obj.sContraction(:));
mLClass = double(obj.MuscleLength(:));
dLClass = double(obj.deltaL(:));
slackC = logical(obj.PulleySlack(:));
infeasC = logical(obj.PulleyInfeasible(:));
FmagC = double(obj.FsingleBPA(:));
ptC = double(obj.PulleyTravel(:));
kSprC = double(obj.kSpr);

% ---- my port side
mL = zeros(N, 1);
for ii = 1:N
    tot = 0;
    for j = 1:size(Location, 1) - 1
        pA = Location(j, :, ii);
        pB = Location(j + 1, :, ii);
        if j + 1 == ctx.CrossPoint
            pB = RowVecTrans(ctx.T_Pam(:, :, ii), pB);
        end
        tot = tot + norm(pA - pB);
    end
    mL(ii) = tot;
end
span = zeros(N, 1);
for ii = 1:N
    exitT1 = RowVecTrans(ctx.T_Pam(:, :, ii) \ eye(4), ...
        Location(ctx.CrossPoint - 1, :, ii));
    span(ii) = norm(exitT1 - Location(ctx.CrossPoint, :, ii));
end
iRef = find(span == max(span), 1);
DeltaL = span(iRef) - span;

% cb (computeForceVector + Lok frame)
cb = zeros(N, 1);
Pbr = [9.48, -36.21, 30.86] / 1000;
phbrA = p1Rec - Pbr;
th1 = atan2(phbrA(2), phbrA(1));
RhbrZ = [cos(th1), -sin(th1), 0; sin(th1), cos(th1), 0; 0, 0, 1];
pbrhA = RhbrZ' * phbrA.';
th2 = atan2(pbrhA(3), pbrhA(1));
Ry = [cos(th2), 0, sin(th2); 0, 1, 0; -sin(th2), 0, cos(th2)];
Thbr = RpToTrans(RhbrZ * Ry', Pbr');
CbM = diag(1 ./ [Xi(2), Xi(3), Xi(3)]);
for ii = 1:N
    pt1 = Location(1, :, ii);
    pt2 = pt1;
    for k = 2:size(Location, 1)
        d = norm(Location(k, :, ii) - Location(1, :, ii));
        if d > 1e-6
            if k == ctx.CrossPoint
                pt2 = RowVecTrans(ctx.T_Pam(:, :, ii), Location(k, :, ii));
            else
                pt2 = Location(k, :, ii);
            end
            break
        end
    end
    u = pt2 - pt1;
    if norm(u) >= 1e-12
        u = u / norm(u);
        ub = RowVecTrans(Thbr \ eye(4), u);
        cb(ii) = (ub * CbM) * ub.';
    end
end

% class cb EXACT via the psi identity at the solved root:
% s - cb*F - F*cSpr - DeltaL - s0ref = 0 (G = nB = 1)
s0refRigid = restRec - mLClass(iRef) + Xi(1) + tendonRec + 2 * fitting;
cSprE = 1 / kSprC;
cbC = nan(N, 1);
okc = FmagC > 0 & ~slackC;
cbC(okc) = (sClass(okc) - FmagC(okc) * cSprE - dLClass(okc) - s0refRigid) ...
    ./ FmagC(okc);
fprintf('rigid s0ref class-side %.6f\n', s0refRigid);

s0 = restRec - mL + Xi(1) + tendonRec + 2 * fitting;
s0ref = s0(iRef);
cSpr = 1 / kSprC;
Fm = maxBPAforce(restRec, '20', Pres);
Ffun = @(s) festo4(Dia, (s / restRec) / KMAX, Pres) * Fm;
sHi = KMAX * restRec;
s = zeros(N, 1); F = zeros(N, 1);
slack = false(N, 1); infeas = false(N, 1);
for i = 1:N
    psi = @(qq) qq - cb(i) * Ffun(qq) - Ffun(qq) * cSpr - DeltaL(i) - s0ref;
    if psi(0) > 0
        s(i) = 0; F(i) = Ffun(0); slack(i) = true;
    elseif psi(sHi) < 0
        s(i) = sHi; F(i) = Ffun(sHi); infeas(i) = true;
    else
        s(i) = max(fzero(psi, [0, sHi]), 0);
        F(i) = Ffun(s(i));
    end
end
ftMine = F;

fprintf('mL diff max: %.3g\n', max(abs(mL - mLClass)));
fprintf('DeltaL diff max: %.3g\n', max(abs(DeltaL - dLClass)));
fprintf('iRef: mine %d, class refindex? (span max at %d)\n', iRef, iRef);
fprintf('cb diff max (valid poses): %.3g\n', max(abs(cb(okc) - cbC(okc))));
fprintf('kSpr class %.6g, mine %.6g\n', kSprC, 1 / cSpr);
fprintf('Fm: %.6g (class Fmax %.6g)\n', Fm, double(obj.Fmax));
fprintf('s0ref mine %.6g\n', s0ref);
fprintf('slack: mine %d class %d ; infeas: mine %d class %d\n', ...
    nnz(slack), nnz(slackC), nnz(infeas), nnz(infeasC));
[~, iWorst] = max(abs(ftMine - ftClass) ./ max(1, abs(ftClass)));
fprintf('worst pose %d: mine %.6f class %.6f  (cb m %.6g c %.6g; s m %.6g c %.6g)\n', ...
    iWorst, ftMine(iWorst), ftClass(iWorst), cb(iWorst), cbC(iWorst), ...
    s(iWorst), sClass(iWorst));
relErr = max(abs(ftMine - ftClass) ./ max(1, abs(ftClass)));
fprintf('relErr max %.4g ; poses > 1e-6: %d\n', relErr, nnz(abs(ftMine - ftClass) ./ max(1, abs(ftClass)) > 1e-6));
bad = find(abs(ftMine - ftClass) ./ max(1, abs(ftClass)) > 1e-6);
for q = 1:min(8, numel(bad))
    i = bad(q);
    % which row does the direction scan pick, and the u vector
    kPick = 0; uRow = [0 0 0];
    for k = 2:size(Location, 1)
        d = norm(Location(k, :, i) - Location(1, :, i));
        if d > 1e-6
            kPick = k;
            if k == ctx.CrossPoint
                uRow = RowVecTrans(ctx.T_Pam(:, :, i), Location(k, :, i)) - Location(1, :, i);
            else
                uRow = Location(k, :, i) - Location(1, :, i);
            end
            break
        end
    end
    fprintf(['  pose %d: ft %.4f/%.4f cb %.6g/%.6g s %.6g/%.6g kPick %d ' ...
        'u [%s] d12 %g\n'], i, ftMine(i), ftClass(i), cb(i), cbC(i), ...
        s(i), sClass(i), kPick, num2str(uRow / norm(uRow), 4), ...
        norm(Location(2, :, i) - Location(1, :, i)));
end
fprintf('PROBE2 DONE\n');
end
