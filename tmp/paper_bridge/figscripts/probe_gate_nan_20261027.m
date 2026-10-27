% probe_gate_nan_20261027.m -- reproduce the study1 gate (1,1,aware)
% per-pose equilibrium pieces and find which term goes NaN at the fzero
% call (never-executed draft code; the canonical run needs the culprit
% named before any fix is chosen).
function probe_gate_nan_20261027()
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

[Location, ~, info] = buildDistalRingLocation20mm(p1Rec, pEndRec, ...
    tendonRec, ctx);
fprintf('route built: %d poses; tendonMax min %g max %g\n', ...
    size(Location, 3), min(info.tendonMax), max(info.tendonMax));

% routeLengths
N = size(Location, 3);
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
% tendonSpan
span = zeros(N, 1);
for ii = 1:N
    T = ctx.T_Pam(:, :, ii);
    pPrev = RowVecTrans(inv(T), Location(ctx.CrossPoint - 1, :, ii));
    tot = 0;
    for j = ctx.CrossPoint:size(Location, 1) - 1
        tot = tot + norm(pPrev - Location(j, :, ii));
        pPrev = Location(j, :, ii);
    end
    span(ii) = tot;
end
% bracketCb
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
    u = Location(2, :, ii) - Location(1, :, ii);
    un = norm(u);
    ub = RowVecTrans(Thbr \ eye(4), u / un);
    cb(ii) = (ub * CbM) * ub.';
end
fprintf('mL: anynan %d (idx %s)\n', any(isnan(mL)), mat2str(find(isnan(mL))));
fprintf('span: anynan %d (idx %s)\n', any(isnan(span)), mat2str(find(isnan(span))));
fprintf('cb: anynan %d (idx %s)\n', any(isnan(cb)), mat2str(find(isnan(cb))));
d12 = zeros(N, 1);
for ii = 1:N
    d12(ii) = norm(Location(2, :, ii) - Location(1, :, ii));
end
fprintf('row1-2 dist: min %g at pose %d; poses with dist < 1e-9: %s\n', ...
    min(d12), find(d12 == min(d12), 1), mat2str(find(d12 < 1e-9)));

% psi loop (nB=1, G=1, aware)
Fm = maxBPAforce(restRec, '20', Pres);
Ffun = @(s) festo4(Dia, (s / restRec) / KMAX, Pres) * Fm;
fprintf('Ffun(0)=%g Ffun(KMAX*rest/2)=%g Ffun(KMAX*rest)=%g Fm=%g\n', ...
    Ffun(0), Ffun(KMAX * restRec / 2), Ffun(KMAX * restRec), Fm);
s0 = restRec - mL + Xi(1) + tendonRec + 2 * fitting;
DeltaL = max(span) - span;
cSpr = 1 / (wraps * 1.51e-6 * 193e9 / tendonRec);
sHi = KMAX * restRec;
nNaN = 0;
for i = 1:N
    s0ref = s0(i) - DeltaL(i);
    psi = @(qq) qq - cb(i) * Ffun(qq) - Ffun(qq) * cSpr - DeltaL(i) - s0ref;
    p0 = psi(0); pH = psi(sHi);
    if ~isfinite(p0) || ~isfinite(pH)
        nNaN = nNaN + 1;
        if nNaN <= 5
            fprintf(['pose %d: psi(0)=%g psi(sHi)=%g cb=%g mL=%g span=%g ' ...
                'DeltaL=%g s0ref=%g\n'], i, p0, pH, cb(i), mL(i), span(i), ...
                DeltaL(i), s0ref);
        end
    end
end
fprintf('poses with non-finite psi endpoints: %d of %d\n', nNaN, N);
fprintf('PROBE DONE\n');
end
