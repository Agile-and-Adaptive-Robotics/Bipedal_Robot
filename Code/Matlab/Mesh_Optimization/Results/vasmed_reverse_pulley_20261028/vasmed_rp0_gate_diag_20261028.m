function vasmed_rp0_gate_diag_20261028()
% vasmed_rp0_gate_diag_20261028 -- gate diagnostics (2026-10-28)
%
% Re-runs the (1,1) equivalence gates with per-pose diagnostics and a
% force-scaled metric: the rp1 gate used max(1,|F|) as the denominator,
% which explodes at near-zero-force poses (the blocked-force endpoint of
% the sweep) and reports ~1.0 "errors" where the absolute force error is
% millinewtons. Here: relative error on poses with |F| >= 5% of Fmax,
% absolute error everywhere, worst-pose indices, and the Spr-convention
% audit (X3 mult*wraps = 18*Aeff*E/L vs pulley wraps = 3*Aeff*E/L).

R = 'D:\GitHub\Bipedal_Robot';
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));
here = fileparts(mfilename('fullpath'));

C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), 'xFullBest', 'ctxCamBest', 'XiUsed');
xRec = C.xFullBest(:).';
p1 = xRec(1:3); pEnd = xRec(4:6); rest = xRec(7); tendon = xRec(8);
ctx = C.ctxCamBest;
Xi = C.XiUsed(:).';
Dia = 20; Pres = 620; KMAX = 0.255; fitting = 0.021; wraps = 3;
phiD = ctx.phiD(:);

[Location, bendM, ~] = buildDistalRingLocation20mm(p1, pEnd, tendon, ctx);
mL = zeros(numel(phiD), 1);
for ii = 1:numel(mL)
    tot = 0;
    for j = 1:size(Location, 1) - 1
        pA = Location(j, :, ii); pB = Location(j + 1, :, ii);
        if j + 1 == ctx.CrossPoint
            pB = RowVecTrans(ctx.T_Pam(:, :, ii), pB);
        end
        tot = tot + norm(pA - pB);
    end
    mL(ii) = tot;
end

Fmax = maxBPAforce(rest, '20');
kmaxL = rest * (1 - KMAX);
pulCfg = struct('nPulleyBPA', 1, 'tackleLineParts', 1, 'gain', 1, ...
    'routingMode', 'moving_via');
obj = MonoPam_pulley('diag', Location, ctx.CrossPoint, Dia, ctx.T_Pam, ...
    rest, kmaxL, tendon, fitting, Pres, Xi(1), Xi(2), Xi(3), wraps, pulCfg);
Fpul = double(obj.FsingleBPA(:));
sPul = double(obj.sContraction(:));

x3_0 = MonoPamDataExplicit_balanceX3('d0', Location, ctx.CrossPoint, Dia, ...
    ctx.T_Pam, rest, kmaxL, tendon, fitting, Pres, Xi(1), Xi(2), Xi(3), ...
    0, wraps, phiD, 1, bendM);
Fx3_0 = vecnorm(double(x3_0.F_p), 2, 2);

x3_1 = MonoPamDataExplicit_balanceX3('d1', Location, ctx.CrossPoint, Dia, ...
    ctx.T_Pam, rest, kmaxL, tendon, fitting, Pres, Xi(1), Xi(2), Xi(3), ...
    Xi(4), wraps, phiD, 1, bendM);
Fx3_1 = vecnorm(double(x3_1.F_p), 2, 2);

Lm0 = mL - tendon - 2 * fitting - Xi(1);
strain0 = (rest - Lm0) / rest;
comp = max(0, 1 - strain0 / KMAX);
deltaL = Xi(4) .* bendM(:) .* comp .^ 2;
Fcorr = max(festo4(Dia, (sPul + deltaL) / rest / KMAX, Pres) * Fmax, 0);

fprintf('Fmax(rest=%.4f) = %.1f N\n', rest, Fmax);
fprintf('kSpr: pulley(3) = %.4g N/m ; X3(6*3=18) = %.4g N/m (Aeff*E/L = %.4g)\n', ...
    3 * 1.51e-6 * 193e9 / tendon, 18 * 1.51e-6 * 193e9 / tendon, ...
    1.51e-6 * 193e9 / tendon);

reportGate('GATE A  Xi3=0 : pulley vs X3', Fpul, Fx3_0, Fmax, phiD);
reportGate('GATE B  Xi3=ad : posthoc vs X3', Fcorr, Fx3_1, Fmax, phiD);

% also X3 Xi3=0 vs pulley torque-z
tzP = double(obj.Torque_p(:, 3));
tzX = double(x3_0.Torque_p(:, 3));
m = abs(tzP - tzX);
fprintf('GATE A torque-z: abs max %.4g N*m at pose %d (%.1f deg); rel(on |tz|>2%%max) %.3g\n', ...
    max(m), find(m == max(m), 1), phiD(find(m == max(m), 1)), ...
    max(abs(tzP(abs(tzX) > 0.02 * max(abs(tzX), [], 'omitnan')) - ...
    tzX(abs(tzX) > 0.02 * max(abs(tzX), [], 'omitnan')))) / ...
    max(abs(tzX(abs(tzX) > 0.02 * max(abs(tzX), [], 'omitnan')))));

save(fullfile(here, 'vasmed_rp0_gatediag_20261028.mat'), 'Fpul', 'Fx3_0', ...
    'Fx3_1', 'Fcorr', 'deltaL', 'sPul', 'phiD', 'Fmax');
end

function reportGate(tag, Fa, Fb, Fmax, phiD)
d = abs(Fa - Fb);
big = abs(Fb) >= 0.05 * Fmax;
relBig = max(d(big) ./ abs(Fb(big)));
[~, iW] = max(d);
fprintf(['%s: poses |F|>=5%%Fmax: %d/%d, max rel err %.3g ; ' ...
    'abs err: max %.4g N (pose %d, %.1f deg, Fpul %.2f Fx3 %.2f)\n'], ...
    tag, nnz(big), numel(Fb), relBig, max(d), iW, phiD(iW), Fa(iW), Fb(iW));
end
