function vasmed_rp0_probe_20261028()
% Numerical probe: WHERE do MonoPam_pulley (1,1) and
% MonoPamDataExplicit_balanceX3 disagree on the distal-ring route?
% Checks the closure identity s0(i) == s0ref + DeltaL(i) (mL = C + span).

R = 'D:\GitHub\Bipedal_Robot';
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));

C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), 'xFullBest', 'ctxCamBest', 'XiUsed');
xRec = C.xFullBest(:).';
p1 = xRec(1:3); pEnd = xRec(4:6); rest = xRec(7); tendon = xRec(8);
ctx = C.ctxCamBest; Xi = C.XiUsed(:).';
Dia = 20; Pres = 620; KMAX = 0.255; fitting = 0.021; wraps = 3;
phiD = ctx.phiD(:);

[Location, bendM, ~] = buildDistalRingLocation20mm(p1, pEnd, tendon, ctx);

% --- mL from the segment formula (both classes' SegmentLengths) ---
mL = zeros(numel(phiD), 1);
for ii = 1:numel(phiD)
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

% --- pulley class internals ---
pulCfg = struct('nPulleyBPA', 1, 'tackleLineParts', 1, 'gain', 1, ...
    'routingMode', 'moving_via');
obj = MonoPam_pulley('probe', Location, ctx.CrossPoint, Dia, ctx.T_Pam, ...
    rest, rest * (1 - KMAX), tendon, fitting, Pres, Xi(1), Xi(2), Xi(3), ...
    wraps, pulCfg);

mLclass = obj.MuscleLength;
span = zeros(numel(phiD), 1);
Erow = ctx.CrossPoint - 1;
for ii = 1:numel(phiD)
    exitT1 = RowVecTrans(ctx.T_Pam(:, :, ii) \ eye(4), Location(Erow, :, ii));
    span(ii) = norm(exitT1 - Location(ctx.CrossPoint, :, ii));
end
iRef = obj.RefIndex;
DeltaL_span = span(iRef) - span;
s0 = rest - mL + Xi(1) + tendon + 2 * fitting;
s0ref = s0(iRef);
identityErr = s0 - (s0ref + DeltaL_span);

fprintf('mL formula vs class: max abs diff %.3g m\n', max(abs(mL - mLclass)));
fprintf('RefIndex %d (%.1f deg). s0ref %.5f\n', iRef, phiD(iRef), s0ref);
fprintf('identity |s0 - (s0ref+DeltaL)|: max %.5f m at pose %d (%.1f deg)\n', ...
    max(abs(identityErr)), find(abs(identityErr) == max(abs(identityErr)), 1), ...
    phiD(find(abs(identityErr) == max(abs(identityErr)), 1)));
fprintf('mL(i)-mL(ref) vs DeltaL_span: corr %.6f, max diff %.5f\n', ...
    corr((mL - mL(iRef))', DeltaL_span'), max(abs((mL - mL(iRef)) + DeltaL_span)));

% --- X3 class state at the same poses ---
kmaxL = rest * (1 - KMAX);
x3 = MonoPamDataExplicit_balanceX3('p', Location, ctx.CrossPoint, Dia, ...
    ctx.T_Pam, rest, kmaxL, tendon, fitting, Pres, Xi(1), Xi(2), Xi(3), ...
    0, wraps, phiD, 1, bendM);
Fx3 = vecnorm(double(x3.F_p), 2, 2);
Fpul = double(obj.FsingleBPA(:));
sPul = double(obj.sContraction(:));

% X3 solved contraction: strain_f * rest
sX3 = double(x3.strain_f) * rest;
fprintf('\npose  angle   mL      span    s0      s0ref+DL  sPul    sX3    Fpul     Fx3\n');
for ip = [1, 4, 7, 20, 50, 78, 92, 100]
    fprintf('%3d %7.1f %+.4f %+.4f %+.4f %+.4f  %+.4f %+.4f %8.1f %8.1f\n', ...
        ip, phiD(ip), mL(ip), span(ip), s0(ip), s0ref + DeltaL_span(ip), ...
        sPul(ip), sX3(ip), Fpul(ip), Fx3(ip));
end
fprintf('strain_f range: %.4f..%.4f ; pulley sCon/rest range: %.4f..%.4f\n', ...
    min(double(x3.strain_f)), max(double(x3.strain_f)), ...
    min(sPul) / rest, max(sPul) / rest);
fprintf('slack poses (pulley): %d ; x3 r==0 poses: %d\n', nnz(obj.PulleySlack), ...
    nnz(Fx3 >= maxBPAforce(rest, '20') * 0.999));
end
