function out = fr_identities(C, cs, kin, ~, ~)
%FR_IDENTITIES Mechanics identities on a case's seed evaluation (gate
% check c): both bracket contributions, parallel load sharing,
% transmission work consistency. All identities are the ones verified by
% mechanics_20261103/mech4 for BiPam_pulley2_20261022; they are re-run
% HERE against this campaign's wrapper.
out = struct('passed', false);

qg = struct('passIdx', 1, 'qList', {cs.passes{1}.qStructs});
[mech, meta] = fr_mechBuild(C, cs, fr_unpackX(cs, cs.xSeed), kin, ...
    'compliant', qg);

% (1) both bracket contributions: the row-1 Xi1/Xi2 bracket deflects AND
% the tackle-mount/Bowden-anchor/tendon-stretch channel is present, and
% each changes the solution vs the rigid branch.
[mechR, ~] = fr_mechBuild(C, cs, fr_unpackX(cs, cs.xSeed), kin, ...
    'rigid', qg);
Ldefl = zeros(mech.nG, 1);
for i = 1:mech.nG
    Ldefl(i) = norm(mech.cell.L_p(1, :, i, 1) - meta.loc(1, :));
end
bracketDeflMax = max(Ldefl);
gamaMax = max(abs(mech.gama));
tauDiff = max(abs(mech.strain_p - mechR.strain_p));
bracketOk = bracketDeflMax > 1e-9 && tauDiff > 1e-9;
lineOk = gamaMax > 1e-12;
fprintf(['    both bracket contributions (%s): row-1 bracket ' 'deflection max %.4f mm, tendon stretch gama max %.4f mm, ' 'dstrain vs rigid %.3e -> %s\n'], cs.name, bracketDeflMax * 1000, gamaMax * 1000, tauDiff, fr_tf(bracketOk && lineOk));

% (2) parallel load sharing: Ftendon == nBPA * Fsingle / G at taut cells
kAct = mech.ActiveCrossing;
nBPA = cs.txConfig.nPulleyBPA(kAct);
G = cs.txConfig.gain(kAct);
taut = ~mech.PulleySlack & ~mech.PulleyInfeasible;
if any(taut)
    res = abs(mech.Ftendon(taut) - nBPA * mech.FsingleBPA(taut) / G);
    relRes = max(res) / max(max(abs(mech.Ftendon(taut))), eps);
else
    relRes = inf;
end
shareOk = relRes < 1e-9;
fprintf(['    parallel load sharing (%s): max |Ftendon - nBPA*F/G| ' '/ max|Ftendon| = %.3e -> %s\n'], cs.name, relRes, fr_tf(shareOk));
% (3) transmission work consistency: nBPA*F*travel == Ft*(deltaL+gama)
% verified on the RIGID build (compliance off: cb = 0, exactly the form
% mechanics_20261103/mech4 [Bi d] verified at 6e-14 J). The COMPLIANT
% build is verified through the class's own closure-with-compliance
% (mech4 [Bi f]): G*travel - (deltaL + gama) + G*cbF == 0 via the
% DeflClosure coefficients.
mechR = fr_mechBuild(C, cs, fr_unpackX(cs, cs.xSeed), kin, 'rigid', ...
    struct('passIdx', 1, 'qList', {cs.passes{1}.qStructs}));
kR = mechR.ActiveCrossing;
taut2 = ~mechR.PulleySlack & ~mechR.PulleyInfeasible & ...
    abs(mechR.PulleyTravel) > 1e-12;
if any(taut2)
    wL = nBPA * mechR.FsingleBPA(taut2) .* mechR.PulleyTravel(taut2);
    wR = mechR.Ftendon(taut2) .* (mechR.deltaL(taut2, kR) + ...
        mechR.gama(taut2));
    scale = max(abs(wL));
    workRes = max(abs(wL - wR)) / max(scale, eps);
else
    workRes = inf;
end
workOk = workRes < 1e-8;
fprintf(['    transmission work, rigid (%s): max |nBPA*F*travel - ' ...
    'Ft*(deltaL+gama)|/scale = %.3e -> %s\n'], cs.name, workRes, ...
    fr_tf(workOk));
% compliant closure WITH the deflection terms: the class's own solve is
% G*travel - deltaL - gama + cAdd*F == 0 (cAdd = the class's
% DeflClosure_cAdd channel, the mount/bowden absorption coefficients;
% mech4 [Bi f] verified the class against its own coefficients at 1.2e-16)
cAdd = mech.cell.DeflClosure_cAdd(1:mech.nG);
cloRes = max(abs(mech.cfg.gain(kAct) * mech.PulleyTravel(taut) - ...
    mech.deltaL(taut, kAct) - mech.gama(taut) + ...
    cAdd(taut) .* mech.FsingleBPA(taut))) / ...
    max(max(abs(mech.deltaL(:, kAct))) + eps, 1e-9);
cloOk = cloRes < 1e-8;
fprintf(['    compliant closure w/ compliance terms (%s): max |G*travel ' ...
    '- deltaL - gama + cAdd*F|/dLscale = %.3e -> %s\n'], cs.name, ...
    cloRes, fr_tf(cloOk));

out.passed = bracketOk && lineOk && shareOk && workOk && cloOk;
out.bracketDeflMax_mm = bracketDeflMax * 1000;
out.gamaMax_mm = gamaMax * 1000;
out.loadShareRelRes = relRes;
out.workRelRes = workRes;
end

function s = fr_tf(b)
if b
    s = 'PASS';
else
    s = 'FAIL';
end
end
