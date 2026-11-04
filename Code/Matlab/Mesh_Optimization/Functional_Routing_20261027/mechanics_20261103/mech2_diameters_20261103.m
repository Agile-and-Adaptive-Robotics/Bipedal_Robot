function mech2_diameters_20261103()
%MECH2_DIAMETERS_20261103 MechanCheck row 2: the 20 mm sleeve diameter
% endpoints (20 mm at resting length -> 38 mm at fully contracted length)
% and the diameter-length relation the repo actually carries, plus the
% 10 mm sleeve's OWN calibration (bpaR diameter-specific braid-angle
% correction) evaluated at the corrected 16.5% contraction limit.
%
% Evaluated:
%   Code/Matlab/Functions/bpaR.m         -- braid-geometry law
%       Din = D0*sqrt(1-(1-k)^2*cos(theta0)^2)/sin(theta0); Dout = Din + 2*H0
%       theta0 = 28.6 deg + per-diameter Martens/Boblan correction
%       (10 mm: -0.1 deg at bpaR.m:48; 20 mm: -0.073 deg at bpaR.m:51)
%   Code/Matlab/Functions/BPAdiameters.m -- same law, same corrections
%   Constant fully-inflated envelopes:
%       buildKneeExtContext20mm.m:195-197  geo.bpaRadius = 0.0385/2 (38.5 mm)
%       buildGeoExclusion.m:69             circleRadius = 0.0381/2 (38.1 mm)

repo = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(repo, 'Code', 'Matlab')));

fprintf('=== MechanCheck 2: sleeve diameter endpoints and laws ===\n\n')

fprintf('--- bpaR.m law endpoints, 20 mm (own braid correction -0.073 deg) ---\n')
for km = [0.24, 0.255]
    Drest  = 2 * bpaR(0,    20, km) * 1000;
    Dfull  = 2 * bpaR(km,   20, km) * 1000;
    fprintf(['KMAX = %.3f: Dout at k=0 = %.3f mm (rest), at k=KMAX = %.3f mm ' ...
        '(fully contracted)\n'], km, Drest, Dfull)
end
% What contraction does the braid law need to reach the 38 mm endpoint?
Dtarget = 38;
f = @(k) 2 * bpaR(k, 20, 1.0) * 1000 - Dtarget;
k38 = fzero(f, [0.05, 0.60]);
fprintf('braid law reaches 38 mm Dout at k = %.4f (%.2f%% contraction)\n', k38, 100 * k38)
fprintf(['-> handoff endpoint pair 20 mm (rest) -> 38 mm (fully contracted) ' ...
    'is NOT reproduced by the braid law as one relation: it gives ' ...
    '23.6 mm at rest (nominal 20 mm is the INNER tube diameter; Dout adds ' ...
    '2*1.8 mm membrane) and ~%.1f mm at 24%%.\n'], 2 * bpaR(0.24, 20, 0.24) * 1000)
fprintf('-> repo instead carries a measured fully-inflated OUTER envelope as\n')
fprintf('   a CONSTANT: 38.5 mm (buildKneeExtContext20mm.m:195-197) and\n')
fprintf('   38.1 mm circle (buildGeoExclusion.m:69).\n\n')

fprintf('--- 10 mm sleeve OWN calibration (bpaR 10 mm branch, -0.1 deg) ---\n')
for km = [0.165, 0.255]
    Drest10 = 2 * bpaR(0,    10, km) * 1000;
    Dfull10 = 2 * bpaR(km,   10, km) * 1000;
    fprintf('KMAX = %.3f: 10 mm Dout at k=0 = %.3f mm, at k=KMAX = %.3f mm\n', ...
        km, Drest10, Dfull10)
end
D10_165_own = 2 * bpaR(0.165, 10, 0.165) * 1000;
D10_165_20law = 2 * bpaR(0.165, 20, 0.165) * (10 / 20) * 1000;
fprintf('10 mm inflated diameter at the corrected 16.5%% contraction:\n')
fprintf('  OWN 10 mm law  : %.3f mm Dout\n', D10_165_own)
fprintf('  20 mm law scaled (NOT ALLOWED, shown for the delta only): %.3f mm Dout\n', ...
    D10_165_20law)
fprintf('  delta (would-be error if the 20 mm law were reused): %+.3f mm\n', ...
    D10_165_20law - D10_165_own)
fprintf(['-> the 10 mm sleeve HAS its own calibration in the repo: the ' ...
    'per-diameter Martens/Boblan braid-angle correction at bpaR.m:46-59 ' ...
    '(10 mm: dtheta0 = -0.1 deg) and BPAdiameters.m:13. The 20 mm law is ' ...
    'never required for 10 mm candidates.\n\n'])

fprintf('--- consumers in the campaign-reusable path ---\n')
fprintf(['predictKneeFlexor20mm.m:46,131 call bpaR(strain, ctx.Dia, ctx.KMAX) ' ...
    'with ctx.KMAX = 0.255 (buildKneeFlexorContext20mm.m:141) -- the 20 mm\n'])
fprintf(['knee-flexor context. For a 10 mm sleeve the same call site evaluates ' ...
    'the 10 mm branch automatically (bpaR switches on round(d_mm)).\n'])
fprintf('=== MechanCheck 2 done ===\n')
end
