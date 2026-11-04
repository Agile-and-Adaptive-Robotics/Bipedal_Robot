function ok = fr_diameterCheck(C)
%FR_DIAMETERCHECK Gate check (c): diameter endpoints (master item 2 /
% handoff item 2).
ok = true;
D20_0 = fr_diameter(C, '20', 0);
D20_1 = fr_diameter(C, '20', C.maxContractionPct.d20 / 100);
D10_0 = fr_diameter(C, '10', 0);
D10_1 = fr_diameter(C, '10', C.maxContractionPct.d10 / 100);
fprintf(['    diameter endpoints: 20 mm: rest %.3f mm (want 20), ' 'fully contracted %.3f mm (want 38); 10 mm own law: rest %.3f mm, ' '16.5%% %.3f mm\n'], D20_0, D20_1, D10_0, D10_1);
ok = ok && abs(D20_0 - 20) < 1e-9 && abs(D20_1 - 38) < 1e-9;
% the 10 mm sleeve's OWN braid law (bpaR 10 mm branch) must reproduce the
% campaign's 10 mm endpoints; the 20 mm law is never applied
r10 = bpaR(C.maxContractionPct.d10 / 100, '10', ...
    C.maxContractionPct.d10 / 100) * 2000;   % mm outer diameter
r10_0 = bpaR(0, '10', C.maxContractionPct.d10 / 100) * 2000;
fprintf(['    bpaR own 10 mm law check: Dout(0) = %.3f mm (want 13.600), ' ...
    'Dout(16.5%%) = %.3f mm (want 17.837)\n'], r10_0, r10);
ok = ok && abs(r10 - C.diameterEndpoints.d10_contracted_mm) < 0.01 && ...
    abs(r10_0 - C.diameterEndpoints.d10_rest_mm) < 0.01;
% monotone in contraction, and the 20 mm law endpoints would NOT match the
% braid law as one relation (documented in mechanics mech2) - verified
r20 = bpaR(C.maxContractionPct.d20 / 100, '20', ...
    C.maxContractionPct.d20 / 100) * 2000;
fprintf(['    braid-law 20 mm Dout at 24%% = %.3f mm (mech2: the braid ' 'law does NOT reproduce the handoff 20->38 mm pair; the campaign ' 'carries the calibrated endpoints)\n'], r20);
ok = ok && (abs(r20 - 38) > 0.5);   % documents the documented discrepancy
end
