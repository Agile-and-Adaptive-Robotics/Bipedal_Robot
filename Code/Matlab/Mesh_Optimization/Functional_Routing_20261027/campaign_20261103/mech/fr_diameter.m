function [Dmm, info] = fr_diameter(C, diaStr, kContraction)
%FR_DIAMETER Sleeve diameter-contraction relation of THIS campaign.
%
% master item 2 requires (a) diameter = 20 mm at resting length, (b) 38 mm
% at the calibrated fully contracted length, (c) a physically consistent
% relation between the endpoints. mechanics_20261103/mech2 established
% that the repo's braid law (bpaR) does NOT reproduce the handoff endpoint
% pair as one relation for 20 mm (it gives Dout 23.6 mm at k=0 - nominal
% 20 mm is the INNER tube), so this campaign carries the handoff's
% calibrated endpoints for 20 mm with LINEAR interpolation (documented in
% fr_setup.C.diameterLaw), and the 10 mm sleeve's OWN braid law endpoints
% (Dout 13.6 mm at rest, 17.837 mm at the 16.5% corrected limit; the 20 mm
% law is never applied to 10 mm).
%
% Dmm: outer envelope diameter in mm at contraction fraction kContraction.
% info: endpoints + method + (for 10 mm) the own-law check values.
if strcmp(diaStr, '20')
    d0 = C.diameterEndpoints.d20_rest_mm;
    d1 = C.diameterEndpoints.d20_contracted_mm;
    k1 = C.maxContractionPct.d20 / 100;
elseif strcmp(diaStr, '10')
    d0 = C.diameterEndpoints.d10_rest_mm;
    d1 = C.diameterEndpoints.d10_contracted_mm;
    k1 = C.maxContractionPct.d10 / 100;
else
    error('fr_diameter:Diameter', 'unsupported diameter %s', diaStr);
end
kk = min(max(kContraction, 0), k1);
Dmm = d0 + (d1 - d0) * (kk / k1);
info = struct('rest_mm', d0, 'contracted_mm', d1, ...
    'contractionFraction_full', k1, 'law', C.diameterLaw);
if strcmp(diaStr, '10')
    % own-law cross-check values from mechanics_20261103/mech2 log
    info.ownLawRestDout_mm = 13.600;
    info.ownLawContractedDout_mm = 17.837;
    info.twentyMmLawApplied = false;
end
end
