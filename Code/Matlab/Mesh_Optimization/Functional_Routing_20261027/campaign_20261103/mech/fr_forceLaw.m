function [F, dom] = fr_forceLaw(C, L0, diaStr, kContraction, P)
%FR_FORCELAB Verified BPA force law wrapper.
%
% Zero-contraction capacity from the repository law of record
% Code/Matlab/Functions/maxBPAforce.m (mechanics_20261103/mech1 log:
% reproduces the handoff's ~441 / ~1314 N at 620 kPa, 450 mm to 0.04%).
% Strain-dependent force along the sweep comes from festo4 through
% BiPam_pulley2_20261022 (Fmax * festo4(D, relstrain, P)); this wrapper is
% for sizing/prints and for the contraction-domain flag.
%
% dom.extrapolated: logical per kContraction - TRUE where the relative
% strain leaves the campaign contraction band [0, KMAX] (handoff item 2
% forbids a universal 25% limit; extrapolation must be FLAGGED, not
% silently applied).
if nargin < 4
    kContraction = 0;
end
if nargin < 5
    P = C.pressure;
end
F = maxBPAforce(L0, diaStr, P);
KMAX = maxContractionFraction(C, diaStr);
relstrain = kContraction / L0 / KMAX;
dom.extrapolated = (relstrain < 0) | (relstrain > 1);
end

function k = maxContractionFraction(C, diaStr)
switch diaStr
    case '10'
        k = C.maxContractionPct.d10 / 100;
    case '20'
        k = C.maxContractionPct.d20 / 100;
    otherwise
        error('fr_forceLaw:Diameter', 'unsupported diameter %s', diaStr);
end
end
