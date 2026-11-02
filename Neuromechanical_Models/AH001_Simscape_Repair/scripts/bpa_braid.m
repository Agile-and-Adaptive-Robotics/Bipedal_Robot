function out = bpa_braid(L, P, cfg)
%BPA_BRAID  Diameter-expanding BPA (braided pneumatic actuator) static model.
%   Ben's directives (2026-11-01/02):
%     - "10/20/40 mm BPA" = UNINFLATED initial diameter D0.
%     - At 620 kPa: 10 mm contracts ~17%, OD -> ~22.5 mm;
%                   20 mm contracts ~25%, OD -> ~38 mm;
%                   40 mm contracts ~32%, OD -> ~75 mm.
%     - Endcaps stick out ~22.5 mm from EACH end (packaging: mount-to-mount
%       length = L + 45 mm total).
%   Braid law (inextensible braid):  b^2 = L^2 + (pi*N*D)^2
%       D(L) = sqrt(b^2 - L^2)/(pi*N)          % diameter GROWS as L shrinks
%   Static force (virtual work):     F = P/(4*pi*N^2) * (3L^2 - b^2), tension
%   only, active window Lmin <= L <= L0.
%
%   cfg fields:
%     L0   rest length [m]          (required in practice)
%     D0   un-inflated diameter [m] (default 0.020)
%     Kmax max contraction fraction  (default: per-size calibration below)
%     N    braid turns              (default: per-size calibration below)
%     Endcap  axial endcap protrusion per end [m] (default 0.0225);
%             out.mount_length = L + 2*Endcap (packaging)
%
%   Per-size calibration (N chosen so the model reproduces Ben's 620 kPa
%   endpoints exactly, given each actuator's own L0):
%     10 mm: Kmax 0.17, OD 22.5 mm ; 20 mm: Kmax 0.25, OD 38 mm ;
%     40 mm: Kmax 0.32, OD 75 mm.  N is solved per call from (D0, Kmax, OD@Kmax):
%       1/(pi N)^2 = (Dk^2 - D0^2) / (L0^2 - Lk^2),  Lk = L0*(1-Kmax)
%   i.e. N = sqrt((L0^2 - Lk^2)/(pi^2*(Dk^2 - D0^2))) -- replace with the
%   measured braid turn count when available.

if nargin < 3, cfg = struct(); end
D0 = getf(cfg, 'D0', 0.020);
% per-size calibration table
switch D0
    case 0.010, Kmax_d = 0.17; Dk = 0.0225;
    case 0.020, Kmax_d = 0.25; Dk = 0.038;
    case 0.040, Kmax_d = 0.32; Dk = 0.075;
    otherwise,  Kmax_d = 0.25; Dk = 1.9*D0;  % uncalibrated size: extrapolate
end
L0    = getf(cfg, 'L0', 0.2667809648);
Kmax  = getf(cfg, 'Kmax', Kmax_d);
Endcap= getf(cfg, 'Endcap', 0.0225);

Lmin = L0 * (1 - Kmax);
if isfield(cfg, 'N') && ~isempty(cfg.N)
    piN = pi * cfg.N;
else
    % solve pi*N from the calibrated endpoint (D at Lmin equals Dk)
    piN = sqrt((L0^2 - Lmin^2) / max(Dk^2 - D0^2, eps));
end
b2 = L0^2 + (piN*D0)^2;

Lc = min(max(L, Lmin), L0);
D  = sqrt(max(b2 - Lc^2, 0)) / piN;
Fraw = P / (4*piN^2) * (3*Lc^2 - b2);
F = max(Fraw, 0);

out.F = F;
out.D = D;
out.contraction = (L0 - Lc) / L0;
out.active = F > 0;
out.L = Lc;
out.mount_length = Lc + 2*Endcap;   % packaging: includes both endcaps
out.piN = piN;
out.b2 = b2;
end

function v = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
