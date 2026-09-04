function R = bpaR(k,d,KMAX)
%BPAR Outer BPA radius as a function of contraction.
%
% R = bpaR(k,d,KMAX)
%
% k     = BPA contraction (strain), e.g. 0.10 for 10% contraction
% d     = nominal BPA inner diameter. Accepted inputs:
%         10, 0.010, '10', ...
%         20, 0.020, '20', ...
%         40, 0.040, '40', ...
% KMAX  = maximum allowable BPA contraction, e.g. 0.255
%
% The braid geometry follows the geometric relationship described by
% Chou and Hannaford (1996). The initial braid angle and diameter-specific
% correction follow Martens and Boblan (2017).
%
% The braid fiber length and number of turns are eliminated algebraically,
% so the instantaneous inner diameter can be calculated directly from the
% fractional BPA contraction k.
%
% The membrane thickness H0 is added to the calculated inner diameter to
% estimate the physical outer diameter:
%
%       Dout = Din + 2*H0
%
% and the returned value is the outer radius:
%
%       R = Dout/2

%% Diameter input
if ischar(d) || isstring(d)
    d = str2double(d);
end

if d < 1
    D0 = d;               % Input already in m
    d_mm = 1000*d;
else
    d_mm = d;             % Input in mm
    D0 = d/1000;
end

%% Initial braid angle
theta0 = 28.6;            % deg, Chou and Hannaford / Martens and Boblan

switch round(d_mm)
    case 10
        dtheta0 = -0.1;   % deg, Martens and Boblan correction

    case 20
        dtheta0 = -0.073; % deg, Martens and Boblan correction

    case 40
        dtheta0 =  0;       % Not investigated

    otherwise
        error('bpaR:UnsupportedDiameter', ...
            'BPA diameter must be 10, 20, or 40 mm.')
end

theta0 = deg2rad(theta0 + dtheta0);

%% Limit contraction to physical BPA range
k = max(0,min(k,KMAX));

%% Membrane thickness
H0 = 0.0018;              % m

%% Instantaneous BPA diameter
% From:
%
%   Lfiber = L0/cos(theta0)
%   n      = L0*tan(theta0)/(pi*D0)
%   L      = L0*(1-k)
%
% Substitution eliminates L0, Lfiber, and n.

Din = D0 .* ...
    sqrt(1 - (1-k).^2 .* cos(theta0).^2) ./ ...
    sin(theta0);

Dout = Din + 2*H0;

%% Outer BPA radius
R = Dout/2;

end