%% Optimize predicted torque for extensors.
function [f_all, bpa_all] = minimizeFlxPinX3_2brkt(Xi0,Xi1,Xi2,Xi3,idx_val,useB2t)
% minimizeFlxPinX3_2brkt: TWO-BRACKET variant of minimizeFlxPin (Ben 2026-09-15).
% The original single-bracket evaluator is preserved untouched in
% minimizeFlxPinX3.m -- call that one to revert.
% minimizeExt: calculates predicted torque and fit metrics for a given BPA index
%
% Inputs:
%   Xi0 - extra length correction (m)
%   Xi1 - bracket "axial" stiffness (N/m)
%   Xi2 - bracket "bending" stiffness (N/m)
%   whichIdx - index of the BPA to process (scalar)
%   useB2t - (optional, default true) second tibia bracket on/off (Ben 2026-09-14)
%
% Two tibia brackets (Ben 2026-09-14): bracket 1 at the insertion (Pbri) as before;
% bracket 2 at Pbri2 = [30.5, -103.41, 0]/1000 from the knee ICR (z=0 assumed),
% K2t = [X1, X2, X1] sharing Xi1/Xi2. Both see the class's tibia-frame force vector
% (Fk = klass.unitD .* FF) and act as series compliances: c_eff = c_b + c_b2 + cSpr.
% Bracket-2 deflections are stored on bpa (eB2 bracket frame, e_tib2 tibia frame).
% Screw head (Ben 2026-09-15): bracket-2 tibia-frame X deflection is CLAMPED at
% -3.5 mm -- free until the tip bears on the 6 mm screw head, then a contact
% normal force Nc pins X there, only the unconstrained directions keep
% deflecting, and the chain equilibrium is re-solved to consistency (the
% contraction changes once clamped: part of the bracket force transfers into
% normal contact instead of deformation). bpa.screw2_N carries Nc per frame;
% screw2_hit flags any contact.
% Tangency check (reject -> GoF [Inf Inf Inf], active only for finite Xi1/Xi2):
%   - 30 mm-dia circle: on non-colinear frames the force line from the DEFLECTED
%     bracket-2 position must stay tangent-or-outside (d >= 15 mm about knee center);
%     relax line = Location(end-1)-Location(end-2) for >=3-row paths, else the fixed
%     bracket2->insertion line; colinear within 2 deg.
% Xi3 (Ben 2026-09-15): UNITLESS wrap-loss factor in [0,1] -- the extensor meaning,
% NOT the old series stiffness. delta_L = Xi3 * R * theta_wrap * comp^2 per frame,
% R = 15 mm (the 30 mm-dia knee circle), theta_wrap = pi minus the angle between
% the class force direction (toward the origin) and the bracket2->insertion
% direction (colinear = straight-through = no wrap), comp = max(0, 1-relstrain).
% The series element is tendon-only: kspr = Spr(klass) (Inf when ten = 0).
% Escape/baseline: Xi1 = Xi2 = Inf, Xi3 = 0.
% Data: persistent kf cache; default adds the 47 cm (kf(3)) +5.3 deg encoder
% correction on EXPERIMENTAL ANGLES ONLY (matches minimizeFlxPin.m). Env
% FLXPX3_NOSHIFT=1 restores the legacy unshifted X3 behavior.
%
% Outputs:
%   fitvec - [RMSE, FVU, MaxResidual] for the selected BPA
%   bpa    - updated BPA struct with prediction fields filled in

%% load
%kf = knee flexor, kf(1) = specific resting length, 
%ke = knee extensor, same as above
%example: kf(1).Mz z-axis torque for pinned knee, flexor, 46cm length
%'exp' suffix means experimentally measured
%'_h' suffix means hybrid calculated
%'_p' suffix means prime, as in the new prediction values

% Data: persistent cache, default WITH the 47 cm (kf(3)) encoder correction.
% FLXPX3_NOSHIFT=1 -> legacy unshifted X3 behavior.
useShift = ~strcmpi(getenv('FLXPX3_NOSHIFT'), '1');
kf = loadKfCached(useShift); %This loads the following, which was ran and saved:

    % load KneeFlxPin_10mm_48cm.mat phiD
    % load Plot_KneeFlxPin10mm_48cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand G Bifemsh_Pam
    % A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    % kf(1) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
    %               'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
    %               'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
    %               'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
    %               'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
    %               'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
    %               'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
    %               'Lmt_p', [], 'mA_p', [], 'M_p', [], 'F_p', [], 'strain_p', [], 'L_p', [], 'gama', [], 'strain_f', [],'Lm_f', []);
    % clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
    % 
    % 
    % % 46cm length
    % load KneeFlxPin_10mm_46cm.mat phiD
    % load Plot_KneeFlxPin10mm_46cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand G Bifemsh_Pam
    % A = sortrows([Angle, Torque, InflatedLength', ICRtoMuscle', TorqueHand]);
    % kf(2) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
    %               'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
    %               'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
    %               'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
    %               'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
    %               'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
    %               'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
    %               'Lmt_p', [], 'mA_p', [], 'M_p', [], 'F_p', [], 'strain_p', [], 'L_p', [], 'gama', [], 'strain_f', [],'Lm_f', []);
    % clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
    % 
    % % 47cm length
    % load KneeFlxPin_10mm_47cm.mat phiD
    % load Plot_KneeFlxPin10mm_47cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand G Bifemsh_Pam
    % A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    % kf(3) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
    %               'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
    %               'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
    %               'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
    %               'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
    %               'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
    %               'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
    %               'Lmt_p', [], 'mA_p', [], 'M_p', [], 'F_p', [], 'strain_p', [], 'L_p', [], 'gama', [], 'strain_f', [],'Lm_f', []);
    % clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
    % 
    % % 40cm length
    % load KneeFlxPin_10mm_40cm.mat phiD
    % load Plot_KneeFlxPin10mm_40cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand G Bifemsh_Pam
    % A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    % kf(4) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
    %               'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
    %               'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
    %               'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
    %               'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
    %               'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
    %               'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
    %               'Lmt_p', [], 'mA_p', [], 'M_p', [], 'F_p', [], 'strain_p', [], 'L_p', [], 'gama', [], 'strain_f', [],'Lm_f', []);
    % clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
    % 
    % % 42cm length
    % load KneeFlxPin_10mm_42cm.mat phiD
    % load Plot_KneeFlxPin10mm_42cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand G Bifemsh_Pam
    % A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    % kf(5) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
    %               'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
    %               'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
    %               'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
    %               'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
    %               'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
    %               'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
    %               'Lmt_p', [], 'mA_p', [], 'M_p', [], 'F_p', [], 'strain_p', [], 'L_p', [], 'gama', [], 'strain_f', [],'Lm_f', []);
    % clear Bifemsh_Pam phiD G Angle Torque InflatedLength ICRtoMuscle TorqueHand A

%% Initialize output
nBPA = numel(kf);
% Default to all BPAs if none specified
if nargin < 5 || isempty(idx_val)
        idx_val = 1:nBPA;
end
% Second tibia bracket on by default (Ben 2026-09-14); pass false to disable
if nargin < 6 || isempty(useB2t)
        useB2t = true;
end


bpa_all = kf;  % initialize
f_all = NaN(nBPA, 3);


%% Evaluate each BPA
for i = idx_val
%     fprintf('Evaluating BPA #%d with [%.4f, %.2e, %.2e]\n', i, Xi0, Xi1, Xi2, Xi3);
    klass_i = kf(i);
    [bpa_all(i), f_all(i,:)] = evaluateBPA(klass_i, Xi0, Xi1, Xi2, Xi3, useB2t);
    if any(isnan(bpa_all(i).strain_p))
        warning('NaNs in strain_p for BPA #%d', i);
    end
end

end


function [bpa_i, fitvec] = evaluateBPA(klass, Xi0, Xi1, Xi2, Xi3, useB2t)
%% Calculate locations and properties
bpa_i = klass;  %
kspr = Spr(bpa_i); %Calculate spring rate (tendon only; Infinite if no tendon) -- Xi3 is the unitless wrap-loss factor now
strain_Xi0 = Contraction(bpa_i, [],Xi0); %Calculate contraction with constant length offset
[L_p, gemma, eB2, e_tib2, chk, delta_L] = Lok(bpa_i, Xi1, Xi2, kspr, strain_Xi0, Xi0, useB2t, Xi3);   %Bracket deformation changing geometry + Xi3 wrap loss (2 brackets, Ben 2026-09-14/15)
unitD_p = UD(bpa_i, L_p);   %New force direction
sL_p = seg(bpa_i, L_p);   %New segment lengths uses deformation but does not subtract length offset
Lmt_p = LMT(sL_p, Xi0+gemma);     %New musclulotendon length. Uses deformed geometry and constant length offset.
strain_p = Contraction(bpa_i, Lmt_p - delta_L, []);  %*new contraction: deformed geometry, length offset, Xi3 wrap loss
F_p = Force(bpa_i, unitD_p, strain_p);  %new force vector
mA_p = Mom(bpa_i, L_p, unitD_p);   %new moment arm
% M_p = Tor(mA_p, F_p, bpa_i.Fm, strain_p);  %new torque
M_p = Tor(mA_p, F_p, strain_p);  %new torque

%% Package into output struct
% bpa = bpa_i;
bpa_i.Lmt_p = Lmt_p;
bpa_i.mA_p = mA_p;
bpa_i.M_p = M_p;
bpa_i.F_p = F_p;
bpa_i.strain_p = strain_p;
bpa_i.L_p = L_p;
bpa_i.gama = gemma;
bpa_i.delta_L = delta_L;     %Xi3 wrap loss per frame (m)
bpa_i.theta_wrap = chk.theta_wrap;   %wrap angle on the 30 mm circle per frame (rad)
bpa_i.eB2 = eB2;             %bracket-2 deflections [axial bendY bendZ], bracket frame, Nx3
bpa_i.e_tib2 = e_tib2;       %bracket-2 deflections, tibia frame, Nx3
bpa_i.screw2_hit = chk.screw2_hit;     %bracket-2 bore on the screw head somewhere (clamped)
bpa_i.screw2_N = chk.screw2_N;         %screw-head contact normal force per frame (N)
bpa_i.tangency_ok = chk.tangency_ok;   %30 mm-circle tangency held on all non-colinear frames
bpa_i.d_circle2 = chk.d_circle2;       %perp. distance knee center -> force line, per frame
bpa_i.colinear2 = chk.colinear2;       %relax-line colinear flag, per frame

% GoF calculation
fitvec = SSE(bpa_i, M_p);
% Physical-limit rejection for bracket 2 (Ben 2026-09-14): screw head or tangency
if chk.reject
    fitvec = [Inf, Inf, Inf];
end

%% Nested functions, modified from MonoPamExplicit

%% -------------- Contraction of the PAM --------------------------
function contraction = Contraction(klass, L_mt, X0)
rest   = klass.rest;
tendon = klass.ten;
fitting = klass.fitn;

if isempty(X0)
    X0 = 0;
end

if isempty(L_mt)
    L_mt = klass.Lmt;
end

contraction = (rest-(L_mt-tendon-2*fitting-X0))/rest;    %(minus Xi0 is also used in LMT function)
end

%% ------------- Location  ------------------------
function [LOC, gama, eB2, e_tib2, chk, delta_L] = Lok(klass,X1,X2,kSpr,strain,X0,useB2t,Xi3)
% Inputs:
%   bpa class info
%   X1, X2 stiffness
%   kSpr, tendon stiffness
%   Funit, force unit direction in the hip frame
%   strain – N×1 strain vector (e.g., from Xi0 offset effect)
%   X0, constant length offset
L = klass.Loc;      %Location (wrapping, attachment points)
C = klass.CP;       %Cross point (moves from one frame to another)
%             T = klass.Tk;       %Transformation matrix
kmax = klass.Kmax;  %max contracted length
KMAX = (klass.rest-kmax)/klass.rest; %turn it into a percentage

relstrain = strain/KMAX;            %relative strain
FF = festo4(klass.dBPA,relstrain,klass.P).*klass.Fm;        %Force magnitude
FF(FF<0) = 0;
%For muscle insertion
unitD = klass.unitD;            %unit direction of force vector, tibia frame
Fk = unitD.*FF;                  %Force vector, tibia frame
pB = L(C,:,(klass.Ak==0));                  %Distance from knee frame to muscle insertion
% Pbri = [-48.11, -107.81, 13.8]/1000;     %vector from knee ICR to flexor insertion bracket (where it starts to cantilever)
Pbri = [-27.5, -107.81, -0.54]/1000;     %vector from knee ICR to flexor insertion bracket (where it starts to cantilever, but at tibial contact, no z offset)
% Pbri = [-27.5, -125.91, -0.54]/1000;     %vector from knee ICR to upper bolt
pkbrB = pB-Pbri;                  %vector from bracket to point B, in the knee frame
thetabrB = atan2(pkbrB(2),pkbrB(1));   %angle between pbrB and x axis
RkbrZ = [cos(thetabrB) -sin(thetabrB) 0; ...     %Rotation matrix
       sin(thetabrB) cos(thetabrB) 0; ...
       0    0   1];
pbrkB = RkbrZ'*pkbrB';       %Vector in the bracket frame
%             Now calculate angle from x-axis to this vector
thetaY = atan2(pbrkB(3), pbrkB(1));  % z vs x (in bracket frame)
% % Rotation matrix about y-axis (local frame adjustment)
Ry = [cos(thetaY) 0  sin(thetaY);
      0           1  0;
     -sin(thetaY) 0  cos(thetaY)];
Rkbr = RkbrZ*Ry';            %Rotate about y-axis in body frame
Tkbr = RpToTrans(RkbrZ, Pbri');    %Transformation matrix, flexor bracket frame in knee frame

% ---- Second tibia bracket (Ben 2026-09-14): screw-head bracket at Pbri2 ----
% Same force vector as bracket 1 (the class's tibia-frame Fk) seen from its own
% pitch-only frame. K2t = [X1, X2, X1] is applied inside fortz.
if useB2t
    Pbri2 = [30.5, -103.41, 0]/1000;   %vector from knee ICR to 2nd tibia bracket (z=0 assumed)
    pkbrB2 = pB-Pbri2;                 %vector from bracket 2 to point B, in the knee frame
    thetabrB2 = atan2(pkbrB2(2),pkbrB2(1));   %angle between pbrB2 and x axis
    RkbrZ2 = [cos(thetabrB2) -sin(thetabrB2) 0; ...     %Rotation matrix
           sin(thetabrB2) cos(thetabrB2) 0; ...
           0    0   1];
    Tkbr2 = RpToTrans(RkbrZ2, Pbri2'); %Transformation matrix, 2nd bracket frame in knee frame
else
    RkbrZ2 = eye(3); Tkbr2 = eye(4); Pbri2 = zeros(1,3);
end

LOC = L;            %new location matrix
N = size(L,3);
%             M = size(L,1);
Fbrk = zeros(N,3);       %Force vector represented in the tibial bracket frame
Fbrk2 = zeros(N,3);      %Force vector represented in the 2nd tibia bracket frame

parfor ii = 1:N                          %Repeat for each orientation
        Fbrk(ii,:) = RowVecTrans(Tkbr\eye(4),Fk(ii,:)); %Force vector in the tibia frame represented in the lower bracket frame
    if useB2t
        Fbrk2(ii,:) = RowVecTrans(Tkbr2\eye(4),Fk(ii,:)); %same class force vector, 2nd bracket frame
    end
end

if isinf(X1) && isinf(X2)  && isinf(kSpr)
    [epsilon, delta, beta, gama, e_ax2t, e_by2t, e_bz2t, scr2hit, scr2N] = deal(zeros(N,1));
else
    [epsilon, delta, beta, gama, e_ax2t, e_by2t, e_bz2t, scr2hit, scr2N] = fortz(klass,Fbrk,Fbrk2,RkbrZ2,X1,X2,kSpr,X0,useB2t);  %strain from force divided by tensile stiffness
end

% ---- Xi3 wrap loss (unitless, Ben 2026-09-15): cable arc on the 30 mm-dia knee circle
[delta_L, theta_wrap] = wrapLoss(klass, L, C, Pbri2, strain, Xi3, useB2t, N);

eB = [epsilon, delta, beta];
pbrBnew = [norm(pkbrB(1:2)), 0, pkbrB(3)] + eB; %new point B, in the bracket's frame
 % pbrBnew = [norm(pkbrB), 0, 0] + eB; %new point B, in the bracket's frame

pBnew = zeros(N,3);
for ii = 1:N                          %Repeat for each orientation
    pBnew(ii,:) = RowVecTrans(Tkbr, pbrBnew(ii,:));     %New point B, in the tibia frame
    LOC(2,:,ii) = pBnew(ii,:);
end

% ---- bracket-2 deflections in the tibia frame + physical checks ----
eB2 = [e_ax2t, e_by2t, e_bz2t];
e_tib2 = (RkbrZ2 * eB2.').';           %bracket-2 deflection vector in the tibia frame, Nx3
chk = checkBracket2(klass, L, pB, Pbri2, e_tib2, X1, X2, useB2t);
chk.screw2_hit = any(scr2hit);   %true when bracket 2 bore on the screw head somewhere
chk.screw2_N = scr2N;            %screw-head contact normal force per frame (N), 0 when free
chk.theta_wrap = theta_wrap;     %Xi3 wrap angle per frame (rad)
end

%% ------------- Xi3 wrap loss (unitless, Ben 2026-09-15) ----------------------
% delta_L = Xi3 * R * theta_wrap * comp^2 per frame. R = 15 mm (the 30 mm-dia
% circle about the knee center). theta_wrap = pi minus the angle between the
% class force direction (toward the origin) and the bracket2->insertion
% direction: a cable running straight through the bracket wraps nothing.
% comp = max(0, 1 - relstrain) -- the additive complement to relative strain.
function [delta_L, theta_wrap] = wrapLoss(klass, L, C, Pbri2, strain, Xi3, useB2t, N)
theta_wrap = zeros(N,1);
delta_L = zeros(N,1);
if ~useB2t || ~isfinite(Xi3) || Xi3 <= 0
    return;
end
KMAX = (klass.rest - klass.Kmax)/klass.rest;
uin = klass.unitD;                       %force unit direction, tibia frame (class-computed)
for ii = 1:N
    v = L(C,:,ii) - Pbri2;               %bracket2 -> insertion, tibia frame
    n = norm(v);
    if n < 1e-9, continue; end
    d = max(-1, min(1, dot(uin(ii,:), v/n)));
    theta_wrap(ii) = pi - acos(d);       %0 when straight through, pi when folded back
end
comp = max(0, 1 - strain(:)/KMAX);
delta_L = Xi3 * 0.015 * theta_wrap .* comp.^2;
end

%% Force and length reduction due to deformation
function [e_axial, e_bendY, e_bendZ, e_cable, e_ax2t, e_by2t, e_bz2t, screw2_hit, screw2_N] = fortz(klass,Fbr,Fbr2,RkbrZ2,X1,X2,kSpr,X0,useB2t)
% e_axial, bracket axial elongation
% e_bendY, bracket bending displacement y - direction
% e_bendZ, bracket bending displacement z - direction
% e_cable, tendon cable stretch
% e_ax2t/e_by2t/e_bz2t, 2nd tibia bracket deflection (its own frame), Ben 2026-09-14
% screw2_hit/screw2_N, screw-head contact flag + normal force per frame (Ben 2026-09-15)
% total length change
    N = size(Fbr,1);
    % Initialize outputs
    [e_axial, e_bendY, e_bendZ, e_cable, e_ax2t, e_by2t, e_bz2t, screw2_hit, screw2_N] = deal(zeros(N, 1));

    if isinf(X1) && isinf(X2) && isinf(kSpr)
        return
    end
    
    if isempty(X0)
        X0=0;
    end
    D = klass.dBPA;         %BPA diameter
    rest = klass.rest;      %resting length
    tendon = klass.ten;     %tendon length
    fitn = klass.fitn;    %fitting length
    mL = klass.Lmt - X0 -tendon -2*fitn;       %Musculotendon length
    mif = klass.Fm;         %maximum force
    kmax = klass.Kmax;      %maximum contracted length
    KMAX = (rest-kmax)/rest; %turn it into a percentage
    P = klass.P;            %pressure
    
    % Normalize force vectors safely
    norms = vecnorm(Fbr, 2, 2);
    valid = norms > 1e-4 & all(~isnan(Fbr), 2);
    u_hat_all = normalize(Fbr);
    
    % Vectorized k_b computation
    K = [X1, X2, X1];  %bracket stiffness array
    K_bracket = diag(K);       %bracket stiffness matrix
    C_bracket = diag([1/K(1), 1/K(2), 1/K(3)]); %bracket compliance
    K2t = [X1, X2, X1];        %2nd tibia bracket stiffness array (Ben 2026-09-14, shares Xi1/Xi2)
    K_bracket2 = diag(K2t);    %2nd bracket stiffness matrix
    u_hat = permute(u_hat_all, [3, 2, 1]);  % [1x3xN]
    C_rep = repmat(C_bracket, [1, 1, N]);   % [3x3xN]
    c_b = pagemtimes(pagemtimes(u_hat, C_rep), permute(u_hat, [2, 1, 3]));
    c_b = reshape(c_b, [N, 1]);
    % second tibia bracket: same class force vector in its own frame
    u_hat2_all = normalize(Fbr2);
    valid2 = vecnorm(Fbr2, 2, 2) > 1e-4 & all(~isnan(Fbr2), 2);
    if useB2t
        u_hat2 = permute(u_hat2_all, [3, 2, 1]);  % [1x3xN]
        C_rep2 = repmat(diag([1/K2t(1), 1/K2t(2), 1/K2t(3)]), [1, 1, N]);
        c_b2 = pagemtimes(pagemtimes(u_hat2, C_rep2), permute(u_hat2, [2, 1, 3]));
        c_b2 = reshape(c_b2, [N, 1]);
        c_b2(~valid2) = 0;
    else
        c_b2 = zeros(N, 1);
    end
    cSpr = 1/kSpr;
    c_eff = c_b + c_b2 + cSpr;        %effective compliance: both tibia brackets + series spring
    k_eff = 1 ./ c_eff;         % effective stiffness along u
    
    
    % Parallel root solve
   parfor i = 1:N
        if ~valid(i)
            continue;
        end
    
        % Per-instance constants
        keff = k_eff(i);
        unit_vec = u_hat_all(i, :);
        unit_vec2 = u_hat2_all(i, :);
        e_bkt2 = zeros(3,1);
        Lm = mL(i);
            
        contraction0    = ( rest - Lm ) / rest;
        relstrain0      = contraction0 / KMAX;  %relative strain
        if relstrain0 >= 1
            r = 0;
        else
            % Chain equilibrium with the screw-head clamp on bracket 2 (Ben 2026-09-15):
            % the tip deflects freely until tibia-frame X reaches -3.5 mm, then it
            % bears on the screw head -- a contact normal force Nc pins X there and
            % only the unconstrained directions keep deflecting. The projected
            % compliance c_b2 drops, so the equilibrium is re-solved to consistency
            % (fixed point on the secant compliance; piecewise-linear chain).
            c_b2_i = c_b2(i);
            e_bkt2 = zeros(3,1);
            r = 0;              %parfor temporary: define before the iteration loop
            for iter = 1:10
                keff = 1 / (c_b(i) + c_b2_i + cSpr);
                relfun = @(r) ...
                    festo4( D, ...
                        (rest - (Lm -  r)) / rest / KMAX, ...
                        P  ...
                    ) * mif - keff * r;

                try
                    r = fzero(relfun, [0, Lm-kmax]);
                catch
                    r = 0;
                end
                r = max(r,0); %guard against r being slightly negative.
                if r == 0 || isinf(X1) || isinf(X2) || ~useB2t || ~valid2(i)
                    break;   %nothing to clamp-iterate
                end
                F_it = festo4(D, (rest-(Lm-r))/rest/KMAX, P) * mif;  %chain force at r
                e_free2 = K_bracket2 \ (F_it * unit_vec2');          %free bracket-2 deflection
                if RkbrZ2(1,:) * e_free2 >= -3.5e-3
                    c_new = c_b2(i);            %regime A: no contact
                    e_bkt2 = e_free2;
                else
                    % contact: pin tibia-frame X at -3.5 mm via normal force Nc
                    t2 = RkbrZ2(1,:).';                     %tibia-x axis in bracket-2 frame
                    a = RkbrZ2(1,:) * e_free2;
                    b = RkbrZ2(1,:) * (K_bracket2 \ t2);
                    Nc = (-3.5e-3 - a) / b;
                    e_bkt2 = K_bracket2 \ (F_it * unit_vec2' + Nc * t2);
                    c_new = (unit_vec2 * e_bkt2) / F_it;    %secant absorbed compliance
                    screw2_hit(i) = true;
                    screw2_N(i) = Nc;
                end
                converged = abs(c_new - c_b2_i) <= 1e-9 * max(c_b2(i), 1e-12);
                c_b2_i = c_new;
                if converged
                    break;
                end
            end
        end
        
        if r == 0
            continue;
        elseif isinf(X1) && isinf(X2)
            % Rigid body: no bracket deformation
            e_axial(i) = 0;
            e_bendY(i) = 0;
            e_bendZ(i) = 0;
            e_cable(i) = r;  % All elongation goes to cable
            continue;
        else
            % Final force magnitude
            contraction = (rest - (Lm - r)) / rest;
            relstrain = contraction / KMAX;
            F_mag = festo4(D, relstrain, P) * mif;
    
            % Bracket displacement
            e_bkt = K_bracket \ (F_mag * unit_vec');

            e_axial(i) = e_bkt(1);
            e_bendY(i) = e_bkt(2);
            e_bendZ(i) = e_bkt(3);

            % 2nd tibia bracket displacement (Ben 2026-09-14): from the clamped
            % chain solve above -- free deflection, or pinned at the screw head
            % with contact normal force Nc when e_tib2 X would pass -3.5 mm
            if useB2t && valid2(i)
                e_ax2t(i) = e_bkt2(1);
                e_by2t(i) = e_bkt2(2);
                e_bz2t(i) = e_bkt2(3);
            end

            % Cable elongation
                if tendon > 0
                    r_bracket = unit_vec * e_bkt + unit_vec2 * e_bkt2;   %both brackets absorb
                    r_cable = r-r_bracket;
                    e_cable(i) = F_mag/kSpr;

                    if abs(r_cable - e_cable(i)) > 1e-6
                        warning('fortz:LengthBalanceMismatch', ...
                            'Frame %d: e_cable = %.9g, r_cable = %.9g, diff = %.9g', ...
                            i, e_cable(i), r_cable, r_cable - e_cable(i));
                    end
                end
        end
   end
end

%% ------------- Segment Lengths ------------------------
function SL = seg(klass, L)
C = klass.CP;
T = klass.Tk;
N = size(T, 3);
M = size(L, 1);
SL = zeros(N, M-1);

for ii = 1:N                    %Repeat for each orientation
    for i = 1:M-1               %Calculate all segments
        pointA = L(i,:,ii);
        pointB = L(i+1,:,ii);
        if i+1 == C
            pointB = RowVecTrans(T(:,:,ii), pointB);
        end
        SL(ii,i) = norm(pointA - pointB);
    end
end
end
               
        %% ------------- Muscle Length ------------------------
        %Function that calculates the musclutendon length
function Lmt = LMT(sL, X0)
% Compute muscle-tendon length from segment lengths and offset Xi0
% Xi0 can be empty [] to skip correction (i.e., when already applied)
% N = size(sL_p, 1);
Lmt = sum(sL, 2);  % Nx1, sum across segments

if ~isempty(X0)
    Lmt = Lmt - X0;
end
end          

        %% -------------- Force Unit Direction ----------------
        %Calculate the unit direction of the muscle force about the joint (tibia frame).
function unitD = UD(klass, L)
            T = klass.Tk;
            C = klass.CP;
            direction = zeros(size(T, 3), 3);
%             unitD = zeros(size(direction));
            
            for i = 1:size(T, 3)
                pointA = L(C-1, :, i);
                pointB = L(C, :, i);
                direction(i, :) = RowVecTrans(T(:, :, i)\eye(4), pointA) - pointB;
%                 unitD(i, :) = direction(i, :)/norm(direction(i, :));
            end
            unitD = normalize(direction);
end
        
        %% -------------- Moment Arm --------------------------
        %Calculate the moment arm about a joint
        %For every ViaPoint, calculate the moment arm of the muscle about
        %the joint it crosses over
function mA = Mom(klass, L, unitD)
            T = klass.Tk;
            C = klass.CP;
            mA = zeros(size(T, 3), 3);
            
            for i = 1:size(T, 3)
                pointB = L(C, :, i);
                mA(i, :) = pointB - unitD(i, :)*dot(unitD(i, :), pointB);
            end
end               


        %% -------------- Force --------------------------
        %Calculate the direction of the forced applied by the muscle
function F = Force(klass, unitD, strain)
        %Inputs:
        %Lmt == muscle-tendon length, scalar
        %rest == resting length of artificial muscle, "size" from Size function
        %dia == diameter of Festo tube, from Size function
        %pres == measured pressure
        %kmax == maximum contraction length
        %Outputs:
        %F == Force, N           
           rest = klass.rest;
           kmax = klass.Kmax;  
           KMAX = (rest-kmax)/rest; %turn it into a percentage 
            
           rel = strain./KMAX;                    %relative strain        
           
           Fn = festo4(klass.dBPA,rel,klass.P);

           scalarForce = Fn.*klass.Fm;
           scalarForce(scalarForce <= 0) = 0;
%            scalarForce(scalarForce > maxF) = NaN;            
            
            F = scalarForce.*unitD;

end
        
        %% ---------------------- Torque --------------
        %Calculate torque by multiplying the the force along the 
        %Useful information
        % i -> Index for Crossing Points/Joints
        % ii -> Index for every degree of motion
        % iii -> Index for axes of interest to observe Torque about
function Mz = Tor(mA, F, strain)
N = size(F, 1);
Mz = zeros(N, 3);

    for i = 1:N
        % if norm(F(i,:)) > maxF || strain(i,:) < -0.03
        if strain(i,:) < -0.03
            Mz(i,:) = NaN;
        else
            Mz(i,:) = cross(mA(i,:), F(i,:));
        end
    end
end   

%% tendon springrate (tendon only -- Xi3 is the unitless wrap-loss factor now)
function springrate = Spr(klass)
    if klass.ten > 0
        mult = 2;           %Multiplier for number of cables used.
        Aeff = 1.51*10^-6;  %Effective area for 19-strand cable
        E = 193*10^9;       %Young's Modulus
        L = klass.ten;      %tendon length
        springrate = mult*Aeff*E/L;
    else
        springrate = Inf;
    end
end

%% Subfunctions
function t = SSE(klass, M_p)
     [Ak_sorted, idx] = sort(klass.Ak);
     M_sorted = M_p(idx, 3);
     Mpredict2 = griddedInterpolant(Ak_sorted, M_sorted);
     M_opt = Mpredict2(klass.Aexp);
     [RMSE, fvu, maxResid] = Go_OfF(klass.Mexp,M_opt);
     t = [RMSE, fvu, maxResid];
end

function vhat = normalize(v)
    N = size(v,1);
    norms = vecnorm(v,2,2);
    valid = norms > 1e-4 & all(~isnan(v), 2);
    vhat = zeros(N, 3);
    vhat(valid, :) = v(valid, :) ./ norms(valid);
end


end

%% ------------- Cached kf data (Ben 2026-09-14) -------------------------------
% FlxPinBPASet + optional 47 cm (kf(3)) encoder correction: +5.3 deg on
% EXPERIMENTAL ANGLES ONLY (Aexp/A_h), matching minimizeFlxPin.m's build-time
% shift (Angle+5.3 then sortrows is order-preserving, so adding post-hoc to the
% sorted columns is identical). Cached per shift mode.
function kf = loadKfCached(useShift)
persistent kf_ns kf_sh
if useShift
    if isempty(kf_sh)
        S = load('FlxPinBPASet.mat', 'kf');
        kf_sh = S.kf;
        kf_sh(3).Aexp = kf_sh(3).Aexp + 5.3;   %47cm encoder correction, exp angles only
        kf_sh(3).A_h  = kf_sh(3).A_h  + 5.3;
        kf_sh = addKfPlaceholders(kf_sh);
    end
    kf = kf_sh;
else
    if isempty(kf_ns)
        S = load('FlxPinBPASet.mat', 'kf');
        kf_ns = addKfPlaceholders(S.kf);
    end
    kf = kf_ns;
end
end

% placeholder fields so bpa_all(i) = evaluateBPA(...) stays a same-shape struct
% assignment once evaluateBPA fills the bracket-2 diagnostics (Ben 2026-09-14)
function kf = addKfPlaceholders(kf)
for jD = 1:numel(kf)
    kf(jD).eB2 = [];         %bracket-2 deflections [axial bendY bendZ], bracket frame
    kf(jD).e_tib2 = [];      %bracket-2 deflections, tibia frame
    kf(jD).screw2_hit = [];
    kf(jD).screw2_N = [];
    kf(jD).delta_L = [];
    kf(jD).theta_wrap = [];
    kf(jD).tangency_ok = [];
    kf(jD).d_circle2 = [];
    kf(jD).colinear2 = [];
end
end

%% ------------- Bracket-2 physical checks (Ben 2026-09-14/15) -----------------
% Screw head: handled in fortz as a CLAMP (tibia-frame X pinned at -3.5 mm by a
% contact normal force once the tip bears on the 6 mm screw head) -- NOT a
% rejection. screw2_hit/screw2_N are merged into chk by Lok.
% 30 mm-dia circle: on frames where the class force direction is NOT colinear
% (within 2 deg) with the relax line, the force line from the DEFLECTED
% bracket-2 position must stay tangent-or-outside the circle about the knee
% center (tibia origin): perpendicular distance >= 15 mm. Violation rejects the
% candidate (GoF = Inf).
% Relax line: Location(end-1)-Location(end-2) (Ben's definition) when the path
% has >= 3 rows; for the 2-row kf sets (origin + insertion only) the fixed
% bracket2->insertion line is used instead.
% Only active for finite Xi1/Xi2 with the bracket enabled -- escape paths and
% the rigid baseline must stay finite for the driver's baseline normalization.
% Env FLXPX3_TANGENCY=DIAG -> tangency stored (d_circle2/colinear2) but NOT
% enforced: with today's 2-row paths the straight-line force direction cuts the
% circle at ~13 high-flexion frames for every candidate, so enforcement needs the
% wrapped >=3-row path from the new CAD (Ben 2026-09-14).
function chk = checkBracket2(klass, L, pB, Pbri2, e_tib2, X1, X2, useB2t)
N = size(e_tib2, 1);
chk = struct('screw2_hit', false, 'tangency_ok', true, 'reject', false, ...
             'd_circle2', nan(N,1), 'colinear2', true(N,1), 'screw2_N', zeros(N,1));
checksActive = useB2t && ~isinf(X1) && ~isinf(X2);
if ~checksActive
    return;
end
tangEnforce = ~strcmpi(getenv('FLXPX3_TANGENCY'), 'DIAG');

% 30 mm-dia circle about the knee center (tibia origin)
u = klass.unitD;                        %force unit direction, tibia frame (class-computed)
un = vecnorm(u, 2, 2);
valid = un > 1e-4;
p2def = repmat(Pbri2(:).', N, 1) + e_tib2;   %deflected bracket-2 position, tibia frame
d = vecnorm(cross(p2def, u, 2), 2, 2) ./ un; %perp. distance knee center -> force line
chk.d_circle2 = d;

% relax line for the colinear test
M = size(L, 1);
if M >= 3
    vA = squeeze(L(M-1,:,:)).';          %Location(end-1,:,:) per frame, Nx3
    vB = squeeze(L(M-2,:,:)).';          %Location(end-2,:,:) per frame, Nx3
    vref = vA - vB;                      %Ben's line: Location(end-1)-Location(end-2)
else
    vref = repmat(pB - Pbri2(:).', N, 1);   %fallback (2-row kf sets): bracket2 -> insertion
end
vn = vecnorm(vref, 2, 2);
cosang = abs(sum(u.*vref, 2)) ./ max(un.*vn, 1e-12);
chk.colinear2 = acosd(min(1, cosang)) <= 2;  %colinear within 2 deg (line, either sense)

active = valid & ~chk.colinear2;
if tangEnforce
    chk.tangency_ok = ~any(active & (d < 15e-3)); %tangent-or-outside the 15 mm radius circle
else
    chk.tangency_ok = true;                       %diagnostic mode: distances still stored
end
chk.reject = ~chk.tangency_ok;   %screw head no longer rejects -- it clamps (Ben 2026-09-15)
end