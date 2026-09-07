%% Two-bracket flexor evaluator (new method, Sept 2026)
function [f_all, bpa_all] = minimizeFlxPin2brk(Xi0,Xi1,Xi2,idx_val,useB2)
% minimizeFlxPin2brk: minimizeFlxPin with bracket-model changes (d), (e), (f):
%   (d) insertion bracket uses the single (yaw-only) rotation:
%         Tkbr = RpToTrans(RkbrZ, Pbri')   [was RpToTrans(Rkbr, Pbri')]
%         pbrBnew = [norm(pkbrB(1:2)), 0, pkbrB(3)] + eB   [was [norm(pkbrB), 0, 0] + eB]
%   (e) insertion-bracket stiffness array  K  = [X1, X2, X1]   [was [X1, X2, X2]]
%   (f) SECOND (origin-side) bracket at the hip frame, to capture movement at pA
%       that the single-bracket approximation misses at high force:
%         computeForceVector ports the minimizeExtX3 pattern: hip-frame unit
%         direction from Loc geometry (origin -> first distinct point).
%         Pbr2 = [-52.61, 0, 75.06]/1000 (hip origin to flexor origin bracket),
%         yaw-only bracket frame Thbr = RpToTrans(RhbrZ, Pbr2'),
%         K2 = [X2, X1, X2], compliance projected on the path force direction,
%         c_eff = c_b + cSpr + c_b2, pbrAnew = [norm(phbrA(1:2)), 0, phbrA(3)] + eA,
%         LOC row 1 (origin) updated with pAnew.
%   useB2 (optional, default true): false reproduces the (d)+(e)-only variant
%   (bracket-2 compliance and deflection zeroed) for clean ablation.
%
% Inputs:
%   Xi0 - extra length correction (m)
%   Xi1 - bracket "axial" stiffness (N/m)
%   Xi2 - bracket "bending" stiffness (N/m)
%   idx_val - indices of the BPAs to process (default: all)
%   useB2 - enable the second (origin-side) bracket (default: true)
%
% Outputs:
%   f_all   - [RMSE, FVU, MaxResidual] per BPA
%   bpa_all - updated BPA structs (bpa_all(i).eA2 holds origin-bracket deflections)

%% load
load FlxPinBPASet.mat kf

%Declare the new diagnostics field so struct assignment in the evaluate loop
%stays field-order compatible with the saved kf template.
for jD = 1:numel(kf)
    kf(jD).eA2 = [];
end

%% Initialize output
nBPA = numel(kf);
if nargin < 4 || isempty(idx_val)
        idx_val = 1:nBPA;
end
if nargin < 5
    useB2 = true;
end

bpa_all = kf;
f_all = NaN(nBPA, 3);

%% Evaluate each BPA
for i = idx_val
    klass_i = kf(i);
    [bpa_all(i), f_all(i,:)] = evaluateBPA(klass_i, Xi0, Xi1, Xi2, useB2);
    if any(isnan(bpa_all(i).strain_p))
        warning('NaNs in strain_p for BPA #%d', i);
    end
end

end


function [bpa_i, fitvec] = evaluateBPA(klass, Xi0, Xi1, Xi2, useB2)
%% Calculate locations and properties
bpa_i = klass;
kspr = Spr(bpa_i); %Calculate spring rate (Infinite if no tendon is used)
strain_Xi0 = Contraction(bpa_i, [],Xi0); %Calculate contraction with constant length offset
Funit = computeForceVector(bpa_i);   %Force unit direction at the origin, hip frame (bracket 2)
[L_p, gemma, eA2] = Lok(bpa_i, Xi1, Xi2, kspr, strain_Xi0, Xi0, useB2, Funit);   %Bracket deformation changing geometry
unitD_p = UD(bpa_i, L_p);   %New force direction
sL_p = seg(bpa_i, L_p);   %New segment lengths uses deformation but does not subtract length offset
Lmt_p = LMT(sL_p, Xi0+gemma);     %New musclulotendon length. Uses deformed geometry and constant length offset.
strain_p = Contraction(bpa_i, Lmt_p, []);  %new contraction amount includes deformed geometry and constant length offset
F_p = Force(bpa_i, unitD_p, strain_p);  %new force vector
mA_p = Mom(bpa_i, L_p, unitD_p);   %new moment arm
M_p = Tor(mA_p, F_p, strain_p);  %new torque

%% Package into output struct
bpa_i.Lmt_p = Lmt_p;
bpa_i.mA_p = mA_p;
bpa_i.M_p = M_p;
bpa_i.F_p = F_p;
bpa_i.strain_p = strain_p;
bpa_i.L_p = L_p;
bpa_i.gama = gemma;
bpa_i.eA2 = eA2;    %origin-bracket deflections [axial2, bendY2, bendZ2], N x 3

% GoF calculation
fitvec = SSE(bpa_i, M_p);

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

%% ------------- Force unit direction (hip frame, bracket 2) ------
function F_unit = computeForceVector(klass)
%Port of minimizeExtX3.computeForceVector: unit direction from the muscle
%origin (hip frame) to the next real point, per frame. Takes into account
%via points and the hip->tibia frame change at the cross point.
L = klass.Loc;      %Location (wrapping, attachment points)
C = klass.CP;       %Cross point (moves from one frame to another)
T = klass.Tk;       %Transformation matrix

N = size(L, 3);
pt1 = NaN(N, 3);
pt2 = NaN(N, 3);

for i = 1:N
    pt1(i,:) = L(1,:,i);  % muscle origin
    found = false;
    for k = 2:size(L,1)
        if ~isequal(L(k,:,i),L(1,:,i))
            if k == C
                pt2(i,:) = RowVecTrans(T(:,:,i), L(k,:,i));
            else
                pt2(i,:) = L(k,:,i);
            end
            found = true;
            break;
        end
    end
    if ~found
        warning("computeForceVector: Frame %d has no valid second point, using pt1=pt2", i);
        pt2(i,:) = pt1(i,:);
    end
end

F_vec = pt2 - pt1;
F_unit = normalize(F_vec);
end

%% ------------- Location  ------------------------
function [LOC, gama, eA2] = Lok(klass,X1,X2,kSpr,strain,X0,useB2,Funit)
% Inputs:
%   bpa class info
%   X1, X2 stiffness
%   kSpr, tendon stiffness
%   strain – N×1 strain vector (e.g., from Xi0 offset effect)
%   X0, constant length offset
%   useB2, enable second (origin-side) bracket
%   Funit, N×3 force unit direction at the origin, hip frame (for bracket 2)
L = klass.Loc;      %Location (wrapping, attachment points)
C = klass.CP;       %Cross point (moves from one frame to another)
kmax = klass.Kmax;  %max contracted length
KMAX = (klass.rest-kmax)/klass.rest; %turn it into a percentage

relstrain = strain/KMAX;            %relative strain
FF = festo4(klass.dBPA,relstrain,klass.P).*klass.Fm;        %Force magnitude
FF(FF<0) = 0;
%For muscle insertion (bracket 1)
unitD = klass.unitD;            %unit direction of force vector, tibia frame
Fk = unitD.*FF;                  %Force vector, tibia frame
Fh = Funit.*FF;                  %Force vector at the origin, hip frame (bracket 2)
pB = L(C,:,(klass.Ak==0));                  %Distance from knee frame to muscle insertion
% Pbri = [-48.11, -107.81, 13.8]/1000;     %vector from knee ICR to flexor insertion bracket (where it starts to cantilever)
Pbri = [-27.5, -107.81, -0.54]/1000;     %vector from knee ICR to flexor insertion bracket (where it starts to cantilever, but at tibial contact, no z offset)
% Pbri = [-27.5, -125.91, -0.54]/1000;     %vector from knee ICR to upper bolt
pkbrB = pB-Pbri;                  %vector from bracket to point B, in the knee frame
thetabrB = atan2(pkbrB(2),pkbrB(1));   %angle between pbrB and x axis
RkbrZ = [cos(thetabrB) -sin(thetabrB) 0; ...     %Rotation matrix
       sin(thetabrB) cos(thetabrB) 0; ...
       0    0   1];
% (d) single (yaw-only) rotation for the insertion bracket
Tkbr = RpToTrans(RkbrZ, Pbri');    %Transformation matrix, flexor bracket frame in knee frame

% (f) second bracket at the origin (hip frame)
pA = L(1,:,1);                                 %Distance from hip frame to muscle origin
Pbr2 = [-52.61, 0, 75.06]/1000;                 %vector from hip origin to flexor origin bracket
phbrA = pA-Pbr2;                               %vector from bracket to point A (in the hip frame)
thetabrA = atan2(phbrA(2),phbrA(1));           %angle between phbrA and x axis
RhbrZ = [cos(thetabrA) -sin(thetabrA) 0; ...   %Rotation matrix
       sin(thetabrA) cos(thetabrA) 0; ...
       0    0   1];
Thbr = RpToTrans(RhbrZ, Pbr2');    %Transformation matrix, origin bracket frame in hip frame

LOC = L;            %new location matrix
N = size(L,3);
Fbrk = zeros(N,3);       %Force vector represented in the tibial bracket frame
Fbrh2 = zeros(N,3);      %Force vector represented in the origin bracket frame

parfor ii = 1:N                          %Repeat for each orientation
        Fbrk(ii,:) = RowVecTrans(Tkbr\eye(4),Fk(ii,:)); %Force vector in the tibia frame represented in the insertion bracket frame
    if useB2
        Fbrh2(ii,:) = RowVecTrans(Thbr\eye(4),Fh(ii,:)); %Force vector in the hip frame represented in the origin bracket frame
    end
end

if isinf(X1) && isinf(X2)  && isinf(kSpr)
    [epsilon, delta, beta, gama, e_ax2, e_by2, e_bz2] = deal(zeros(N,1));
else
    [epsilon, delta, beta, gama, e_ax2, e_by2, e_bz2] = fortz(klass,Fbrk,Fbrh2,X1,X2,kSpr,X0,useB2);
end

eB = [epsilon, delta, beta];
pbrBnew = [norm(pkbrB(1:2)), 0, pkbrB(3)] + eB; %new point B, in the insertion bracket's frame (d)

pBnew = zeros(N,3);
for ii = 1:N                          %Repeat for each orientation
    pBnew(ii,:) = RowVecTrans(Tkbr, pbrBnew(ii,:));     %New point B, in the tibia frame
    LOC(2,:,ii) = pBnew(ii,:);
end

% (f) displaced muscle origin, in the origin bracket frame
eA = [e_ax2, e_by2, e_bz2];
eA2 = eA;
if useB2
    pbrAnew = [norm(phbrA(1:2)), 0, phbrA(3)] + eA; %Muscle origin location, bracket frame
    pAnew = zeros(N,3);
    for ii = 1:N
        pAnew(ii,:) = RowVecTrans(Thbr, pbrAnew(ii,:)); %New point A, in the hip frame
        LOC(1,:,ii) = pAnew(ii,:);
        for k = 2:C-1        %replace any wrapping points that repeat the origin
            if isequal(L(k,:,ii),L(1,:,ii))
               LOC(k,:,ii) = pAnew(ii,:);
            end
        end
    end
end
end

%% Force and length reduction due to deformation
function [e_axial, e_bendY, e_bendZ, e_cable, e_axial2, e_bendY2, e_bendZ2] = fortz(klass,Fbr,Fbr2,X1,X2,kSpr,X0,useB2)
% e_axial, insertion-bracket axial elongation
% e_bendY, insertion-bracket bending displacement y - direction
% e_bendZ, insertion-bracket bending displacement z - direction
% e_cable, tendon cable stretch
% e_axial2/e_bendY2/e_bendZ2, origin-bracket deflections (bracket 2)
    N = size(Fbr,1);
    % Initialize outputs
    [e_axial, e_bendY, e_bendZ, e_cable, e_axial2, e_bendY2, e_bendZ2] = deal(zeros(N,1));

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

    % Normalize force vectors safely (insertion bracket)
    norms = vecnorm(Fbr, 2, 2);
    valid = norms > 1e-4 & all(~isnan(Fbr), 2);
    u_hat_all = normalize(Fbr);

    % Vectorized k_b computation (insertion bracket, (e) axis order)
    K = [X1, X2, X1];   %bracket stiffness array
    K_bracket = diag(K);       %bracket stiffness matrix
    C_bracket = diag([1/K(1), 1/K(2), 1/K(3)]); %bracket compliance
    u_hat = permute(u_hat_all, [3, 2, 1]);  % [1x3xN]
    C_rep = repmat(C_bracket, [1, 1, N]);   % [3x3xN]
    c_b = pagemtimes(pagemtimes(u_hat, C_rep), permute(u_hat, [2, 1, 3]));
    c_b = reshape(c_b, [N, 1]);
    cSpr = 1/kSpr;

    % (f) origin bracket projected on the path force direction
    if useB2
        norms2 = vecnorm(Fbr2, 2, 2);
        valid2 = norms2 > 1e-4 & all(~isnan(Fbr2), 2);
        u_hat_all2 = normalize(Fbr2);
        K2 = [X2, X1, X2]; %bracket stiffness array, origin bracket
        K_bracket2 = diag(K2);       %bracket stiffness matrix
        C_bracket2 = diag([1/K2(1), 1/K2(2), 1/K2(3)]); %bracket compliance
        u_hat2 = permute(u_hat_all2, [3, 2, 1]);  % [1x3xN]
        C_rep2 = repmat(C_bracket2, [1, 1, N]);   % [3x3xN]
        c_b2 = pagemtimes(pagemtimes(u_hat2, C_rep2), permute(u_hat2, [2, 1, 3]));
        c_b2 = reshape(c_b2, [N, 1]);
        valid = valid & valid2;  %need a sane force direction at both brackets
    else
        u_hat_all2 = zeros(N, 3);
        K_bracket2 = eye(3);
        c_b2 = zeros(N, 1);
    end
    c_eff = c_b + cSpr + c_b2;        %effective compliance (both brackets + tendon)
    k_eff = 1 ./ c_eff;         % effective stiffness along u


    % Parallel root solve
   parfor i = 1:N
        if ~valid(i)
            continue;
        end

        % Per-instance constants
        keff = k_eff(i);
        unit_vec = u_hat_all(i, :);
        unit_vec2 = u_hat_all2(i, :);
        Lm = mL(i);

        contraction0    = ( rest - Lm ) / rest;
        relstrain0      = contraction0 / KMAX;  %relative strain
        if relstrain0 >= 1
            r = 0;
        else
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

            % Insertion-bracket displacement
            e_bkt = K_bracket \ (F_mag * unit_vec');

            e_axial(i) = e_bkt(1);
            e_bendY(i) = e_bkt(2);
            e_bendZ(i) = e_bkt(3);

            % (f) Origin-bracket displacement
            e_bkt2 = zeros(3,1);
            if useB2
                e_bkt2 = K_bracket2 \ (F_mag * unit_vec2');
                e_axial2(i) = e_bkt2(1);
                e_bendY2(i) = e_bkt2(2);
                e_bendZ2(i) = e_bkt2(3);
            end

            % Cable elongation
                if tendon > 0
                    r_bracket = unit_vec * e_bkt;
                    if useB2
                        r_bracket = r_bracket + unit_vec2 * e_bkt2;
                    end
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

            for i = 1:size(T, 3)
                pointA = L(C-1, :, i);
                pointB = L(C, :, i);
                direction(i, :) = RowVecTrans(T(:, :, i)\eye(4), pointA) - pointB;
            end
            unitD = normalize(direction);
end

        %% -------------- Moment Arm --------------------------
        %Calculate the moment arm about a joint
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
           rest = klass.rest;
           kmax = klass.Kmax;
           KMAX = (rest-kmax)/rest; %turn it into a percentage

           rel = strain./KMAX;                    %relative strain

           Fn = festo4(klass.dBPA,rel,klass.P);

           scalarForce = Fn.*klass.Fm;
           scalarForce(scalarForce <= 0) = 0;

            F = scalarForce.*unitD;

end

        %% ---------------------- Torque --------------
function Mz = Tor(mA, F, strain)
N = size(F, 1);
Mz = zeros(N, 3);

    for i = 1:N
        if strain(i,:) < -0.03
            Mz(i,:) = NaN;
        else
            Mz(i,:) = cross(mA(i,:), F(i,:));
        end
    end
end

%% tendon springrate
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
