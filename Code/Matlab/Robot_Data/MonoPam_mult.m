% Pam Data -- two actual BPAs on one shared bracket
% Author: Ben Bolen
% Multi-route extension: 9/2026
% Description: Cell-aware counterpart to MonoPamDataExplicit_balanceX3.
% The two BPAs have separate routes and force directions, but use the same
% BPA parameters, tendon specification, bend measure, force magnitude, and
% shared bracket displacement. Torque is the sum from both actual routes.

classdef MonoPam_mult < handle

    %% ------------ Public Properties ---------------------------
    properties
        Name                        % BPA/muscle name
        Location                    % {2x1} BPA routing-point arrays
        Cross                       % First row in the next reference frame
        Diameter                    % BPA diameter, mm
        TransformationMat           % Reference-frame transforms
        RestingL                    % Resting length of each BPA, m
        Kmax                        % Fully contracted length of each BPA, m
        FittingLength               % Length of one fitting, m
        TendonL                     % Common tendon length for each path, m
        Pressure                    % Common BPA pressure, kPa
        AngleD                      % Joint-angle vector, degrees
        Xi0                         % Constant length offset, m
        Xi1                         % Bracket axial stiffness, N/m
        Xi2                         % Bracket transverse stiffness, N/m
        Xi3                         % Bend-loss scale factor
        BendMeasure                 % Common Nx1 sum(R*alpha), m
        Wraps                       % Common cable wrap count
        BPAcount                    % Number of actual BPA routes; currently 2
        L_p                         % {2x1} deformed routing-point arrays
        Lmt_p                       % {2x1} deformed musculotendon lengths
        delta_L                     % {2x1} bend-related usable-length losses
        strain_p                    % {2x1} contraction excluding delta_L
        strain_f                    % {2x1} contraction including delta_L
        F_p                         % {2x1} equal-magnitude BPA force vectors
        mA_p                        % {2x1} stiffness-aware moment arms
        Torque_p                    % Nx3 sum of both stiffness-aware torques
        gama                        % {2x1} equal tendon stretch histories
        kSpr                        % Stiffness of one physical tendon path
        Fmag                        % Nx1 common force magnitude per BPA, N
        c_bracket                   % Nx2 route projections of shared compliance
        c_eff                       % Nx2 bracket plus one-tendon compliance
        k_eff                       % Nx2 inverse effective compliance
        bracketDeflection           % Nx3 shared displacement in bracket frame
        forceMismatch               % Nx1 route force-prediction difference, N
    end

    %% ------------ Dependent Properties ------------------------
    properties (Dependent)
        SegmentLengths              % {2x1} route segment-length arrays
        LongestSegment              % {2x1} longest segment by route
        MuscleLength                % {2x1} complete route lengths
        Contraction                 % {2x1} rigid/original contractions
        LengthCheck                 % {2x1} usable/unusable labels
        UnitDirection               % {2x1} force directions at the joint
        MomentArm                   % {2x1} moment-arm vectors
        Fmax                        % Maximum force of one BPA, N
        Force                       % {2x1} rigid/original force vectors
        Torque                      % Nx3 summed rigid/original torque
    end

    methods
        %% ------------- Muscle Data Constructor -----------------
        function PD = MonoPam_mult(name,location,cross,diameter,t, ...
                rest,kmax,tendon,fitn,pres,xi0,xi1,xi2,xi3,wraps, ...
                angleD,bpaCount,bendMeasure)
            if nargin == 17
                bendMeasure = [];
            elseif nargin ~= 18
                error('MonoPam_mult:BadInputCount', ...
                    'Expected 17 or 18 inputs, got %d.',nargin)
            end
            if ~iscell(location) || numel(location) ~= 2
                error('MonoPam_mult:Location', ...
                    'Location must be a two-element cell array.')
            end
            if bpaCount ~= 2 || bpaCount ~= numel(location)
                error('MonoPam_mult:BPAcount', ...
                    'BPAcount and the number of Location cells must equal 2.')
            end
            if ~isscalar(tendon) || ~isscalar(wraps)
                error('MonoPam_mult:CommonTendon', ...
                    'TendonL and Wraps must be common scalar inputs.')
            end

            PD.Name = name;                   % BPA/muscle name
            PD.Location = location(:);         % two routing-point arrays
            PD.Cross = cross;                 % first row in the next frame
            PD.Diameter = diameter;           % BPA diameter, mm
            PD.TransformationMat = t;         % frame transforms
            PD.RestingL = rest;               % BPA resting length, m
            PD.Kmax = kmax;                   % fully contracted length, m
            PD.TendonL = tendon;              % tendon length on each path, m
            PD.FittingLength = fitn;           % one fitting length, m
            PD.Pressure = pres;               % BPA pressure, kPa
            PD.Xi0 = xi0;                     % constant length offset, m
            PD.Xi1 = xi1;                     % bracket axial stiffness, N/m
            PD.Xi2 = xi2;                     % transverse stiffness, N/m
            PD.Xi3 = xi3;                     % bend-loss scale factor
            PD.Wraps = wraps;                 % cable wrap count per tendon
            PD.AngleD = angleD(:);            % joint angles, degrees
            PD.BPAcount = bpaCount;            % two actual BPA routes
            PD.BendMeasure = bendMeasure;      % common R*angle history, m

            N = numel(PD.AngleD);
            for j = 1:PD.BPAcount
                L = PD.Location{j};
                if size(L,2) ~= 3 || size(L,3) ~= N
                    error('MonoPam_mult:LocationSize', ...
                        'Location{%d} must be M-by-3-by-%d.',j,N)
                end
            end
            if size(PD.Location{1},1) ~= size(PD.Location{2},1)
                error('MonoPam_mult:RouteSize', ...
                    'Both Location arrays need the same number of rows.')
            end
            if size(t,3) ~= N
                error('MonoPam_mult:TransformSize', ...
                    'TransformationMat needs one page per joint position.')
            end
            if ~isempty(bendMeasure) && numel(bendMeasure) ~= N
                error('MonoPam_mult:BendMeasureSize', ...
                    'BendMeasure needs one value per joint position.')
            end
            PD = PD.updateStiffnessGeometry();
        end

        %% ------------- Segment Lengths ------------------------
        function segLengths = get.SegmentLengths(obj)
            T = obj.TransformationMat;
            C = obj.Cross;
            N = size(T,3);
            segLengths = cell(obj.BPAcount,1);
            for j = 1:obj.BPAcount
                L = obj.Location{j};
                routeLengths = zeros(N,size(L,1)-1);
                for ii = 1:N
                    for i = 1:size(L,1)-1
                        pointA = L(i,:,ii);
                        pointB = L(i+1,:,ii);
                        if i+1 == C
                            pointB = RowVecTrans(T(:,:,ii),pointB);
                        end
                        routeLengths(ii,i) = norm(pointA-pointB);
                    end
                end
                segLengths{j} = routeLengths;
            end
        end

        %% -------------- Longest Segment Calculation -----------
        function longestSeg = get.LongestSegment(obj)
            segLengths = obj.SegmentLengths;
            longestSeg = cellfun( ...
                @(x)x(:,find(mean(x,1)==max(mean(x,1)),1,'first')), ...
                segLengths,'UniformOutput',false);
        end

        %% ------------- Muscle Length --------------------------
        function mL = get.MuscleLength(obj)
            segLengths = obj.SegmentLengths;
            mL = cellfun(@(x)sum(x,2),segLengths,'UniformOutput',false);
        end

        %% -------------- Force Unit Direction ------------------
        function unitD = get.UnitDirection(obj)
            T = obj.TransformationMat;
            C = obj.Cross;
            N = size(T,3);
            unitD = cell(obj.BPAcount,1);
            for j = 1:obj.BPAcount
                L = obj.Location{j};
                direction = zeros(N,3);
                for i = 1:N
                    pointA = L(C-1,:,i);
                    pointB = L(C,:,i);
                    direction(i,:) = ...
                        RowVecTrans(T(:,:,i)\eye(4),pointA)-pointB;
                end
                unitD{j} = normalize(direction);
            end
        end

        %% -------------- Moment Arm ----------------------------
        function mA = get.MomentArm(obj)
            C = obj.Cross;
            N = size(obj.TransformationMat,3);
            unitD = obj.UnitDirection;
            mA = cell(obj.BPAcount,1);
            for j = 1:obj.BPAcount
                L = obj.Location{j};
                routeMomentArm = zeros(N,3);
                for i = 1:N
                    pointB = L(C,:,i);
                    routeMomentArm(i,:) = pointB-unitD{j}(i,:)* ...
                        dot(unitD{j}(i,:),pointB);
                end
                mA{j} = routeMomentArm;
            end
        end

        %% -------------- Contraction of the PAM ----------------
        function contraction = get.Contraction(obj)
            mL = obj.MuscleLength;
            contraction = cell(obj.BPAcount,1);
            for j = 1:obj.BPAcount
                contraction{j} = (obj.RestingL-(mL{j}-obj.TendonL- ...
                    2*obj.FittingLength))/obj.RestingL;
            end
        end

        %% -------------- Length Check --------------------------
        function lengthCheck = get.LengthCheck(obj)
            contraction = obj.Contraction;
            lengthCheck = cell(obj.BPAcount,1);
            for j = 1:obj.BPAcount
                if obj.RestingL >= 0 && max(contraction{j}) <= 0.25 && ...
                        min(contraction{j}) >= -0.1
                    lengthCheck{j} = 'Usable';
                else
                    lengthCheck{j} = 'Unusable';
                end
            end
        end

        %% -------------- Maximum Force -------------------------
        function maxF = get.Fmax(obj)
            if obj.Diameter == 10
                maxF = maxBPAforce(obj.RestingL,620);
            elseif obj.Diameter == 20
                maxF = maxBPAforce(obj.RestingL,'20',620);
            elseif obj.Diameter == 40
                maxF = 6000;
            else
                error('MonoPam_mult:Diameter','Unsupported BPA diameter.')
            end
        end

        %% -------------- Force ---------------------------------
        function F = get.Force(obj)
            unitD = obj.UnitDirection;
            strain = obj.Contraction;
            KMAX = (obj.RestingL-obj.Kmax)/obj.RestingL;
            routeForce = zeros(numel(obj.AngleD),obj.BPAcount);
            for j = 1:obj.BPAcount
                Fn = festo4(obj.Diameter,strain{j}/KMAX,obj.Pressure);
                routeForce(:,j) = max(0,Fn.*obj.Fmax);
            end
            commonForce = mean(routeForce,2);
            F = cell(obj.BPAcount,1);
            for j = 1:obj.BPAcount
                F{j} = commonForce.*unitD{j};
            end
        end

        %% ---------------------- Torque ------------------------
        function tor = get.Torque(obj)
            mA = obj.MomentArm;
            F = obj.Force;
            tor = zeros(size(F{1}));
            for j = 1:obj.BPAcount
                tor = tor+cross(mA{j},F{j},2);
            end
        end

        %% ======================================================
        %% Stiffness-aware pipeline (X3 structure, two routes)
        %% ======================================================
        function obj = updateStiffnessGeometry(obj)
            % kSpr is one physical tendon stiffness; do not multiply by two.
            obj.kSpr = Spr(obj,obj.Wraps);
            Funit = computeForceVector(obj);
            [strain_Xi3,delta_L_i] = Contraction_k( ...
                obj,[],obj.Xi0,[],obj.Xi3);
            obj.delta_L = delta_L_i;

            deltaAndXi0 = cellfun(@(x)x+obj.Xi0,delta_L_i, ...
                'UniformOutput',false);
            [L_p_i,gama_i,Fmag_i,cBracket_i,cEff_i,dBr_i,mismatch_i] = ...
                Lok(obj,obj.Xi1,obj.Xi2,obj.kSpr,Funit, ...
                strain_Xi3,deltaAndXi0);
            obj.L_p = L_p_i;
            obj.gama = gama_i;
            obj.Fmag = Fmag_i;
            obj.c_bracket = cBracket_i;
            obj.c_eff = cEff_i;
            obj.k_eff = 1./cEff_i;
            obj.bracketDeflection = dBr_i;
            obj.forceMismatch = mismatch_i;

            % Segment lengths after deformation; deliberately inline.
            T = obj.TransformationMat;
            C = obj.Cross;
            N = size(T,3);
            sL_p = cell(obj.BPAcount,1);
            for j = 1:obj.BPAcount
                L = obj.L_p{j};
                routeLengths = zeros(N,size(L,1)-1);
                for ii = 1:N
                    for i = 1:size(L,1)-1
                        pointA = L(i,:,ii);
                        pointB = L(i+1,:,ii);
                        if i+1 == C
                            pointB = RowVecTrans(T(:,:,ii),pointB);
                        end
                        routeLengths(ii,i) = norm(pointA-pointB);
                    end
                end
                sL_p{j} = routeLengths;
            end

            Lmt_p_i = cellfun(@(x)LMT(x,obj.Xi0),sL_p, ...
                'UniformOutput',false);
            obj.Lmt_p = Lmt_p_i;

            [strain_f_i,~] = Contraction_k( ...
                obj,Lmt_p_i,[],obj.gama,obj.Xi3);
            obj.strain_f = strain_f_i;

            % Unit directions after deformation; deliberately inline.
            unitD_p = cell(obj.BPAcount,1);
            for j = 1:obj.BPAcount
                L = obj.L_p{j};
                direction = zeros(N,3);
                for i = 1:N
                    pointA = L(C-1,:,i);
                    pointB = L(C,:,i);
                    direction(i,:) = ...
                        RowVecTrans(T(:,:,i)\eye(4),pointA)-pointB;
                end
                unitD_p{j} = normalize(direction);
            end

            [F_p_i,~,finalMismatch] = Force_p( ...
                obj,unitD_p,obj.strain_f,obj.Fmag);
            obj.F_p = F_p_i;
            obj.forceMismatch = finalMismatch;

            % Moment arms after deformation; deliberately inline.
            mA_p_i = cell(obj.BPAcount,1);
            for j = 1:obj.BPAcount
                L = obj.L_p{j};
                routeMomentArm = zeros(N,3);
                for i = 1:N
                    pointB = L(C,:,i);
                    routeMomentArm(i,:) = pointB-unitD_p{j}(i,:)* ...
                        dot(unitD_p{j}(i,:),pointB);
                end
                mA_p_i{j} = routeMomentArm;
            end
            obj.mA_p = mA_p_i;

            [strain_p_i,~] = Contraction_k( ...
                obj,Lmt_p_i,[],obj.gama,[]);
            obj.strain_p = strain_p_i;
            obj.Torque_p = Tor(obj,obj.mA_p,obj.F_p,obj.strain_p);
        end
    end
end

%% =====================================================================
%% Helper functions (same order and roles as the X3 class)
%% =====================================================================

%% ------------- Force unit direction -------------------------
function F_unit = computeForceVector(obj)
T = obj.TransformationMat;
C = obj.Cross;
N = size(T,3);
F_unit = cell(obj.BPAcount,1);
for j = 1:obj.BPAcount
    L = obj.Location{j};
    pt1 = zeros(N,3);
    pt2 = NaN(N,3);
    for i = 1:N
        pt1(i,:) = L(1,:,i);
        for k = 2:size(L,1)
            if norm(L(k,:,i)-L(1,:,i)) > 1e-6
                pt2(i,:) = L(k,:,i);
                if k == C
                    pt2(i,:) = RowVecTrans(T(:,:,i),pt2(i,:));
                end
                break
            end
        end
        if any(isnan(pt2(i,:)))
            warning('MonoPam_mult:RepeatedRoute', ...
                'Route %d, frame %d has no valid second point.',j,i)
            pt2(i,:) = pt1(i,:);
        end
    end
    F_unit{j} = normalize(pt2-pt1);
end
end

%% -------------- Contraction of the PAM ----------------------
function [contraction,delta_L] = Contraction_k(obj,Lmt,X0,gama,X3)
rest = obj.RestingL;
tendon = obj.TendonL;
fitn = obj.FittingLength;
theta_k = obj.AngleD(:);
N = numel(theta_k);
KMAX = (rest-obj.Kmax)/rest;
if isempty(Lmt)
    Lmt = obj.MuscleLength;
end
if isempty(X0)
    X0 = 0;
end
if isempty(gama)
    gama = repmat({zeros(N,1)},obj.BPAcount,1);
elseif ~iscell(gama)
    gama = repmat({gama(:)},obj.BPAcount,1);
end

contraction = cell(obj.BPAcount,1);
delta_L = cell(obj.BPAcount,1);
for j = 1:obj.BPAcount
    deltaRoute = zeros(N,1);
    if ~isempty(X3)
        Lm0 = Lmt{j}-tendon-2*fitn-X0-gama{j};
        strain0 = (rest-Lm0)/rest;
        comp = max(0,1-strain0/KMAX);
        if ~isempty(obj.BendMeasure)
            bendMeasure = obj.BendMeasure(:);
            deltaRoute = X3.*bendMeasure.*comp.^2;
        else
            ang = -9.19;
            angleRad = deg2rad((ang-theta_k)*80/(ang+120));
            idx = angleRad > 0;
            deltaRoute = X3*0.022*deg2rad(28).*comp.^2;
            deltaRoute(idx) = deltaRoute(idx)+ ...
                X3*0.176.*angleRad(idx).*comp(idx).^2;
        end
    end
    Lm_adj = Lmt{j}-tendon-2*fitn-X0-gama{j}-deltaRoute;
    contraction{j} = (rest-Lm_adj)/rest;
    delta_L{j} = deltaRoute;
end
end

%% ------------- Location ------------------------------------
function [LOC,gema,Fmag,cBracket,cEff,dBr,forceMismatch] = ...
        Lok(obj,X1,X2,kSpr,Funit,strain_predef,deltaL)
L = obj.Location;
N = size(obj.TransformationMat,3);
iRef = min(92,N);
pA = (L{1}(1,:,iRef)+L{2}(1,:,iRef))/2;
switch obj.Diameter
    case 20
        Pbr = [9.48,-33.38,30.86]/1000;
    case 10
        Pbr = [-19,22,27.6]/1000;
    otherwise
        Pbr = [0,0,0];
end
phbrA = pA-Pbr;
thetaZ = atan2(phbrA(2),phbrA(1));
RhbrZ = [cos(thetaZ),-sin(thetaZ),0; ...
          sin(thetaZ), cos(thetaZ),0; 0,0,1];
pbrhA = RhbrZ'*phbrA';
thetaY = atan2(pbrhA(3),pbrhA(1));
Ry = [cos(thetaY),0,sin(thetaY); 0,1,0; ...
     -sin(thetaY),0,cos(thetaY)];
Rhbr = RhbrZ*Ry';

[dBr,commonGama,Fmag,cBracket,cEff,forceMismatch] = ...
    fortz(obj,Funit,Rhbr,X1,X2,kSpr,deltaL,strain_predef);

LOC = L;
for j = 1:obj.BPAcount
    for ii = 1:N
        dBody = dBr(ii,:)*Rhbr';
        originalPoint = L{j}(1,:,ii);
        LOC{j}(1,:,ii) = originalPoint+dBody;
        for jj = 2:obj.Cross-1
            if norm(L{j}(jj,:,ii)-originalPoint) < 1e-10
                LOC{j}(jj,:,ii) = L{j}(jj,:,ii)+dBody;
            else
                break
            end
        end
    end
end
gema = repmat({commonGama},obj.BPAcount,1);
end

%% ------------- Shared force and length equilibrium ----------
function [dBr,e_cable,Fmag,cBracket,cEff,forceMismatch] = ...
        fortz(obj,Funit,Rhbr,X1,X2,kSpr,deltaL,strain_predef)
% fzero solves one common BPA force at each joint position. The shared
% bracket load is F*(u1+u2); each independent tendon carries F.
N = numel(obj.AngleD);
rest = obj.RestingL;
tendon = obj.TendonL;
fitn = obj.FittingLength;
KMAX = (rest-obj.Kmax)/rest;
Fm = obj.Fmax;
D = obj.Diameter;
P = obj.Pressure;
if isinf(X1), c1 = 0; else, c1 = 1/X1; end
if isinf(X2), c2 = 0; else, c2 = 1/X2; end
Cbr = diag([c1,c2,c2]);
if isinf(kSpr), cCable = 0; else, cCable = 1/kSpr; end

uBr = cell(obj.BPAcount,1);
for j = 1:obj.BPAcount
    uBr{j} = Funit{j}*Rhbr;
end
uSum = uBr{1}+uBr{2};
cBracket = zeros(N,obj.BPAcount);
for j = 1:obj.BPAcount
    cBracket(:,j) = sum((uSum*Cbr).*uBr{j},2);
end
cEff = cBracket+cCable;

mL = obj.MuscleLength;
for j = 1:obj.BPAcount
    mL{j} = mL{j}-deltaL{j}-tendon-2*fitn;
end
mL1 = mL{1};
mL2 = mL{2};
Fmag = zeros(N,1);
dBr = zeros(N,3);
e_cable = zeros(N,1);
forceMismatch = zeros(N,1);
initialForce = zeros(N,obj.BPAcount);
for j = 1:obj.BPAcount
    initialForce(:,j) = max(0, ...
        festo4(D,strain_predef{j}/KMAX,P).*Fm);
end

% Only the independent, per-angle fzero solves are parallelized.
parfor i = 1:N
    Lm1 = mL1(i);
    Lm2 = mL2(i);
    ce1 = cEff(i,1);
    ce2 = cEff(i,2);
    relfun = @(F) F-0.5*Fm*( ...
        max(0,festo4(D,(rest-(Lm1-F*ce1))/rest/KMAX,P))+ ...
        max(0,festo4(D,(rest-(Lm2-F*ce2))/rest/KMAX,P)) );
    fLo = relfun(0);
    fHi = relfun(Fm);
    if isfinite(fLo) && isfinite(fHi) && sign(fLo) ~= sign(fHi)
        try
            commonForce = fzero(relfun,[0,Fm]);
        catch
            commonForce = mean(initialForce(i,:));
        end
    else
        commonForce = mean(initialForce(i,:));
    end
    commonForce = min(Fm,max(0,commonForce));
    pred1 = Fm*max(0,festo4(D, ...
        (rest-(Lm1-commonForce*ce1))/rest/KMAX,P));
    pred2 = Fm*max(0,festo4(D, ...
        (rest-(Lm2-commonForce*ce2))/rest/KMAX,P));
    Fmag(i) = commonForce;
    dBr(i,:) = commonForce*(uSum(i,:)*Cbr);
    e_cable(i) = commonForce*cCable;
    forceMismatch(i) = abs(pred1-pred2);
end
end

%% ------------- Muscle Length -------------------------------
function Lmt = LMT(sL,X0)
Lmt = sum(sL,2);
if ~isempty(X0), Lmt = Lmt-X0; end
end

%% -------------- Force --------------------------------------
function [F,commonForce,forceMismatch] = ...
        Force_p(obj,unitD_p,strain,solvedForce)
KMAX = (obj.RestingL-obj.Kmax)/obj.RestingL;
N = numel(obj.AngleD);
routeForce = zeros(N,obj.BPAcount);
for j = 1:obj.BPAcount
    Fn = festo4(obj.Diameter,strain{j}/KMAX,obj.Pressure);
    routeForce(:,j) = max(0,Fn.*obj.Fmax);
end
forceMismatch = abs(routeForce(:,1)-routeForce(:,2));
if nargin < 4 || isempty(solvedForce)
    commonForce = mean(routeForce,2);
else
    commonForce = solvedForce(:);
end
F = cell(obj.BPAcount,1);
for j = 1:obj.BPAcount
    F{j} = commonForce.*unitD_p{j};
end
end

%% ---------------------- Torque ------------------------------
function Mz = Tor(obj,mA_p,F_p,strain_p)
Mz = zeros(size(F_p{1}));
if obj.Diameter == 20, ss = -0.03; else, ss = -0.02; end
for j = 1:obj.BPAcount
    routeTorque = cross(mA_p{j},F_p{j},2);
    routeTorque(strain_p{j} < ss,:) = NaN;
    Mz = Mz+routeTorque;
end
end

%% -------------- Tendon spring rate -------------------------
function springrate = Spr(obj,wraps)
if obj.TendonL == 0
    springrate = Inf;
    return
end
if obj.Diameter == 20, mult = 6; else, mult = 2; end
if nargin >= 2 && ~isempty(wraps), mult = mult*wraps; end
Aeff = 1.51e-6;
E = 193e9;
springrate = mult*Aeff*E/obj.TendonL;
end

%% -------------- Normalize row vectors ----------------------
function vhat = normalize(v)
norms = vecnorm(v,2,2);
valid = norms > 1e-3 & all(~isnan(v),2);
vhat = zeros(size(v));
vhat(valid,:) = v(valid,:)./norms(valid);
end
