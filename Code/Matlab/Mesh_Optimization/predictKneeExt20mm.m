function pred = predictKneeExt20mm(x, ctx)
% predictKneeExt20mm
%
% Builds the 20 mm extensor distal-ring path, evaluates the X3 stiffness
% model, and returns torque/strain/path-length diagnostics.
%
% x = [p1(1:3), p8(1:3), rest, tendon]
%
% p1 is a 1x3 point in the femur/parent frame.
% p8 is a 1x3 point in the theta1/t1 frame.
%
% buildDistalRingLocation20mm constructs p1-p8. Points p5-p8 are
% transformed from the t1 frame to the ICR frame using T_ICR_t1(:,:,i)
% before they are stored in Location(:,:,i).
% MonoPamDataExplicit_balanceX3 uses T_Pam(:,:,i) for the crossing segment.

    p1       = x(1:3);       % parent-frame / femur-side attachment
    p8       = x(4:6);       % theta1-frame insertion design variable
    rest     = x(7);         % active BPA rest length
    tendon   = x(8);         % physical tendon length
    
    % Fixed identified stiffness parameters
    Xi0 = ctx.Xi0;
    Xi1 = ctx.Xi1;
    Xi2 = ctx.Xi2;
    Xi3 = ctx.Xi3;

    BPAcount = ctx.BPAcount;

    KMAX = ctx.KMAX;
    kmax = rest*(1-KMAX);      % measured free-contracted BPA length at 620 kPa

    pred.ok = true;
    pred.failReason = "";

    if rest <= 0 || kmax <= 0 || kmax >= rest || tendon < 0
        pred.ok = false;
        pred.failReason = "Invalid length parameters";
        return
    end
    %% Build distal-ring Location matrix
    [Location, bendMeasure, routeInfo] = buildDistalRingLocation20mm( p1, p8, tendon, ctx);

    %% Normal, undeformed path length calculation
    % This is useful for sanity checking and is the same length calculation
    % used for the initial guess in buildKneeExtContext20mm.
    Lmt0 = muscleLengthNormal(Location, ctx.CrossPoint, ctx.T_Pam);


   %% Evaluate X3/stiffness-aware class

    try
        bpa = MonoPamDataExplicit_balanceX3( ...
                ctx.Name, ...
                Location, ...
                ctx.CrossPoint, ...
                ctx.Dia, ...
                ctx.T_Pam, ...
                rest, ...
                kmax, ...
                tendon, ...
                ctx.fitting, ...
                ctx.targetPressure, ...
                Xi0, ...
                Xi1, ...
                Xi2, ...
                Xi3, ...
                ctx.wraps, ...
                ctx.phiD, ...
                BPAcount, ...
                bendMeasure);
    catch ME
        pred.ok = false;
        pred.failReason = string(ME.message);
        return
    end

        %% Store model outputs
    pred.Location = Location;
    pred.bpa = bpa;

    pred.Torque  = bpa.Torque_p;
    pred.TorqueZ = bpa.Torque_p(:,3);

    % strain_f includes Xi3 and is used for force.
    pred.strain_f = bpa.strain_f(:);

    % strain_p excludes Xi3 and is used for measured-strain / stretch checks.
    pred.strain_p = bpa.strain_p(:);

    % Default strain field remains the physical/measured-comparison strain.
    pred.strain = pred.strain_p;

    pred.activeLength_effective = rest .* (1 - pred.strain_f);
    pred.activeLength_actual  = rest .* (1 - pred.strain_p);

    pred.delta_L = bpa.delta_L(:);
    pred.Lmt_p = bpa.Lmt_p(:);
    pred.Lmt0 = Lmt0(:);

    %% Design vaiables

    % Store design variables for optimizer output.
    pred.p1 = p1;
    pred.p8 = p8;
    pred.rest = rest;
    pred.tendon = tendon;
    pred.KMAX = KMAX;
    pred.kmax = kmax;
    pred.BPAcount = BPAcount;
    pred.bendMeasure = bendMeasure;
    pred.routeInfo = routeInfo;
    pred.tendonMax = routeInfo.tendonMax;
    pred.endcap_t1 = routeInfo.endcap_t1;

    % Backwards-compatible alias so existing scripts do not break immediately.
    pred.p2 = p8;
    
%% Path-length constraint
%
% In MonoPamDataExplicit_balanceX3:
%   Lmt_p = deformed_geometric_path_length - Xi0
%
% Therefore:
%   deformed_geometric_path_length = Lmt_p + Xi0
%
% Zero-strain model length:
%   rest + tendon + 2*fitting + Xi0
%
% Feasible if:
%   max(pathLength0) <= rest + tendon + 2*fitting + Xi0

% Deformed/stiffness-aware geometric path length diagnostic.
pred.pathLength = pred.Lmt_p + Xi0;

% Normal undeformed geometric path length used for rest-length constraint.
pred.pathLength0 = pred.Lmt0;

[pred.maxPathLength0, pred.idxMaxPath0] = max(pred.pathLength0);
pred.maxPathAngleD0 = ctx.phiD(pred.idxMaxPath0);

pred.restLmt = rest + tendon + 2*ctx.fitting + Xi0;

% <= 0 is feasible.
pred.cRestLength = pred.maxPathLength0 - pred.restLmt;

% Backward-compatible field names for existing scripts.
% These refer to the undeformed path-length constraint values.
pred.maxPathLength = pred.maxPathLength0;
pred.idxMaxPath = pred.idxMaxPath0;
pred.maxPathAngleD = pred.maxPathAngleD0;

% Deformed-path diagnostic max.
[pred.maxPathLengthDeformed, pred.idxMaxPathDeformed] = max(pred.pathLength);
pred.maxPathAngleDDeformed = ctx.phiD(pred.idxMaxPathDeformed);

% Backward-compatible field name for old print statements.
% This is NOT an extension-frame two-point distance anymore.
pred.extensionDistance = pred.maxPathLength;

%% Extra tendon diagnostic
% Use the minimum knee angle, because this is where the extensor path is longest.

[~, idxLongest] = min(ctx.phiD);

pred.tenEndSeg = routeInfo.tendonMax(idxLongest);
pred.tenEndSegAngleD = ctx.phiD(idxLongest);

end

%% =====================================================================
%% Local helper functions
%% =====================================================================

% function raw = distalRingRawLocation(phiD_i)
% % Location matrix from Knee_Extensor_10mm.
% %
% % Rows 1:4 are femur / parent frame.
% % Rows 5:8 are theta1 / tibia frame before conversion to ICR.
% %
% % Original Knee_Extensor_10mm points:
% % p1 = [0.040, 0.035, 0];                 %Origin
% % p2 = [0.0759, -0.27476, 0];             %BPA contacts mounting base
% % p3 = [0.05982, -0.37427, 0.000];        %femur channel contact, updated
% % p4 = [0.03955, -0.42183, 0.000];        %femoral condyle contact, updated
% % p5 = [0.05871, 0.025, 0];               %Tibia contact initial
% % p6 = [0.05871, 0.01228, 0];             %Tibia contact, updated
% % p7 = [0.054, -0.00612, 0];              %Tibia tendon contact, updated
% % p8 = [0.03604, -0.02844, 0];            %patellar ligament ring
% 
% p1 = [0.040,   0.035,    0.000];          %Origin
% p2 = [0.0759, -0.27476,  0.000];          %BPA contacts mounting base
% p3 = [0.05982,-0.37427,  0.000];          %femur channel contact, updated
% p4 = [0.03955,-0.42183,  0.000];          %femoral condyle contact, updated
% p5 = [0.05871, 0.025,    0.000];          %Tibia contact initial
% p6 = [0.05871, 0.01228,  0.000];          %Tibia contact, updated
% p7 = [0.054,  -0.00612,  0.000];          %Tibia tendon contact, updated
% p8 = [0.03604,-0.02844,  0.000];          %patellar ligament ring
% 
% %Set up angle limits (degrees)
% % AA = -68;
% BB = -69.8;   %add p4 when flexion reaches this value
% CC = -9.19;   %add p3 and p5 when flexion reaches this value
% DD = 10;
% 
%     if phiD_i > CC
%         raw = [p1;
%                p2;
%                p2;
%                p2;
%                p6;
%                p6;
%                p7;
%                p8];
% 
%     elseif phiD_i <= CC && phiD_i > BB
%         raw = [p1;
%                p2;
%                p3;
%                p3;
%                p5;
%                p6;
%                p7;
%                p8];
% 
%     elseif phiD_i <= BB
%         raw = [p1;
%                p2;
%                p3;
%                p4;
%                p5;
%                p6;
%                p7;
%                p8];
% 
%     else
%         raw = zeros(8,3);
%     end
% 
% end
% function Location = buildDistalRingLocation(p1, p8, phiD, T_ICR_t1)
% % Build the extensor Location matrix using the Knee_Extensor_10mm route.
% %
% % p1 replaces row 1.
% % p8 replaces row 8.
% % Rows 5:8 are converted from theta1 -> ICR with T_ICR_t1(:,:,i).
% 
% N = numel(phiD);
% Location = zeros(8, 3, N);
% 
%     for i = 1:N
%         raw = distalRingRawLocation(phiD(i));
% 
%         raw(1,:) = p1;
%         raw(8,:) = p8;
% 
%         for j = 5:8
%             raw(j,:) = RowVecTrans(T_ICR_t1(:,:,i), raw(j,:));
%         end
% 
%         Location(:,:,i) = raw;
%     end
% 
% end

function Lmt = muscleLengthNormal(Location, CrossPoint, T)
% Normal undeformed musculotendon length calculation, matching
% MonoPamDataExplicit SegmentLengths/MuscleLength behavior.

N = size(Location, 3);
M = size(Location, 1);
Lmt = zeros(N, 1);

for ii = 1:N
    for j = 1:M-1
        pointA = Location(j,:,ii);
        pointB = Location(j+1,:,ii);

        if j+1 == CrossPoint
            pointB = RowVecTrans(T(:,:,ii), pointB);
        end

        Lmt(ii) = Lmt(ii) + norm(pointA - pointB);
    end
end

end