function pred = predictKneeExt20mm(x, ctx)
% predictKneeExt20mm
%
% Builds the 20 mm extensor distal-ring path, evaluates the X3 stiffness
% model, and returns torque/strain/path-length diagnostics.
%
% x = [p1(1:3), p2(1:3), rest, tendon]
%
% p1 is row 1 of the distal-ring matrix in femur/parent frame.
% p2 is row 6 of the distal-ring matrix in theta1/tibia frame.
%
% Rows 4:6 are converted theta1 -> ICR with T_ICR_t1(:,:,i), then
% MonoPamDataExplicit_balanceX3 uses T_Pam(:,:,i) for the crossing segment.

    p1       = x(1:3);       % parent-frame / femur-side attachment
    p2       = x(4:6);       % theta1-frame insertion design variable
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
    Location = buildDistalRingLocation(p1, p2, ctx.phiD, ctx.T_ICR_t1);

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
                BPAcount);
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
    pred.p2 = p2;
    pred.rest = rest;
    pred.tendon = tendon;
    pred.KMAX = KMAX;
    pred.kmax = kmax;
    pred.BPAcount = BPAcount;
    
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
    %   max(pathLength) <= rest + tendon + 2*fitting + Xi0

    pred.pathLength = pred.Lmt_p + Xi0;
    pred.pathLength0 = pred.Lmt0;

    [pred.maxPathLength, pred.idxMaxPath] = max(pred.pathLength);
    pred.maxPathAngleD = ctx.phiD(pred.idxMaxPath);

    [pred.maxPathLength0, pred.idxMaxPath0] = max(pred.pathLength0);
    pred.maxPathAngleD0 = ctx.phiD(pred.idxMaxPath0);

    pred.restLmt = rest + tendon + 2*ctx.fitting + Xi0;

    % <= 0 is feasible.
    pred.cRestLength = pred.maxPathLength - pred.restLmt;

    % Backward-compatible field name for old print statements.
    % This is NOT an extension-frame two-point distance anymore.
    pred.extensionDistance = pred.maxPathLength;

    %% Extra tendon diagnostic
    rawHome = distalRingRawLocation(ctx.phiD(ctx.pos));
    rawHome(1,:) = p1;
    rawHome(6,:) = p2;
    pred.tendonRow56Home = norm(rawHome(5,:) - rawHome(6,:));

end

%% =====================================================================
%% Local helper functions
%% =====================================================================

function raw = distalRingRawLocation(phiD_i)
% Distal ring Location matrix from Knee_Extensor_40mm.
%
% Rows 1:3 are femur / parent frame.
% Rows 4:6 are theta1 / tibia frame before conversion to ICR.

if phiD_i < -80
    raw = [ ...
        0.040,   0.035,   0.000;
        0.099,  -0.275,   0.000;
        0.05546,-0.44305, 0.000;
        0.06594,-0.010,   0.000;
        0.06594,-0.03716, 0.000;
        0.01969,-0.11115, 0.000];

elseif phiD_i >= -80 && phiD_i < -40
    raw = [ ...
        0.040,   0.035,   0.000;
        0.099,  -0.240,   0.000;
        0.08317,-0.385,   0.000;
        0.07094,-0.0129,  0.000;
        0.07094,-0.03716, 0.000;
        0.01969,-0.11115, 0.000];

else
    raw = [ ...
        0.040,   0.035,   0.000;
        0.099,  -0.219,   0.000;
        0.099,  -0.30252, 0.000;
        0.07094,-0.0129,  0.000;
        0.07094,-0.03716, 0.000;
        0.01969,-0.11115, 0.000];
end

end

function Location = buildDistalRingLocation(p1, p2, phiD, T_ICR_t1)
% Build the 6-row extensor distal-ring Location matrix.
%
% p1 replaces row 1.
% p2 replaces row 6.
% Rows 4:6 are converted from theta1 -> ICR with T_ICR_t1(:,:,i).

N = numel(phiD);
Location = zeros(6, 3, N);

for i = 1:N
    raw = distalRingRawLocation(phiD(i));

    raw(1,:) = p1;
    raw(6,:) = p2;

    for j = 4:6
        raw(j,:) = RowVecTrans(T_ICR_t1(:,:,i), raw(j,:));
    end

    Location(:,:,i) = raw;
end

end

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