function pred = predictKneeExt20mm(x, ctx)
% predictKneeExt20mm
%
% Builds the 20 mm extensor distal-ring path, evaluates the X3 stiffness
% model, and returns torque/strain/path-length diagnostics.
%
% Optimization vector:
% x = [p1(1:3), pEnd(1:3), rest, tendon]
%
% Route rows:
%   p1:p5 femur frame
%   p6:p9 t1 frame before t1 -> ICR conversion in the route builder
%   p9 = pEnd
%   CrossPoint = 6

    p1       = x(1:3);       % femur-side attachment
    pEnd     = x(4:6);       % distal t1 insertion design variable, row p9
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
    try
        [Location, bendMeasure, routeInfo] = buildDistalRingLocation20mm(p1,pEnd,tendon,ctx);
    catch ME
        pred.ok = false;
        pred.failReason = string(ME.message);
        return
    end

    %% Normal, undeformed path length calculation
    Lmt0 = muscleLengthNormal(Location, ctx.CrossPoint, ctx.T_Pam);

    %% Evaluate X3/stiffness-aware class
    try
        bpa = MonoPamDataExplicit_balanceX3(ctx.Name,Location,ctx.CrossPoint,ctx.Dia,ctx.T_Pam,rest,kmax,tendon,ctx.fitting,ctx.targetPressure,Xi0,Xi1,Xi2,Xi3,ctx.wraps,ctx.phiD,BPAcount,bendMeasure);
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

    pred.momentArmVector = bpa.mA_p;
    pred.momentArm = hypot(bpa.mA_p(:,1), bpa.mA_p(:,2));

    pred.delta_L = bpa.delta_L(:);
    pred.Lmt_p = bpa.Lmt_p(:);
    pred.Lmt0 = Lmt0(:);

    %% Design variables
    pred.p1 = p1;

    pred.pEnd = pEnd;
    pred.p9 = pEnd;

    pred.rest = rest;
    pred.tendon = tendon;
    pred.KMAX = KMAX;
    pred.kmax = kmax;
    pred.BPAcount = BPAcount;

    pred.bendMeasure = bendMeasure;
    pred.routeInfo = routeInfo;

    pred.tendonMax = routeInfo.tendonMax;
    pred.endcap_t1 = routeInfo.endcap_t1;

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
    pred.maxPathLength = pred.maxPathLength0;
    pred.idxMaxPath = pred.idxMaxPath0;
    pred.maxPathAngleD = pred.maxPathAngleD0;

    % Deformed-path diagnostic max.
    [pred.maxPathLengthDeformed, pred.idxMaxPathDeformed] = max(pred.pathLength);
    pred.maxPathAngleDDeformed = ctx.phiD(pred.idxMaxPathDeformed);

    % Backward-compatible field name for old print statements.
    pred.extensionDistance = pred.maxPathLength;

    %% Extra tendon diagnostic
    [~, idxLongest] = min(ctx.phiD);

    pred.tenEndSeg = routeInfo.tendonMax(idxLongest);
    pred.tenEndSegAngleD = ctx.phiD(idxLongest);

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
