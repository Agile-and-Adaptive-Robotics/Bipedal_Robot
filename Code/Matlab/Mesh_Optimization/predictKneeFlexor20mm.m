function pred = predictKneeFlexor20mm(x, ctx)

    p1       = x(1:3);       % parent-frame / femur-side attachment
    p2       = x(4:6);       % theta1-frame end/insertion design variable
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

    geoUsed = ctx.geo;
    radiusMode = string(geoUsed.bpaRadiusMode);

    if radiusMode == "bpaR"
        % Start from a physical zero-contraction radius, not either scalar
        % guess. The update below then uses the modeled physical strain_p.
        radius0 = bpaR(zeros(ctx.N,1), ctx.Dia, KMAX);
        geoUsed.bpaRb = radius0;
        geoUsed.bpaRs = radius0;
        maxRadiusIterations = geoUsed.bpaRMaxIterations;
    else
        maxRadiusIterations = 1;
    end

    radiusConverged = radiusMode == "scalar";
    radiusChange = 0;

    try
        for radiusIteration = 1:maxRadiusIterations
            ctxUsed = ctx;
            ctxUsed.geo = geoUsed;

            [Location, bendMeasure, routeInfo] = ...
                buildKneeFlexorRoute20mm(p1, p2, tendon, ctxUsed);

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

            if radiusMode == "scalar"
                break
            end

            % strain_p excludes Xi3 and is the physical BPA contraction.
            % bpaR therefore returns the physical outer radius at each pose.
            radiusModel = bpaR(bpa.strain_p(:), ctx.Dia, KMAX);

            if numel(radiusModel) ~= ctx.N || ...
                    any(~isfinite(radiusModel)) || any(radiusModel <= 0)
                error('predictKneeFlexor20mm:InvalidBpaR', ...
                    'bpaR must return %d finite positive radii.', ctx.N)
            end

            radiusModel = radiusModel(:);
            radiusChange = max([ ...
                abs(radiusModel - geoUsed.bpaRb(:)); ...
                abs(radiusModel - geoUsed.bpaRs(:))]);

            if radiusChange <= geoUsed.bpaRTolerance
                radiusConverged = true;
                break
            end

            if radiusIteration == maxRadiusIterations
                break
            end

            % With a physical radius model there is one radius per pose.
            % The bpaRb/bpaRs names remain separate only because they serve
            % different routing and collision roles in scalar mode.
            geoUsed.bpaRb = radiusModel;
            geoUsed.bpaRs = radiusModel;
        end
    catch ME
        pred.ok = false;
        pred.failReason = string(ME.message);
        return
    end

    pred.bpa = bpa;
    pred.Location = Location;
    pred.bendMeasure = bendMeasure;
    pred.routeInfo = routeInfo;
    pred.geo = geoUsed;
    pred.bpaRadiusMode = radiusMode;
    pred.bpaRadius = routeInfo.bpaRs;
    pred.bpaRadiusIteration = radiusIteration;
    pred.bpaRadiusConverged = radiusConverged;
    pred.bpaRadiusChange = radiusChange;
    pred.Torque   = bpa.Torque_p;
    pred.TorqueX  = bpa.Torque_p(:,1);
    pred.TorqueY  = bpa.Torque_p(:,2);
    pred.TorqueZ  = bpa.Torque_p(:,3);
    pred.offAxisTorque = hypot(bpa.Torque_p(:,1), bpa.Torque_p(:,2));
    pred.Contraction = bpa.Contraction(:);
    pred.strain_f = bpa.strain_f(:);
    pred.strain_p = bpa.strain_p(:);
    pred.strain = pred.strain_p;  % compatibility with existing checks
    pred.relativeContraction = pred.Contraction ./ KMAX;
    pred.relativeStrainF = pred.strain_f ./ KMAX;
    pred.relativeStrainP = pred.strain_p ./ KMAX;
    pred.relativeStrain = pred.relativeStrainP;
    pred.activeLength = rest .* (1 - pred.strain_p);
    pred.pathLength0 = bpa.MuscleLength(:);
    pred.pathLength = bpa.Lmt_p(:) + Xi0;
    pred.delta_L = bpa.delta_L(:);
    pred.gama = bpa.gama(:);
    pred.momentArmVector = bpa.mA_p;
    pred.momentArm = hypot(bpa.mA_p(:,1), bpa.mA_p(:,2));

    % Store design variables for optimizer output.
    pred.p1 = p1;
    pred.p2 = p2;
    pred.pEnd = p2;
    pred.pWrap = routeInfo.pWrapT1(ctx.idxExtension,:);
    pred.rest = rest;
    pred.tendon = tendon;
    pred.KMAX = KMAX;
    pred.kmax = kmax;
    pred.BPAcount = BPAcount;
    
    % Extension-frame geometry check.
    idx = ctx.idxExtension;
    
    % v2 is the final insertion point in the knee/ICR frame.
    v2 = Location(3,:,idx);
    
    % w2 is that same point transformed into the femur frame.
    w2 = RowVecTrans(ctx.T_Pam(:,:,idx), v2);
    
    pred.v2 = v2;
    pred.w2 = w2;
    
    pred.extensionDistance = pred.pathLength0(idx);
    
    % Physical no-load musculotendon length available between attachment points.
    % Xi0 is not included here because this is a physical packaging constraint.
    % pred.restLmt = rest + tendon + 2*ctx.fitting;

    % Diagnostic only: this is the modeled zero-strain musculotendon length.
    % Positive Xi0 makes the model behave as if the required Lmt is longer.
    pred.restLmt = rest + tendon + 2*ctx.fitting + Xi0;

    % Constraint value <= 0 is feasible:
    % distance(p1, v2 at extension) <= rest + tendon + 2*fitting
    pred.cRestLength = pred.extensionDistance - pred.restLmt;

end
