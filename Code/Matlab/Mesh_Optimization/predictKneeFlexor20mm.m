function pred = predictKneeFlexor20mm(x, ctx)
%PREDICTKNEEFLEXOR20MM Evaluate one flexor design over the full angle grid.
% ctx supplies values fixed for the run; x supplies the eight optimizer
% variables. pred collects the route, class object, torque, lengths,
% strains, radius iteration diagnostics, and constraint-related values so
% the objective, nonlinear constraints, reports, and plots stay consistent.

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
    % Saved contexts from before radius offsets remain usable.
    if ~isfield(geoUsed,'bpaRbOffset'), geoUsed.bpaRbOffset = 0; end
    if ~isfield(geoUsed,'bpaRsOffset'), geoUsed.bpaRsOffset = 0; end
    % radiusMode comes from ctx.bpaRadiusMode in the context builder.
    % Use "scalar" there for fixed geo.bpaRb/geo.bpaRs. Do not replace the
    % bpaR call with zero: bpaR(0,...) means the physical radius at zero
    % contraction, whereas scalar mode bypasses bpaR for the entire motion.
    radiusMode = string(geoUsed.bpaRadiusMode);

    if radiusMode == "bpaR"
        % Start from a physical zero-contraction radius, not either scalar
        % guess. The update below then uses the modeled physical strain_p.
        radius0 = bpaR(zeros(ctx.N,1), ctx.Dia, KMAX);
        geoUsed.bpaRb = radius0 + geoUsed.bpaRbOffset;
        geoUsed.bpaRs = radius0 + geoUsed.bpaRsOffset;
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

            % routeInfo.pWrapT1 is the moving t1-frame wrap-point array.
            % When active, its transformed value becomes Location(2,:,i).
            [Location1, bendMeasure1, routeInfo1] = ...
                buildKneeFlexorRoute20mm(p1, p2, tendon, ctxUsed);

            if BPAcount == 2
                % Ben, 2026-09-21: BPA 2 keeps its mirrored distal
                % attachment (pEnd{2}z = -pEnd{1}z) but its origin is no
                % longer the xy-plane mirror of BPA 1.  Instead the two
                % BPAs share the same p1-to-pEnd z offset:
                %   p1{2}z = p1{1}z - (pEnd{1}z - pEnd{2}z)
                % (see flexorBpa2Endpoints20mm), so both origins sit off
                % the same side of the knee like the proposed CAD routing.
                % The second route is SOLVED by buildKneeFlexorRoute20mm,
                % not mirrored: its wrap point and tibia/femur collision
                % state are resolved for this asymmetric path.
                [p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);
                [Location2, bendMeasure2, routeInfo2] = ...
                    buildKneeFlexorRoute20mm(p1B, p2B, tendon, ctxUsed);

                Location = {Location1; Location2};
                bendMeasure = {bendMeasure1; bendMeasure2};
                bpa = MonoPam_mult( ...
                    ctx.Name, Location, ctx.CrossPoint, ctx.Dia, ctx.T_Pam, ...
                    rest, kmax, tendon, ctx.fitting, ctx.targetPressure, ...
                    Xi0, Xi1, Xi2, Xi3, ctx.wraps, ctx.phiD, BPAcount, ...
                    bendMeasure);
            elseif BPAcount == 1
                Location = Location1;
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
                bendMeasure1);
            else
                error('predictKneeFlexor20mm:BPAcount', ...
                    'BPAcount must be 1 or 2 for the current flexor model.')
            end

            if radiusMode == "scalar"
                break
            end

            % strain_p excludes Xi3 and is the physical BPA contraction.
            % bpaR therefore returns the physical outer radius at each pose.
            if BPAcount == 2
                % The two routes now wrap differently (asymmetric
                % routing).  Average the two strain arrays for the one
                % shared radius used to rebuild both routes;
                % forceMismatch remains the diagnostic of the mechanical
                % asymmetry between the two paths.
                strainRadius = mean([bpa.strain_p{1},bpa.strain_p{2}],2);
            else
                strainRadius = bpa.strain_p(:);
            end
            radiusModel = bpaR(strainRadius, ctx.Dia, KMAX);

            if numel(radiusModel) ~= ctx.N || ...
                    any(~isfinite(radiusModel)) || any(radiusModel <= 0)
                error('predictKneeFlexor20mm:InvalidBpaR', ...
                    'bpaR must return %d finite positive radii.', ctx.N)
            end

            radiusModel = radiusModel(:);
            radiusChange = max([ ...
                abs(radiusModel + geoUsed.bpaRbOffset - geoUsed.bpaRb(:)); ...
                abs(radiusModel + geoUsed.bpaRsOffset - geoUsed.bpaRs(:))]);

            if radiusChange <= geoUsed.bpaRTolerance
                radiusConverged = true;
                break
            end

            if radiusIteration == maxRadiusIterations
                break
            end

            % Keep your separate wrap/collision offsets in every iteration.
            geoUsed.bpaRb = radiusModel + geoUsed.bpaRbOffset;
            geoUsed.bpaRs = radiusModel + geoUsed.bpaRsOffset;
        end
    catch ME
        pred.ok = false;
        pred.failReason = string(ME.message);
        return
    end

    pred.bpa = bpa;
    pred.Location = Location1;  % BPA 1 route; nonlcon checks both routes
    pred.LocationAll = Location;
    pred.bendMeasure = bendMeasure1;         % BPA 1 bend history
    if BPAcount == 2
        pred.Location2 = Location2;          % BPA 2 (same-side origin)
        pred.bendMeasure2 = bendMeasure2;
        pred.routeInfo2 = routeInfo2;
        pred.p1B = p1B;                      % femur frame
        pred.p2B = p2B;                      % t1 frame
    end
    pred.routeInfo = routeInfo1;
    pred.geo = geoUsed;
    pred.bpaRadiusMode = radiusMode;
    pred.bpaRadius = routeInfo1.bpaRs;
    pred.bpaRadiusIteration = radiusIteration;
    pred.bpaRadiusConverged = radiusConverged;
    pred.bpaRadiusChange = radiusChange;
    pred.Torque   = bpa.Torque_p;
    pred.TorqueX  = bpa.Torque_p(:,1);
    pred.TorqueY  = bpa.Torque_p(:,2);
    pred.TorqueZ  = bpa.Torque_p(:,3);
    pred.offAxisTorque = hypot(bpa.Torque_p(:,1), bpa.Torque_p(:,2));
    if BPAcount == 2
        pred.ContractionAll = bpa.Contraction;
        pred.strain_fAll = bpa.strain_f;
        pred.strain_pAll = bpa.strain_p;
        pred.Contraction = bpa.Contraction{1}(:);
        pred.strain_f = bpa.strain_f{1}(:);
        pred.strain_p = bpa.strain_p{1}(:);
        pred.forceMismatch = bpa.forceMismatch;
        pred.momentArmVectorAll = bpa.mA_p;
        pred.momentArmVector = bpa.mA_p{1};
        pred.pathLength0 = bpa.MuscleLength{1}(:);
        pred.pathLength = bpa.Lmt_p{1}(:) + Xi0;
        pred.delta_L = bpa.delta_L{1}(:);
        pred.gama = bpa.gama{1}(:);
    else
        pred.Contraction = bpa.Contraction(:);
        pred.strain_f = bpa.strain_f(:);
        pred.strain_p = bpa.strain_p(:);
        pred.forceMismatch = zeros(ctx.N,1);
        pred.momentArmVector = bpa.mA_p;
        pred.pathLength0 = bpa.MuscleLength(:);
        pred.pathLength = bpa.Lmt_p(:) + Xi0;
        pred.delta_L = bpa.delta_L(:);
        pred.gama = bpa.gama(:);
    end
    pred.strain = pred.strain_p;  % compatibility with existing checks
    pred.relativeContraction = pred.Contraction ./ KMAX;
    pred.relativeStrainF = pred.strain_f ./ KMAX;
    pred.relativeStrainP = pred.strain_p ./ KMAX;
    pred.relativeStrain = pred.relativeStrainP;
    pred.activeLength = rest .* (1 - pred.strain_p);
    pred.momentArm = hypot(pred.momentArmVector(:,1),pred.momentArmVector(:,2));

    % Store design variables for optimizer output.
    pred.p1 = p1;
    pred.p2 = p2;
    pred.pEnd = p2;
    pred.pWrap = routeInfo1.pWrapT1(ctx.idxExtension,:);
    pred.rest = rest;
    pred.tendon = tendon;
    pred.KMAX = KMAX;
    pred.kmax = kmax;
    pred.BPAcount = BPAcount;
    
    % Extension-frame geometry check.
    idx = ctx.idxExtension;
    
    % v2 is the final insertion point in the knee/ICR frame.
    v2 = Location1(3,:,idx);
    
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
