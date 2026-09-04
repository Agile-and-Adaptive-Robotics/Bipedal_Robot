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

    KMAX = ctx.KMAX;
    kmax = rest*(1-KMAX);      % measured free-contracted BPA length at 620 kPa

    pred.ok = true;
    pred.failReason = "";

    if rest <= 0 || kmax <= 0 || kmax >= rest || tendon < 0
        pred.ok = false;
        pred.failReason = "Invalid length parameters";
        return
    end

    [Location, bendMeasure, routeInfo] = ...
        buildKneeFlexorRoute20mm(p1, p2, tendon, ctx);

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
            ctx.BPAcount, ...
            bendMeasure);
    catch ME
        pred.ok = false;
        pred.failReason = string(ME.message);
        return
    end

    pred.bpa = bpa;
    pred.Location = Location;
    pred.bendMeasure = bendMeasure;
    pred.routeInfo = routeInfo;
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
