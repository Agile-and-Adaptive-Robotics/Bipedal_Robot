function pred = predictKneeExt20mm(x, ctx)

    p1       = x(1:3);       % parent-frame / femur-side attachment
    p2       = x(4:6);       % theta1-frame insertion design variable
    rest     = x(7);         % active BPA rest length
    tendon   = x(8);         % physical tendon length
    
    % Fixed identified stiffness parameters
    Xi0 = ctx.Xi0;
    Xi1 = ctx.Xi1;
    Xi2 = ctx.Xi2;

    KMAX = ctx.KMAX;
    kmax = rest*(1-KMAX);      % measured free-contracted BPA length at 620 kPa

    pred.ok = true;
    pred.failReason = "";

    if rest <= 0 || kmax <= 0 || kmax >= rest || tendon < 0
        pred.ok = false;
        pred.failReason = "Invalid length parameters";
        return
    end

    N = ctx.N;
    Location = zeros(2, 3, N);

    for i = 1:N
        % p2 is defined in the theta1 frame.
        % Convert it into the knee/ICR frame before storing it in Location.
        v2 = RowVecTrans(ctx.T_ICR_t1(:, :, i), p2);
        Location(:, :, i) = [p1; v2];
    end

    try
        bpa = MonoPamDataExplicit_balance( ...
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
            ctx.wraps);
    catch ME
        pred.ok = false;
        pred.failReason = string(ME.message);
        return
    end

    pred.Location = Location;
    pred.Torque   = bpa.Torque_p;
    pred.TorqueZ  = bpa.Torque_p(:,3);
    pred.strain   = bpa.strain_p(:);
    pred.activeLength = rest .* (1 - pred.strain);

    % Store design variables for optimizer output.
    pred.p1 = p1;
    pred.p2 = p2;
    pred.rest = rest;
    pred.tendon = tendon;
    pred.KMAX = KMAX;
    pred.kmax = kmax;
    
    % Extension-frame geometry check.
    idx = ctx.idxExtension;
    
    % v2 is the Location row-2 point in the knee/ICR frame.
    v2 = Location(2,:,idx);
    
    % w2 is that same point transformed into the femur frame.
    w2 = RowVecTrans(ctx.T_Pam(:,:,idx), v2);
    
    pred.v2 = v2;
    pred.w2 = w2;
    
    pred.extensionDistance = norm(p1 - w2);
    
    % Physical no-load musculotendon length available between attachment points.
    % Xi0 is not included here because this is a physical packaging constraint.
    pred.restLmt = rest + tendon + 2*ctx.fitting;

    % Diagnostic only: this is the modeled zero-strain musculotendon length.
    % Positive Xi0 makes the model behave as if the required Lmt is longer.
    % pred.restLmt = rest + tendon + 2*ctx.fitting + Xi0;

    % Constraint value <= 0 is feasible:
    % distance(p1, v2 at extension) <= rest + tendon + 2*fitting
    pred.cRestLength = pred.extensionDistance - pred.restLmt;

end