function J = objective_KneeFlexor20mm(x, ctx)

    pred = predictKneeFlexor20mm(x, ctx);

    if ~pred.ok
        J = 1e12;
        return
    end

    idxOperating = ctx.phiD >= -120 & ctx.phiD <= 10;

    if any(~isfinite(pred.TorqueZ(idxOperating))) || any(~isfinite(pred.strain(idxOperating)))
        J = 1e11;
        return
    end

    humanAbs = interp1(ctx.humanAngleD, ctx.humanTorqueAbs, ctx.phiD, 'pchip', 'extrap');

    humanAbs = humanAbs(idxOperating);
    robotAbs = abs(pred.TorqueZ(idxOperating));

    requiredAbs = 1.02*humanAbs(:);
    humanScale = max(humanAbs(:), 1);

    shortfallFraction = max(0, (requiredAbs - robotAbs(:))./humanScale);

    Jworst = max(shortfallFraction);
    Jmean = mean(shortfallFraction.^2);

    shapeError = (robotAbs(:) - requiredAbs)./humanScale;
    Jshape = mean(shapeError.^2);

    % Resting length constraint:
    % rest + tendon + 2*fitting >= distance(p1, w2)
    % pred.cRestLength <= 0 is feasible.
    JrestLength = 1e6 * max(0, pred.cRestLength).^2;

    % Strain feasibility.
    % KMAX is the measured free-contraction strain at 620 kPa.
    % maxRelStrain = 1 allows strain up to KMAX.
    maxStrainAllowed = ctx.maxRelStrain * pred.KMAX;
    
    operatingStrain = pred.strain(idxOperating);

    JstrainHi = 1e4*max(0, max(operatingStrain) - maxStrainAllowed).^2;
    JstrainLo = 1e4*max(0, ctx.minStrain - min(operatingStrain)).^2;

    % Keep solution near practical geometry unless torque requires otherwise.
    x0 = ctx.x0(:);
    dx = x(:) - x0;

    geomScale = [ ...
    0.030;   % p1 x
    0.300;   % p1 y -- large movement allowed
    0.030;   % p1 z
    0.030;   % p2 x
    0.030;   % p2 y
    0.030];  % p2 z

    Jgeom = 1e-2 * sum((dx(1:6)./geomScale).^2);
    Jlen  = 1e-3 * ((x(7) - x0(7))/0.040).^2;

    J = 1e5*Jworst + 1e3*Jmean + 1e-2*Jshape + JrestLength + JstrainHi + JstrainLo + Jgeom + Jlen;

    if ~isfinite(J)
        J = 1e12;
    end
end