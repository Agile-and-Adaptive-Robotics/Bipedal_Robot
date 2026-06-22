function J = objective_KneeExt20mm(x, ctx)

    pred = predictKneeFlexor20mm(x, ctx);

    if ~pred.ok
        J = 1e12;
        return
    end

    if any(~isfinite(pred.TorqueZ)) || any(~isfinite(pred.strain))
        J = 1e11;
        return
    end

    humanAbs = interp1( ...
        ctx.humanAngleD, ...
        ctx.humanTorqueAbs, ...
        ctx.phiD, ...
        'pchip', ...
        'extrap');

    robotAbs = abs(pred.TorqueZ);

    % Main requirement: robot >= human.
    torqueDeficit = max(0, humanAbs(:) - robotAbs(:));
    Jtorque = mean((torqueDeficit ./ ctx.torqueScale).^2);

    % Optional tiny overshoot penalty, only to avoid absurd geometry.
    torqueOvershoot = max(0, robotAbs(:) - humanAbs(:));
    Jovershoot = 1e-3 * mean((torqueOvershoot ./ ctx.torqueScale).^2);

    % Resting length constraint:
    % rest + tendon + 2*fitting >= distance(p1, w2)
    % pred.cRestLength <= 0 is feasible.
    JrestLength = 1e6 * max(0, pred.cRestLength).^2;

    % Strain feasibility.
    % KMAX is the measured free-contraction strain at 620 kPa.
    % maxRelStrain = 1 allows strain up to KMAX.
    maxStrainAllowed = ctx.maxRelStrain * pred.KMAX;
    
    JstrainHi = 1e4 * max(0, max(pred.strain) - maxStrainAllowed).^2;
    JstrainLo = 1e4 * max(0, ctx.minStrain - min(pred.strain)).^2;

    % Keep solution near practical geometry unless torque requires otherwise.
    x0 = ctx.x0(:);
    dx = x(:) - x0;

    Jgeom = 1e-2 * sum(dx(1:6).^2) / (0.060^2);
    Jlen  = 1e-3 * ((x(7) - x0(7))/0.040).^2;

    J = 100*Jtorque + Jovershoot + JrestLength + JstrainHi + JstrainLo + Jgeom + Jlen;

    if ~isfinite(J)
        J = 1e12;
    end
end