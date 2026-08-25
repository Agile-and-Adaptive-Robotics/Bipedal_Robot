function J = objective_KneeExt20mm(x, ctx)

    %% New 20 mm routing constraints

    [cRoute, ~] = nonlconExt20mm(x, ctx);

    if any(~isfinite(cRoute))
        J = 1e12;
        return
    end


    %% Prediction

    pred = predictKneeExt20mm(x, ctx);

    if ~pred.ok
        J = 1e12;
        return
    end

    if any(~isfinite(pred.TorqueZ)) || any(~isfinite(pred.strain))
        J = 1e11;
        return
    end


    %% Human target

    humanAbs = interp1( ...
        ctx.humanAngleD, ...
        ctx.humanTorqueAbs, ...
        ctx.phiD, ...
        'pchip', ...
        'extrap');

    robotAbs = abs(pred.TorqueZ);


    %% Torque requirement

    torqueDeficit = max(0, humanAbs(:) - robotAbs(:));

    Jtorque = mean((torqueDeficit ./ ctx.torqueScale).^2);


    %% Small overshoot penalty

    torqueOvershoot = max(0, robotAbs(:) - humanAbs(:));

    Jovershoot = 1e-3 * ...
        mean((torqueOvershoot ./ ctx.torqueScale).^2);


    %% Rest-length constraint

    JrestLength = ...
        1e6 * max(0, pred.cRestLength).^2;


    %% Strain constraints

    maxStrainAllowed = ctx.maxRelStrain * pred.KMAX;

    JstrainHi = ...
        1e4 * max(0, ...
        max(pred.strain_f) - maxStrainAllowed).^2;

    JstrainLo = ...
        1e4 * max(0, ...
        ctx.minStrain - min(pred.strain_p)).^2;


    %% 20 mm routing geometry penalty
    %
    % Hard constraints are also supplied to the optimizers.
    % This penalty keeps direct calls to the objective well behaved.

    Jroute = 1e6 * sum(max(0,cRoute).^2);


    %% Mild design-change penalties

    x0 = ctx.x0(:);
    dx = x(:) - x0;

    Jgeom = 1e-2 * ...
        sum(dx(1:6).^2) / (0.060^2);

    Jlen = 1e-3 * ...
        ((x(7) - x0(7))/0.040).^2;


    %% Total objective

    J = ...
        100*Jtorque + ...
        Jovershoot + ...
        JrestLength + ...
        JstrainHi + ...
        JstrainLo + ...
        Jroute + ...
        Jgeom + ...
        Jlen;


    if ~isfinite(J)
        J = 1e12;
    end

end