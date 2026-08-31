function J = objective_KneeFlexor20mm(x, ctx)

    pred = predictKneeFlexor20mm(x, ctx);

    if ~pred.ok
        J = 1e12;
        return
    end

    if any(~isfinite(pred.Torque(:))) || ...
            any(~isfinite(pred.strain_f)) || ...
            any(~isfinite(pred.strain_p)) || ...
            any(~isfinite(pred.Contraction))
        J = 1e11;
        return
    end

    humanAbs = interp1( ...
        ctx.humanAngleD, ...
        ctx.humanTorqueAbs, ...
        ctx.phiD, ...
        'pchip', ...
        'extrap');

    requiredTorque = (1 + ctx.requiredTorqueMargin)*humanAbs(:);
    robotAbs = abs(pred.TorqueZ(:));

    % Pointwise torque requirement with worst-position priority.
    torqueShortfall = max(0, requiredTorque - robotAbs);
    shortfallFraction = torqueShortfall ./ ctx.torqueScale;
    Jworst = max(shortfallFraction);
    Jtorque = mean(shortfallFraction.^2);

    % Very small shape term; meeting the requirement takes precedence.
    shapeError = (robotAbs - requiredTorque) ./ ctx.torqueScale;
    Jshape = mean(shapeError.^2);

    % Off-axis torque is allowed up to the pointwise value produced by x0.
    % Only excess resultant torque sqrt(Tx^2 + Ty^2) is penalized.
    offAxisExcess = max(0, ...
        pred.offAxisTorque(:) - ctx.offAxisTorqueX0(:));
    JoffAxis = ctx.offAxisPenaltyWeight * ...
        mean((offAxisExcess ./ ctx.offAxisTorqueScale).^2);

    % Resting length constraint:
    % rest + tendon + 2*fitting >= distance(p1, w2)
    % pred.cRestLength <= 0 is feasible.
    JrestLength = 1e6 * max(0, pred.cRestLength).^2;

    % Strain and no-load contraction feasibility.  strain_f includes Xi3
    % and drives force; strain_p excludes Xi3 and remains the stretch check.
    maxStrainAllowed = ctx.maxRelStrain * pred.KMAX;

    JstrainHi = 1e5 * max(0, ...
        max([pred.strain_f; pred.Contraction]) - maxStrainAllowed).^2;
    JstrainLo = 1e5 * max(0, ...
        ctx.minStrain - min(pred.strain_p)).^2;

    % The new contact must release during the extension-to-flexion sweep.
    % The constant term makes a no-release route unacceptable.
    if pred.routeInfo.releaseFound
        JwrapRelease = 0;
    else
        remainingTurn = max(0, -pred.routeInfo.finalSignedTurnD)/180;
        JwrapRelease = ctx.wrapReleasePenaltyWeight*(1 + remainingTurn).^2;
    end

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

    J = 1e5*Jworst + 1e3*Jtorque + 1e-2*Jshape + ...
        JoffAxis + JrestLength + JstrainHi + JstrainLo + ...
        JwrapRelease + Jgeom + Jlen;

    if ~isfinite(J)
        J = 1e12;
    end
end
