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

    % Enforce the off-axis requirement separately for x and y.  At each
    % knee position, penalize only the amount by which the new component
    % magnitude exceeds the corresponding original no-wrap component.
    offAxisExcessX = max(0, ...
        abs(pred.TorqueX(:)) - abs(ctx.originalTorqueX(:)));
    offAxisExcessY = max(0, ...
        abs(pred.TorqueY(:)) - abs(ctx.originalTorqueY(:)));

    normalizedExcessX = offAxisExcessX ./ ctx.offAxisTorqueScaleX;
    normalizedExcessY = offAxisExcessY ./ ctx.offAxisTorqueScaleY;

    JoffAxis = ctx.offAxisPenaltyWeight * ( ...
        max(normalizedExcessX).^2 + mean(normalizedExcessX.^2) + ...
        max(normalizedExcessY).^2 + mean(normalizedExcessY.^2));

    % Resting length constraint:
    % rest + tendon + 2*fitting >= distance(p1, w2)
    % pred.cRestLength <= 0 is feasible.
    JrestLength = 1e6 * max(0, pred.cRestLength).^2;

    % Strain and no-load contraction feasibility.  strain_f includes Xi3
    % and drives force; strain_p excludes Xi3 and remains the stretch check.
    maxStrainAllowed = ctx.maxRelStrain * pred.KMAX;

    JstrainHi = 1e7 * max(0, ...
        max([pred.strain_f; pred.Contraction]) - maxStrainAllowed).^2;
    JstrainLo = 1e7 * max(0, ...
        -min(pred.Contraction)/pred.KMAX).^2;

    % The new contact must release during the extension-to-flexion sweep.
    % The constant term makes a no-release route unacceptable.
    if ~any(pred.routeInfo.active) || pred.routeInfo.releaseFound
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

    % Soft BPA-2 anchor (Ben, 2026-09-21): prefer keeping the second
    % BPA's DERIVED endpoints near their targets -- the Q-angle / aLDFA
    % build line of the loaded design.  Score term only, never a
    % constraint.  Weight 10 with a 1 cm norm: 1 cm drift on both ends
    % costs 20, while 1 N m of worst-angle torque shortfall costs ~4300,
    % so torque always wins when they conflict.
    if isfield(ctx, 'bpa2TargetP1') && ~isempty(ctx.bpa2TargetP1) && ...
            isfield(pred, 'p1B')
        anchorScale1 = norm(pred.p1B(:) - ctx.bpa2TargetP1(:))/ctx.bpa2AnchorNorm;
        anchorScale2 = norm(pred.p2B(:) - ctx.bpa2TargetP2(:))/ctx.bpa2AnchorNorm;
        Jbpa2Anchor = ctx.bpa2AnchorWeight*(anchorScale1^2 + anchorScale2^2);
    else
        Jbpa2Anchor = 0;
    end

    J = 1e5*Jworst + 1e3*Jtorque + 1e-2*Jshape + ...
        JoffAxis + JrestLength + JstrainHi + JstrainLo + ...
        JwrapRelease + Jgeom + Jlen + Jbpa2Anchor;

    if ~isfinite(J)
        J = 1e12;
    end
end
