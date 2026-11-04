function res = fr_validateKin(C, kin)
%FR_VALIDATEKIN Validate the parsed kinematics against the geometry
% extraction's ground-truth bodyTransformsInGround (geometry_20261103,
% generated 2026-11-03 with OpenSim 4.4.1 and verified 0.0 mm per-point
% there). Comparison is exact transform agreement (R and p), NOT
% bounding-box agreement (master directive item 1 forbids bbox "proof").
%
% The extraction's realized_deg stores rotational coordinates in degrees;
% pelvis_ty (translational) is stored as meters*pi/180 (0.95 m -> 54.43),
% so dividing every stored value by rad2deg(1) recovers radians for
% rotations and meters for translations in one consistent rule.
%
% res.maxErrMm: max over poses and bodies of |p_mine - p_opensim|;
% res.maxRotErrDeg: max rotation-matrix angle error.

ex = jsondecode(fileread(fullfile(C.geom, 'extraction_robotbody.json')));
poses = fieldnames(ex.bodyTransformsInGround);
res.maxErrMm = 0;
res.maxRotErrDeg = 0;
res.perPose = struct();
for ip = 1:numel(poses)
    pn = poses{ip};
    rd = ex.poses.(pn).realized_deg;
    q = struct();
    fn = fieldnames(rd);
    for i = 1:numel(fn)
        q.(fn{i}) = rd.(fn{i}) / rad2deg(1);
    end
    worst = 0; worstRot = 0;
    btg = ex.bodyTransformsInGround.(pn);
    bn = fieldnames(btg);
    for ib = 1:numel(bn)
        body = bn{ib};
        T = fr_bodyT(kin, q, body);
        dp = norm(T(1:3, 4) - btg.(body).p(:)) * 1000;
        dR = rotErrDeg(T(1:3, 1:3), btg.(body).R);
        worst = max(worst, dp);
        worstRot = max(worstRot, dR);
    end
    res.perPose.(pn).maxErrMm = worst;
    res.perPose.(pn).maxRotErrDeg = worstRot;
    res.maxErrMm = max(res.maxErrMm, worst);
    res.maxRotErrDeg = max(res.maxRotErrDeg, worstRot);
end
end

function d = rotErrDeg(R1, R2)
dR = R1 * R2';
ang = acos(min(1, max(-1, (trace(dR) - 1) / 2)));
d = rad2deg(ang);
end
