function pred = predictOriginalKneeFlexor20mm(ctx)
%PREDICTORIGINALKNEEFLEXOR20MM Original two-point route with zero X3 bend.

p1 = ctx.originalP1(:).';
pEnd = ctx.originalPEnd(:).';

Location = zeros(2,3,ctx.N);
for i = 1:ctx.N
    pEndICR = RowVecTrans(ctx.T_ICR_t1(:,:,i), pEnd);
    Location(:,:,i) = [p1; pEndICR];
end

rest = ctx.originalRest;
kmax = rest*(1 - ctx.KMAX);
tendon = ctx.originalTendon;

pred.ok = true;
pred.failReason = "";
bendMeasure = zeros(ctx.N,1);

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
        ctx.Xi0, ...
        ctx.Xi1, ...
        ctx.Xi2, ...
        ctx.Xi3, ...
        ctx.originalWraps, ...
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
pred.Torque = bpa.Torque_p;
pred.TorqueX = bpa.Torque_p(:,1);
pred.TorqueY = bpa.Torque_p(:,2);
pred.TorqueZ = bpa.Torque_p(:,3);
pred.Contraction = bpa.Contraction(:);
pred.strain_p = bpa.strain_p(:);
pred.activeLength = rest.*(1 - pred.strain_p);
pred.pathLength0 = bpa.MuscleLength(:);
pred.momentArmVector = bpa.mA_p;
pred.momentArm = hypot(bpa.mA_p(:,1), bpa.mA_p(:,2));
pred.p1 = p1;
pred.pEnd = pEnd;
pred.rest = rest;
pred.tendon = tendon;
pred.KMAX = ctx.KMAX;
pred.kmax = kmax;

end
