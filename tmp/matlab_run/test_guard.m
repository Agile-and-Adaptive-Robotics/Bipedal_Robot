% Verify app fixes + reader script + torque-guard math.
here = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Code/Matlab/HX711 v3.0/HX711 v3.0';
td  = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Testing_Data/2026_06_Festo';

% 1) lint both files
r1 = checkcode(fullfile(here,'HX711_BPA.m'));
r2 = checkcode(fullfile(td,'readserialnumbers2.m'));
hard1 = r1(contains(lower({r1.message}),'parse error|unterminated|might be missing'));
hard2 = r2(contains(lower({r2.message}),'parse error|unterminated|might be missing'));
assert(isempty(hard1) && isempty(hard2), 'lint hard errors');
fprintf('lint: app %d msgs, reader %d msgs (no hard errors)\n', numel(r1), numel(r2));

% 2) app regression
cd(here);
test_HX711_BPA_offline;

% 3) torque-guard math (same as computeTorqueZ in readserialnumbers2)
root = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(root,'Code','Matlab')));
addpath(fullfile(root,'Code','Matlab','Mesh_Optimization'));
addpath(fullfile(root,'Testing_Data','2022_02_Festo'), '-end');
ctx = buildKneeExtContext20mm();
dLC = 292.9/1000; angLC = -90.83;
pRF = [dLC*cosd(angLC), dLC*sind(angLC), 0];
phiV = ctx.phi(:);
txV = squeeze(ctx.T_t1_ICR(1,4,:)) - pRF(1); txV = txV(:);
tyV = squeeze(ctx.T_t1_ICR(2,4,:)) - pRF(2); tyV = tyV(:);

Kdeg = -30; Ldeg = 30; F = 50;   % typical row
Kr = deg2rad(Kdeg); Lr = deg2rad(Ldeg);
tx = interp1(phiV, txV, Kr, 'pchip');
ty = interp1(phiV, tyV, Kr, 'pchip');
Trk = RpToTrans(eye(3), [tx; ty; 0]);
Fr = -[0; 0; 0; F*cos(pi - Lr); F*sin(pi - Lr); 0];
Fk = Adjoint(Trk)'*Fr;
tz = Fk(3);
tgt = interp1(ctx.humanAngleD(:), ctx.humanTorque(:), Kdeg, 'pchip');
planar = F*cosd(Ldeg - 2.83) * (292.9)/1000;  % Excel row-16 style approx
fprintf('guard math: torque %.3f N*m (planar approx %.3f), human target %.3f N*m\n', tz, planar, tgt);
assert(isfinite(tz) && tz > 0, 'torque not finite/positive');
assert(abs(tz - planar)/planar < 0.25, 'adjoint far from planar sanity');
assert(isfinite(tgt) && tgt > 0, 'target not finite');
fprintf('GUARD MATH PASS\n');
