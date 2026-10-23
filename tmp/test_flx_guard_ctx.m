function test_flx_guard_ctx()
%TEST_FLX_GUARD_CTX Headless check of the readserialnumbers2 flexor
% torque-guard wiring (context swap + LC convention + verdict signs).
% Does NOT call readserialnumbers2 and opens no serial port -- it only
% replicates that function's model-load block and torque math.

functionFolder = fullfile('D:', 'GitHub', 'Bipedal_Robot', ...
    'Testing_Data', '2026_06_Festo');
root = functionFolder;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break;
    end
    if strcmp(parent, root)
        error('Could not locate the Bipedal_Robot repo root.');
    end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

evalc('ctx = buildKneeFlexorContext20mm();');

% Torque-guard geometry as edited into readserialnumbers2 (42 cm defaults)
dLC   = 320/1000;
angLC = -82.97;
pRF = [dLC*cosd(angLC), dLC*sind(angLC), 0];
phiV = ctx.phi(:);
txV = squeeze(ctx.T_t1_ICR(1, 4, :)) - pRF(1);  txV = txV(:);
tyV = squeeze(ctx.T_t1_ICR(2, 4, :)) - pRF(2);  tyV = tyV(:);

% 1) Human target present and negative (flexion) across the flexion RoM
tst = interp1(ctx.humanAngleD(:), ctx.humanTorque(:), ...
    [-120 -90 -60 -30 0], 'pchip');
fprintf('human bifemsh target, N*m, at K = -120/-90/-60/-30/0 deg:\n');
fprintf('  %.3f  %.3f  %.3f  %.3f  %.3f\n', tst);
assert(tst(1) < 0 && tst(2) < 0 && tst(3) < 0 && tst(4) < 0, ...
    'human flexor target expected negative through flexion');

% 2) computeTorqueZ mirror: LC entered from the tibia axis, 90 - LC
pts = [ -80, 10, 100;    % Kdeg, Ldeg(tibia), F(N)
        -60, 15, 150;
       -110,  8,  60];
tz = zeros(size(pts, 1), 1);
for i = 1:size(pts, 1)
    Kr = deg2rad(pts(i, 1));
    Lr = deg2rad(90 - pts(i, 2));
    tx = interp1(phiV, txV, Kr, 'pchip');
    ty = interp1(phiV, tyV, Kr, 'pchip');
    Trk = RpToTrans(eye(3), [tx; ty; 0]);
    Fr = [0; 0; 0; pts(i,3)*cos(Lr + pi); pts(i,3)*sin(Lr + pi); 0];
    Fk = Adjoint(Trk)' * Fr;
    tz(i) = Fk(3);
end
fprintf('live-guard tz, N*m, at (K,L,F) = (-80,10,100) (-60,15,150) (-110,8,60):\n');
fprintf('  %.3f  %.3f  %.3f\n', tz);
assert(all(isfinite(tz)), 'nonfinite guard torque');
assert(all(tz < 0), 'flexor guard torque expected negative (flexion)');

% 3) Verdict semantics: MET iff tz <= tgt (both negative)
tgt80 = interp1(ctx.humanAngleD(:), ctx.humanTorque(:), -80, 'pchip');
if tz(1) <= tgt80
    v = 'MET';
else
    v = 'BELOW';
end
fprintf('at K = -80 deg: target %.3f N*m, tz %.3f N*m -> %s\n', ...
    tgt80, tz(1), v);

fprintf('ALL CHECKS PASSED\n');
end
