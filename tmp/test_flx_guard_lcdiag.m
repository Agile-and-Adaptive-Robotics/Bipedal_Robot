function test_flx_guard_lcdiag()
%TEST_FLX_GUARD_LCDIAG Diagnose the live-guard LC angle convention.
% Ben reports: 80 N at ~1/3 m arm should give ~-27 N*m; the guard shows
% -12. This script (a) reproduces the 42 cm flexor measured torques with
% the exact old-code math from Knee_Flexor_Data_20mm.m, and (b) sweeps
% the live guard's force-vector convention over entered LC values to
% find which one matches the rig's protractor zero and flexion sign.
% No serial port is opened.

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
T_t1_ICR = ctx.T_t1_ICR;
phi = ctx.phi(:);

%% (a) 42 cm measured data, EXACT old-code math (of-record values)
c = pi/180;
Load = [18.5 18.05 32 40.97 44.4 50.84 62.45 69.4 64.3 70.6 90.15 70];
K_ang = [-125 -114 -98 -83 -75.5 -69 -55.5 -53.001 3 7 -6.5 -32]*c;
LC_ang = [32.5 30 28 26 24.5 17 20 24 7 5.5 10 13.5]*c;
TorqueZ1 = oldAdjoint(Load, K_ang, LC_ang, 320/1000, -82.97, ...
    T_t1_ICR, phi);

Load2 = [91.7 114.1 129.8 78 92.78 63.6 97.28 71.92 84.5 93.9];
K_ang2 = [-53.01 -41 -30 -26.01 -26.001 -18.5 -18 -7 -9 0]*c;
LC_ang2 = [-5 -7 -12 -12 -12.5 -16.5 -15 -19.5 -15 -17]*c;
TorqueZ2 = oldAdjoint(Load2, K_ang2, LC_ang2, 218.29/1000, -79.84, ...
    T_t1_ICR, phi);

fprintf('\n=== 42 cm tibia-1 (old code, entries 5.5..32.5, +pi) ===\n');
fprintf('   K_deg   Load_N   LC_deg   TorqueZ   |arm| m\n');
for i = 1:numel(Load)
    fprintf('%8.1f %8.2f %8.1f %9.3f %8.4f\n', K_ang(i)/c, Load(i), ...
        LC_ang(i)/c, TorqueZ1(i), abs(TorqueZ1(i))/Load(i));
end
fprintf('\n=== 42 cm tibia-2 (old code, entries -5..-19.5, +pi) ===\n');
for i = 1:numel(Load2)
    fprintf('%8.1f %8.2f %8.1f %9.3f %8.4f\n', K_ang2(i)/c, Load2(i), ...
        LC_ang2(i)/c, TorqueZ2(i), abs(TorqueZ2(i))/Load2(i));
end

%% (b) live-guard convention sweep
% Guard kinematics with the 42 cm reaction point (what the live guard
% currently uses): tx/ty = ICR - pRF, pchip interp on phi.
dLC = 320/1000; angLC = -82.97;
pRF = [dLC*cosd(angLC), dLC*sind(angLC), 0];
txV = squeeze(T_t1_ICR(1, 4, :)) - pRF(1);  txV = txV(:);
tyV = squeeze(T_t1_ICR(2, 4, :)) - pRF(2);  tyV = tyV(:);

convs = { ...
    'current (90-LC, +pi)',   @(Ldeg) deg2rad(90 - Ldeg), true; ...
    'E-2.83, no pi',          @(Ldeg) deg2rad(Ldeg - 2.83), false; ...
    'E-2.83, +pi',            @(Ldeg) deg2rad(Ldeg - 2.83), true; ...
    'E, no pi',               @(Ldeg) deg2rad(Ldeg), false; ...
    'E, +pi',                 @(Ldeg) deg2rad(Ldeg), true; ...
    '90-LC, no pi',           @(Ldeg) deg2rad(90 - Ldeg), false};

Es = [0 2.83 5 10 15 20 27 30 45 60 75 90 105 120 -10 -27 -30];
Ks = [-100 -80 -60 -40];
F = 80;

for cv = 1:size(convs, 1)
    fprintf(['\n=== conv: %s  (F = %d N, |r| ~ 0.34 m at K = -80) ' ...
        '===\n'], convs{cv, 1}, F);
    fprintf('  K\\LC');
    fprintf(' %7.1f', Es);
    fprintf('\n');
    for iK = 1:numel(Ks)
        fprintf('%7.0f', Ks(iK));
        for iE = 1:numel(Es)
            tz = guardTz(Ks(iK), Es(iE), F, convs{cv, 2}, convs{cv, 3}, ...
                phi, txV, tyV);
            fprintf(' %7.2f', tz);
        end
        fprintf('\n');
    end
end

%% (c) What entry gives Ben's reported -12, and what gives ~-27?
fprintf('\n=== current conv: entry that reproduces the live -12 ===\n');
for iK = 1:numel(Ks)
    for iE = 1:numel(Es)
        tz = guardTz(Ks(iK), Es(iE), F, convs{1, 2}, convs{1, 3}, ...
            phi, txV, tyV);
        if abs(tz - (-12)) < 1.5
            fprintf('K = %d deg, LC entry = %.2f -> tz = %.2f\n', ...
                Ks(iK), Es(iE), tz);
        end
    end
end
end

function tz = guardTz(Kdeg, Ldeg, F, toLr, usePi, phi, txV, tyV)
    Kr = deg2rad(Kdeg);
    Lr = toLr(Ldeg);
    tx = interp1(phi, txV, Kr, 'pchip');
    ty = interp1(phi, tyV, Kr, 'pchip');
    if ~isfinite(tx) || ~isfinite(ty)
        tz = NaN;
        return;
    end
    Trk = RpToTrans(eye(3), [tx; ty; 0]);
    if usePi
        Fr = [0; 0; 0; F*cos(Lr + pi); F*sin(Lr + pi); 0];
    else
        Fr = [0; 0; 0; F*cos(Lr); F*sin(Lr); 0];
    end
    Fk = Adjoint(Trk)' * Fr;
    tz = Fk(3);
end

function TorqueZ = oldAdjoint(Load, K_ang, LC_ang, d, ang, T_t1_ICR, phi)
% EXACT block from Knee_Flexor_Data_20mm.m "Compare to results" section.
p_rf = [d*cosd(ang), d*sind(ang), 0]';
T_t1_rf = RpToTrans(eye(3), p_rf);
Trk = pagemtimes(TransInv(T_t1_rf), T_t1_ICR);
fcnX = fit(phi, squeeze(Trk(1, 4, :)), 'cubicspline');
fcnY = fit(phi, squeeze(Trk(2, 4, :)), 'cubicspline');
n = numel(Load);
Fk = zeros(6, 1, n);
for i = 1:n
    Trki = RpToTrans(eye(3), [fcnX(K_ang(i)), fcnY(K_ang(i)), 0]');
    Fri = [0; 0; 0; Load(i)*cos(LC_ang(i) + pi); Load(i)*sin(LC_ang(i) + pi); 0];
    Fk(:, :, i) = Adjoint(Trki)' * Fri;
end
TorqueZ = squeeze(Fk(3, 1, :));
end
