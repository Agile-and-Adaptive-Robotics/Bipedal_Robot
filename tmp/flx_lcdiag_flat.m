%FLX_LCDIAG_FLAT Flat LC-convention diagnostic (no subfunctions).
% Uses base-workspace ctx (buildKneeFlexorContext20mm) if present.
if ~exist('ctx', 'var') || ~isfield(ctx, 'T_t1_ICR')
    evalc('ctx = buildKneeFlexorContext20mm();');
end
T_t1_ICR = ctx.T_t1_ICR;
phi = ctx.phi(:);
c = pi/180;

%% (a) 42 cm measured data, exact old-code math from the data file
Load = [18.5 18.05 32 40.97 44.4 50.84 62.45 69.4 64.3 70.6 90.15 70];
K_ang = [-125 -114 -98 -83 -75.5 -69 -55.5 -53.001 3 7 -6.5 -32]*c;
LC_ang = [32.5 30 28 26 24.5 17 20 24 7 5.5 10 13.5]*c;
d = 320/1000; ang = -82.97;
p_rf = [d*cosd(ang), d*sind(ang), 0]';
T_t1_rf = RpToTrans(eye(3), p_rf);
Trk1 = pagemtimes(TransInv(T_t1_rf), T_t1_ICR);
fcnX1 = fit(phi, squeeze(Trk1(1,4,:)), 'cubicspline');
fcnY1 = fit(phi, squeeze(Trk1(2,4,:)), 'cubicspline');
Tz1 = zeros(1, numel(Load));
for i = 1:numel(Load)
    px = fcnX1(K_ang(i)); py = fcnY1(K_ang(i));
    Trki = RpToTrans(eye(3), [px, py, 0]');
    Fri = [0; 0; 0; Load(i)*cos(LC_ang(i)+pi); Load(i)*sin(LC_ang(i)+pi); 0];
    Fkv = Adjoint(Trki)'*Fri; Tz1(i) = Fkv(3);
end

Load2 = [91.7 114.1 129.8 78 92.78 63.6 97.28 71.92 84.5 93.9];
K_ang2 = [-53.01 -41 -30 -26.01 -26.001 -18.5 -18 -7 -9 0]*c;
LC_ang2 = [-5 -7 -12 -12 -12.5 -16.5 -15 -19.5 -15 -17]*c;
d2 = 218.29/1000; ang2 = -79.84;
p_rf2 = [d2*cosd(ang2), d2*sind(ang2), 0]';
Trk2 = pagemtimes(TransInv(RpToTrans(eye(3), p_rf2)), T_t1_ICR);
fcnX2 = fit(phi, squeeze(Trk2(1,4,:)), 'cubicspline');
fcnY2 = fit(phi, squeeze(Trk2(2,4,:)), 'cubicspline');
Tz2 = zeros(1, numel(Load2));
for i = 1:numel(Load2)
    px = fcnX2(K_ang2(i)); py = fcnY2(K_ang2(i));
    Trki = RpToTrans(eye(3), [px, py, 0]');
    Fri = [0; 0; 0; Load2(i)*cos(LC_ang2(i)+pi); Load2(i)*sin(LC_ang2(i)+pi); 0];
    Fkv = Adjoint(Trki)'*Fri; Tz2(i) = Fkv(3);
end

fprintf('\n=== 42 cm tibia-1 (old code, entries 5.5..32.5) ===\n');
fprintf('   K_deg   Load_N   LC_deg   TorqueZ   |arm| m\n');
for i = 1:numel(Load)
    fprintf('%8.1f %8.2f %8.1f %9.3f %8.4f\n', K_ang(i)/c, Load(i), ...
        LC_ang(i)/c, Tz1(i), abs(Tz1(i))/Load(i));
end
fprintf('\n=== 42 cm tibia-2 (old code, entries -5..-19.5) ===\n');
for i = 1:numel(Load2)
    fprintf('%8.1f %8.2f %8.1f %9.3f %8.4f\n', K_ang2(i)/c, Load2(i), ...
        LC_ang2(i)/c, Tz2(i), abs(Tz2(i))/Load2(i));
end

%% (b) live-guard sweep, flat
dLC = 320/1000; angLC = -82.97;
pRF = [dLC*cosd(angLC), dLC*sind(angLC), 0];
txV = squeeze(T_t1_ICR(1,4,:)) - pRF(1); txV = txV(:);
tyV = squeeze(T_t1_ICR(2,4,:)) - pRF(2); tyV = tyV(:);

Es = [0 2.83 5 10 15 20 27 30 45 60 75 90 105 120 -10 -27 -30];
Ks = [-100 -80 -60 -40];
F = 80;
convNames = {'current 90-LC,+pi', 'E-2.83,no pi', 'E-2.83,+pi', ...
    'E,no pi', 'E,+pi', '90-LC,no pi'};

convPi = [true false true false true false];

for cv = 1:6
    fprintf('\n=== conv %d: %s (F = %d N) ===\n', cv, convNames{cv}, F);
    fprintf('  K\\LC');
    fprintf(' %7.1f', Es);
    fprintf('\n');
    for iK = 1:numel(Ks)
        Kr = deg2rad(Ks(iK));
        tx = interp1(phi, txV, Kr, 'pchip');
        ty = interp1(phi, tyV, Kr, 'pchip');
        Trk0 = RpToTrans(eye(3), [tx; ty; 0]);
        fprintf('%7.0f', Ks(iK));
        for iE = 1:numel(Es)
            Lr = switchLr(cv, Es(iE));
            if convPi(cv)
                Frv = [0; 0; 0; F*cos(Lr+pi); F*sin(Lr+pi); 0];
            else
                Frv = [0; 0; 0; F*cos(Lr); F*sin(Lr); 0];
            end
            Fkv = Adjoint(Trk0)'*Frv; tz = Fkv(3);
            fprintf(' %7.2f', tz);
        end
        fprintf('\n');
    end
end

function Lr = switchLr(cv, Ldeg)
switch cv
    case 1
        Lr = deg2rad(90 - Ldeg);
    case 2
        Lr = deg2rad(Ldeg - 2.83);
    case 3
        Lr = deg2rad(Ldeg - 2.83);
    case 4
        Lr = deg2rad(Ldeg);
    case 5
        Lr = deg2rad(Ldeg);
    case 6
        Lr = deg2rad(90 - Ldeg);
end
end
