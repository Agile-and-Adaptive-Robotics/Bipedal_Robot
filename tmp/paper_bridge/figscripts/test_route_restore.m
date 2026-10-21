%test_route_restore.m — does the commented 10mm routing (with the script's
%current BB/CC activation logic) reproduce the ARCHIVED mat's Location?
root = 'D:/GitHub/Bipedal_Robot';
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
festo = fullfile(root, 'Testing_Data/2022_02_Festo');
A = load(fullfile(festo, 'KneeExt_10mm_52cm_stale_20261019.mat'));
phiD = A.phiD(:)'; positions = numel(phiD);
T_ICR_t1 = A.T_ICR_t1;

% commented 10 mm routing from Knee_Extensor_10mm.m
p1 = [0.040, 0.035, 0];
p2 = [0.0759, -0.27476, 0];
p3 = [0.05982, -0.37427, 0.000];
p4 = [0.03955, -0.42183, 0.000];
p5 = [0.05871, 0.025, 0];
p6 = [0.05871, 0.01228, 0];
p7 = [0.054, -0.00612, 0];
p8 = [0.03604, -0.02844, 0];

BB = -69.8; CC = -9.19;
Location = zeros(8,3,positions);
for i = 1:positions
    v5 = RowVecTrans(T_ICR_t1(:,:,i), p5);
    v6 = RowVecTrans(T_ICR_t1(:,:,i), p6);
    v7 = RowVecTrans(T_ICR_t1(:,:,i), p7);
    v8 = RowVecTrans(T_ICR_t1(:,:,i), p8);
    if phiD(i) > CC
        Location(:,:,i) = [p1; p2; p2; p2; v6; v6; v7; v8];
    elseif phiD(i) <= CC && phiD(i) > BB
        Location(:,:,i) = [p1; p2; p3; p3; v5; v6; v7; v8];
    else
        Location(:,:,i) = [p1; p2; p3; p4; v5; v6; v7; v8];
    end
end
d = max(abs(Location(:) - A.Location(:)));
fprintf('ROUTE max|rebuilt - archived| = %.3e ; isequal = %d\n', d, isequal(Location, A.Location));
fprintf('OLD adj TendonL = %.4f (plain-class tendon hack?)\n', A.Vas_Pam_adj.TendonL);
disp('DONE_routetest');
