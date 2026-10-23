%probe_ext20_variants_20261023.m -- find the adjoint convention/reaction-point
% combination that reproduces the task-stated six ExtTest20mm_2 torques.
root = 'D:/GitHub/Bipedal_Robot';
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
addpath(fullfile(root, 'Code/Matlab/Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data/2022_02_Festo'), '-end');
warning('off', 'all');
cdeg = pi/180;

Load6 = [180.2 181.5 168 181 116 100];
K6    = [-26.2 -30 -61 -43 -119 5];
LC6   = [28 26 23 27 10.5 21];
taskT = [37.23 38.85 37.99 38.31 30.13 22.86];   % at angles [-119 -61 -43 -26.2 -30 5]
angTask = [-119 -61 -43 -26.2 -30 5];

C = load(fullfile(root, 'Code/Matlab/Mesh_Optimization/Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), 'ctxCamBest');
ctxA = C.ctxCamBest;
ctxB = []; %#ok<NASGU>
try
    evalc('ctxB = buildKneeExtContext20mm();');
catch
end
ctxs = {ctxA};
if exist('ctxB', 'var') && ~isempty(ctxB)
    ctxs{end+1} = ctxB;
end

prfs = {[-0.065 -0.255 0]', [-0.050 -0.255 0]'};
convs = {'compare (cosLC,-sinLC)', 'flexor (-cosLC,-sinLC)', 'adjoint22 (sinLC,-cosLC)'};

K_ang = K6 * cdeg;
for ci = 1:numel(ctxs)
    ctx = ctxs{ci};
    phi = ctx.phi(:);
    for pi_ = 1:numel(prfs)
        p_rf = prfs{pi_};
        T_t1_rf = RpToTrans(eye(3), p_rf);
        Trk = pagemtimes(TransInv(T_t1_rf), ctx.T_t1_ICR);
        f15 = fit(phi, squeeze(Trk(1,4,:)), 'cubicspline');
        f16 = fit(phi, squeeze(Trk(2,4,:)), 'cubicspline');
        tx = arrayfun(@(a) f15(a), K_ang);
        ty = arrayfun(@(a) f16(a), K_ang);
        for cv = 1:3
            LC = LC6 * cdeg;
            switch cv
                case 1, Fx = Load6.*cos(LC);  Fy = -Load6.*sin(LC);
                case 2, Fx = -Load6.*cos(LC); Fy = -Load6.*sin(LC);
                case 3, Fx = Load6.*sin(LC);  Fy = -Load6.*cos(LC);
            end
            tau = ty .* Fx - tx .* Fy;   % pure-translation adjoint tau_z
            [~, iOrd] = sort(K6);
            tauTask = tau(arrayfun(@(a) find(K6 == a, 1), angTask));
            dRaw = max(abs(tauTask - taskT));
            dFl = max(abs(-tauTask - taskT));
            fprintf('ctx%d prf(%.3f,%.3f) %-26s maxdiff raw %7.3f flipped %7.3f\n', ...
                ci, p_rf(1), p_rf(2), convs{cv}, dRaw, dFl);
            if min(dRaw, dFl) < 0.5
                fprintf('   -> MATCH-CLASS: '); fprintf('%8.3f', tauTask); fprintf('\n');
            end
        end
    end
    % also print the t positions once
    if ci == 1
        fprintf('tx (m): '); fprintf('%9.4f', tx); fprintf('\n');
        fprintf('ty (m): '); fprintf('%9.4f', ty); fprintf('\n');
    end
end
disp('DONE_variants');
