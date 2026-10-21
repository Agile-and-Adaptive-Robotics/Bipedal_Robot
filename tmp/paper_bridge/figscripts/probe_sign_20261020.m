% probe_sign_20261020.m -- hand-verify the BiPam_pulley crossing-1 torque
% sign for soleus_r vs tib_ant_r from the recorded mats. Read-only.
root = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');

sn = {'soleus_r', 'tib_ant_r', 'bifemsh_r'};
for s = 1:numel(sn)
    f = dir(fullfile(resDir, sprintf('BiPulley_opensim_%s_*.mat', sn{s})));
    S = load(fullfile(resDir, f(1).name), 'spec', 'T', 'points', 'finalObj');
    sp = S.spec;
    obj = BiPam_pulley(sp.name, S.points, sp.cross, 20, S.T, ...
        S.finalObj.RestingL, S.finalObj.Kmax, S.finalObj.TendonL, 0.025, 620, ...
        0.005, 5e5, 1e4, 6, struct('nPulleyBPA', 1, ...
        'tackleLineParts', [1 1], 'routingMode', 'moving_via'));
    ax1 = sp.kin.axes(1, :) / norm(sp.kin.axes(1, :));
    fprintf('\n== %s: rows of points (frame order), cross=[%d %d]\n', sn{s}, sp.cross);
    for r = 1:size(S.points, 1)
        fprintf('   row %d: [%8.4f %8.4f %8.4f]  frame %d\n', r, S.points(r, :), ...
            biRowFrame(sp.cross, r));
    end
    % mid-grid cell: hand-compute crossing-1 torque in the child frame
    ii = 3;
    T1 = S.T(:, :, ii, 1);
    pE = S.points(1, :);                       % exit row (proximal frame)
    pI = S.points(end, :);                     % insertion row
    % map both into crossing-1 child frame (frame 2)
    pE2 = RowVecTrans(inv(T1), pE);
    pI2 = S.points(end, :);                    % frame 3? check row frame
    if biRowFrame(sp.cross, size(S.points, 1)) == 3
        T2 = T1;   % theta2 grid may be identity; use T(:,:,1,2)
        % (distal rows on frame 3 -> map 3->2 via T2 inverse of the joint-2 slab)
        T2s = S.T(:, :, 1, 2);
        pI2 = RowVecTrans(inv(T2s), S.points(end, :));
        % careful: T(:,:,1,2) maps frame3->frame2; a frame-3 point into frame 2
        % needs the INVERSE of that map.
    end
    u = pE2 - pI2;  u = u / norm(u);
    tauHand = S.finalObj.Ftendon(ii, 1) * dot(cross(pI2, u), ax1);
    tauClass = sum(squeeze(obj.Torque_ins(ii, :, 1, 1)) .* ax1, 2);
    fprintf('cell ii=3: hand tau = %.3f | class tau = %.3f | Ftendon=%.1f N\n', ...
        tauHand, tauClass, S.finalObj.Ftendon(ii));
    % also: dL/dtheta1 sign vs tau (OpenSim: moment arm = -dL/dq)
    mL = obj.MuscleLength;          % N1 x N2 rigid
    dLdth = (mL(2, 1) - mL(1, 1)) / (S.T(1, 1, 2, 1) - S.T(1, 1, 1, 1));
    fprintf('dL/dtheta1(1) = %.5f m/rad ; class tau(1) = %.3f ; -F*dL/dq = %.3f\n', ...
        dLdth, sum(squeeze(obj.Torque_ins(1, :, 1, 1)) .* ax1, 2), ...
        -S.finalObj.Ftendon(1) * dLdth);
end

function f = biRowFrame(cross, r)
if r < cross(1); f = 1; elseif cross(1) == cross(2); f = 3; elseif r < cross(2); f = 2; else; f = 3; end
end
