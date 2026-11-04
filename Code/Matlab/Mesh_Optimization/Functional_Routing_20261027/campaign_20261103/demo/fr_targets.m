function target = fr_targets(C, kin, cs)
%FR_TARGETS Human target torque curve from the STOCK human model (master
% directive item 8: targets come from the stock model's routes, never from
% the robot's optimized routes).
%
% tau_target(q) = sum over the target stock muscles of
%     Fmax_m * ( -dL_m/dq_primary(q) )
% with L_m(q) the muscle's STOCK route length computed by fr_evalRoute
% (declared conditional/moving behavior implemented) and dL/dq_primary the
% PARTIAL derivative wrt the case's primary swept coordinate (central
% differences, h = 0.5 deg) - the same virtual-work convention the BPA
% side is scored with, so the meet-or-beat comparison is sign-consistent.
% Fmax_m = the osim's own max_isometric_force.
%
% Also reports the peak vs the target contract's mifN (informational).

qList = cs.passes{1}.qStructs;
nG = numel(qList);
coord = cs.passes{1}.sweptCoord;
h = 0.5 * pi / 180;
tau = zeros(nG, 1);
FmaxUsed = zeros(numel(cs.stockMuscles), 1);
for m = 1:numel(cs.stockMuscles)
    rt = fr_muscleroute(C.osimStock, cs.stockMuscles{m});
    FmaxUsed(m) = fr_muscleProp(C.osimStock, cs.stockMuscles{m}, ...
        'max_isometric_force');
    for i = 1:nG
        qp = qList{i}; qp.(coord) = qp.(coord) + h;
        qm = qList{i}; qm.(coord) = qm.(coord) - h;
        Lp = fr_evalRoute(rt, kin, qp).lengthM;
        Lm = fr_evalRoute(rt, kin, qm).lengthM;
        tau(i) = tau(i) - FmaxUsed(m) * (Lp - Lm) / (2 * h);
    end
end
target = struct();
target.tau = tau;
target.qList = qList;
target.FmaxUsed = FmaxUsed;
target.muscles = cs.stockMuscles;
target.peakAbs = max(abs(tau));
% informational: contract mifN (stock group sum basis)
tc = jsondecode(fileread(fullfile(C.targetsDir, 'target_contract.json')));
target.contractMifN = NaN;
target.mapId = cs.groupMapId;
for g = 1:numel(tc.groups)
    if tc.groups(g).map_id == cs.groupMapId
        target.contractMifN = tc.groups(g).mifN;
    end
end
target.peakVsContractPct = 100 * (target.peakAbs / ...
    max(target.contractMifN, eps) - 1);
end
