function [margin,witness] = fr_tibialLeadMargin(C, kin, pts, frames, owning, poses, clearance)
% Margin of the free cable from the last tibial point to the first foot point.
% Use the tibia field in its own frame, including throughout ankle motion.
margin = inf;
witness = [];
foot = find(ismember(owning, {'talus_r', 'calcn_r', 'toes_r'}), 1);
if isempty(foot) || foot < 2 || ~strcmp(owning{foot-1}, 'tibia_r')
    return
end
for i = 1:numel(poses)
    q = poses{i};
    A = fr_chainT(kin, q, frames{foot-1}, 'tibia_r');
    B = fr_chainT(kin, q, frames{foot}, 'tibia_r');
    a = (A(1:3,1:3)*pts(foot-1,:)' + A(1:3,4))';
    b = (B(1:3,1:3)*pts(foot,:)' + B(1:3,4))';
    L = norm(b-a);
    % The stop and the attachment may touch their mounting surfaces.
    % Check the free span outside the first and last 2 mm.
    if L <= 0.004
        continue
    end
    n = ceil((L-0.004)/0.001)+1;
    s = linspace(0.002/L, 1-0.002/L, n)';
    samples=a+s.*(b-a);
    [conservative, nominal] = fr_surfaceDistance(C, 'tibia_r', samples);
    if clearance > 0
        [d,closest] = min(conservative);
        d = d - (L-0.004)/(2*(n-1)) - clearance;
    else
        [d,closest] = min(nominal);
    end
    if nargout>1 && d<margin
        witness=struct('poseIndex',i,'segment_m',[a;b], ...
            'closest_m',samples(closest,:),'nominalDistance_m',nominal(closest), ...
            'margin_m',d);
    end
    margin = min(margin, d);
end
end
