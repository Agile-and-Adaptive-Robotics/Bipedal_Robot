function result = fr_fitTravelLengths(C, cs, kin, design, branch, opts, maxGain)
% Refit two lengths with fixed attachments, size, count and series split.
% For fixed rest length, the BiPulley residual decreases with tendon
% length. Find the shortest tendon that clears slack, then check stroke.
% Only promising fits pay for the full torque and clearance evaluation.
assert(strcmp(cs.xi3Mode, 'perPose'), 'Length search assumes the current Xi3-free cable closure');
assert(max(cs.txConfig.gain)==1, 'Supply direct transmission as the starting topology');
if nargin < 7; maxGain = 3; end
n = numel(cs.movable);
lb = cs.lb(3*n+(1:2))/10;
ub = cs.ub(3*n+(1:2))/10;
optsFit = opts;
best = [];
attempts = struct([]);
passGrid = cell(1,numel(cs.passes));
lowerTendon = lb(2);
for p = 1:numel(cs.passes)
    [P,B] = fr_locForPass(C,cs,design,kin,p);
    poses = cs.passes{p}.qStructs;
    route = struct('points',P,'frames',{B});
    passGrid{p} = struct('passIdx',p,'qList',{poses},'resolvedRoute',route);
    if ~strcmp(cs.spanKind,'none')
        sp = fr_spans(C,kin,P,B,poses,cs.spanKind);
        lowerTendon = max(lowerTendon,sp.requiredSpan_m);
    end
end
coupled = [];
if isfield(opts,'coupledTravelN')
    coupled = fr_ankleCoupledSetup(C,kin,cs,opts.coupledTravelN);
    [P,B] = fr_locForPass(C,coupled.cs,design,kin,1);
    coupled.qg.resolvedRoute = struct('points',P,'frames',{B});
    sp = fr_spans(C,kin,P,B,coupled.cs.passes{1}.qStructs,'foot');
    lowerTendon = max(lowerTendon,sp.requiredSpan_m);
end
for gain = 1:maxGain
    candidate = cs;
    reeving = fr_travelReeving(gain);
    candidate.txConfig.gain = [gain 1];
    candidate.txConfig.tackleLineParts = [reeving.movingFalls 1];
    candidate.txConfig.reeving = reeving;
    candidate.hasTackle = gain > 1;
    if candidate.hasTackle
        candidate.txConfig.tackleFamily = 'L';
    else
        candidate.txConfig.tackleFamily = 'none';
    end
    if ~isempty(coupled)
        coupled.cs.txConfig = candidate.txConfig;
        coupled.cs.hasTackle = candidate.hasTackle;
    end
    % Include the warm rest length, the limits, and intermediate values.
    rests = unique([min(ub(1),max(lb(1),design.L0)), linspace(lb(1),ub(1),7)]);
    viable = zeros(0,2);
    trial = repmat(struct('rest_m',NaN,'tendon_m',NaN,'slack',NaN,'stroke',NaN),1,numel(rests));
    for j = 1:numel(rests)
        d = design; d.L0 = rests(j);
        lo = lowerTendon; hi = ub(2);
        trial(j).rest_m = d.L0;
        if lo > hi
            continue
        end
        [slack,stroke] = assess(lo);
        if slack > 0
            [slackHi,~] = assess(hi);
            if slackHi > 0
                trial(j).slack = slackHi;
                continue
            end
            while hi-lo > 0.00025
                mid = (hi+lo)/2;
                if assess(mid)>0
                    lo = mid;
                else
                    hi = mid;
                end
            end
            lo = hi;
            [slack,stroke] = assess(lo);
        end
        trial(j).tendon_m = lo;
        trial(j).slack = slack;
        trial(j).stroke = stroke;
        if slack == 0 && stroke == 0
            viable(end+1,:) = [d.L0,lo]; %#ok<AGROW>
        end
    end
    attempts(gain).gain = gain;
    attempts(gain).reeving = reeving;
    attempts(gain).lengthTrials = trial;
    attempts(gain).viableCount = size(viable,1);
    fprintf('%s %s gain %d: %d/%d sampled rest lengths close across all passes\n', ...
        cs.name,branch,gain,size(viable,1),numel(rests));
    for j = 1:size(viable,1)
        d = design; d.L0 = viable(j,1); d.tendon = viable(j,2);
        [f,rep] = fr_objective(C,candidate,kin,fr_packX(candidate,d),branch,optsFit);
        if isempty(best) || f < best.f
            best = struct('caseDef',candidate,'design',d,'rep',rep,'f',f,'reeving',reeving);
        end
    end
    if ~isempty(best)
        break
    end
    % Escalate only if a taut candidate reached the stroke boundary.
    % Slack-only and tendon-bound failures require geometry or length work.
    if ~any([trial.slack]==0 & [trial.stroke]>0)
        break
    end
end
result = struct('best',best,'attempts',attempts);

    function [slack,stroke] = assess(tendon)
        probe = d;
        probe.tendon = tendon;
        slack = 0; stroke = 0;
        for pass = 1:numel(candidate.passes)
            m = fr_mechBuild(C,candidate,probe,kin,branch,passGrid{pass});
            slack = slack + nnz(m.PulleySlack);
            stroke = stroke + nnz(m.PulleyInfeasible | ...
                m.strain_p > 1-m.kmax/m.rest+1e-9 | ~isfinite(m.strain_p));
        end
        if ~isempty(coupled)
            m = fr_mechBuild(C,coupled.cs,probe,kin,branch,coupled.qg);
            cellStrain = m.cell.sContraction(:)/m.rest;
            slack = slack + nnz(m.cell.PulleySlack);
            stroke = stroke + nnz(m.cell.PulleyInfeasible(:) | ...
                cellStrain > 1-m.kmax/m.rest+1e-9 | ~isfinite(cellStrain));
        end
    end
end
