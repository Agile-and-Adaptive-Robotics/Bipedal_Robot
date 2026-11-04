function sched = fr_layers(C, availableCases)
%FR_LAYERS Layer schedule from the target contract's layerOrder.
%
% target_contract.json layerOrder (deep -> superficial):
%   trunk/deep hip: 27, 26, 1, 2, 3
%   thigh:          11, 10, 25, 5, 4, 8, 9, 6, 7, 24
%   lower leg:      15, 20, 21, 17, 18, 12, 13, 14, 16, 22, 23, 19
% Rule (contract): the deepest layer optimizes against bone; every later
% group treats ALL earlier layers' occupied envelopes as hard obstacles;
% coordinated refinement runs after the first pass.
%
% Stage-1 demonstrator subset: only the wired demonstrator cases are
% scheduled (documented in the gate log); the full 27-group campaign
% consumes the same schedule when the group cases land.
tc = jsondecode(fileread(fullfile(C.targetsDir, 'target_contract.json')));
order = [];
regions = fieldnames(tc.layerOrder);
regionOf = containers.Map();
for r = 1:numel(regions)
    ids = tc.layerOrder.(regions{r}).order;
    for i = 1:numel(ids)
        order(end+1) = ids(i); %#ok<AGROW>
        regionOf(ids(end)) = regions{r};
    end
end
sched = struct();
sched.ids = order;
sched.region = repmat({''}, 1, numel(order));
for i = 1:numel(order)
    sched.region{i} = regionOf(order(i));
end
% attach the demonstrator case names by map id
sched.caseByMapId = containers.Map('KeyType', 'double', 'ValueType', 'char');
for c = 1:numel(availableCases)
    sched.caseByMapId(availableCases(c).groupMapId) = availableCases(c).name;
end
sched.scheduledIds = order(ismember(order, cell2mat( ...
    sched.caseByMapId.keys)));
end
