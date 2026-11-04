function [layerIds, regionName] = fr_layerIds(tc, CASES)
%FR_LAYERIDS Contract layerOrder -> {region id lists, display names}.
%
% jsondecode maps 'trunk/deep hip' to the field trunk_deepHip and keeps
% the NON-region string fields (basis, rule) in the same struct, so the
% regions are the STRUCT-VALUED fields carrying an .order array, in JSON
% order. Display names come from the groups' own region field (all ids
% within one region share it).
lo = tc.layerOrder;
fn = fieldnames(lo);
regOf = containers.Map('KeyType', 'double', 'ValueType', 'char');
for c = 1:numel(CASES)
    regOf(CASES(c).groupMapId) = CASES(c).region;
end
layerIds = {};
regionName = {};
for r = 1:numel(fn)
    v = lo.(fn{r});
    if ~isstruct(v) || ~isfield(v, 'order'); continue; end
    % jsondecode returns JSON arrays as COLUMN vectors - force row
    layerIds{end + 1} = reshape(double(v.order), 1, []);
    if isKey(regOf, v.order(1))
        regionName{end + 1} = regOf(v.order(1));
    else
        regionName{end + 1} = fn{r};
    end
end
end
