function positions=fr_parallelRows(count,pitch)
% Flat, centered rows of lengthwise vessels, in a cross-section of the limb.
width=min(count,ceil(sqrt(2*count)));
positions=zeros(count,2);
for first=1:width:count
    ids=first:min(count,first+width-1);
    positions(ids,1)=pitch*((0:numel(ids)-1)-(numel(ids)-1)/2);
    positions(ids,2)=pitch*floor((first-1)/width);
end
positions=positions-mean(positions,1);
end
