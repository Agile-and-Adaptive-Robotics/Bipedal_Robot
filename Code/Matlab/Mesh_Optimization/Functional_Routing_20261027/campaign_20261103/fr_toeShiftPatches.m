function s=fr_toeShiftPatches(C,kin,cs)
% Prepare an optional distal shift for the toe insertion and its toe guides.
assert(strcmp(cs.spanKind,'foot'));
s.rows=find(strcmp(cs.rowOwning,'toes_r'));
assert(~isempty(s.rows) && s.rows(end)==size(cs.seedLoc,1));
T=fr_chainT(kin,kin.qDefault,'toes_r','pelvis');
s.distal=T(1,1:3);
s.points=cs.seedLoc(s.rows,:);
s.patch=cell(size(s.rows));
for i=1:numel(s.rows)
    s.patch{i}=fr_patch(C,'Toes','toes_r',s.points(i,:),15);
end
end
