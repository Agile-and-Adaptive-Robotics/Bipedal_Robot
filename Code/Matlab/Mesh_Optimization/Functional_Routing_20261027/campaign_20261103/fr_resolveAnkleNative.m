function G=fr_resolveAnkleNative(C,kin,cs,p,n,angles)
% Resolve explicit independent cable paths in the current crossing frames.
if nargin<6
    G=fr_ankleCoupledSetup(C,kin,cs,n);
else
    G=fr_ankleCoupledSetup(C,kin,cs,n,angles);
end
frames=G.cs.rowFrame{1};
assert(size(p,1)==numel(frames));
for r=1:numel(frames)
    T=fr_chainT(kin,kin.qDefault,cs.rowOwning{r},frames{r});
    p(r,:)=(T(1:3,1:3)*p(r,:)'+T(1:3,4))';
end
G.qg.resolvedRoute=struct('points',p,'frames',{frames});
end
