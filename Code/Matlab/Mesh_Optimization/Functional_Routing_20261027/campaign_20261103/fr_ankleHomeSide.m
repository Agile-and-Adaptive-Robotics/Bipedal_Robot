function z=fr_ankleHomeSide(kin,r)
% Right-side attachment and guide locations in the home pelvis frame only.
z=inf;
for k=1:size(r.p,1)
    T=fr_chainT(kin,kin.qDefault,r.caseDef.rowOwning{k},'pelvis');
    p=T(1:3,1:3)*r.p(k,:)'+T(1:3,4);
    z=min(z,p(3));
end
end
