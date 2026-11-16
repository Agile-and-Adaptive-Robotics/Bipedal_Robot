function L=fr_ankleShaftLead(C,kin,cs,loc,frames,poses)
% Check both free lead spans and the 3D cable turns at guide and stop.
guide=cs.distalSleeveGuideRow; stop=cs.shankStopRow;
foot=find(ismember(cs.rowOwning,{'talus_r','calcn_r','toes_r'}),1);
assert(guide<stop && stop+1==foot);
P=loc(guide-1:stop,:);
L=struct('shaftCableMargin_m',inf,'guideTurn_deg',0,'stopTurn_deg',0);
for i=1:size(P,1)-1
    if i==1; continue; end
    a=P(i,:); b=P(i+1,:); len=norm(b-a);
    if len<=.010; continue; end
    n=max(2,ceil(len/.002)+1);
    s=linspace(.004/len,1-.004/len,n)';
    d=fr_surfaceDistance(C,'tibia_r',a+s.*(b-a));
    L.shaftCableMargin_m=min(L.shaftCableMargin_m,min(d)-.002-len/(2*(n-1)));
end
L.guideTurn_deg=turn(P(2,:)-P(1,:),P(3,:)-P(2,:));
for j=1:numel(poses)
    T=fr_chainT(kin,poses{j},frames{foot},'tibia_r');
    dest=(T(1:3,1:3)*loc(foot,:)'+T(1:3,4))';
    L.stopTurn_deg=max(L.stopTurn_deg,turn(P(end,:)-P(end-1,:),dest-P(end,:)));
end
end

function a=turn(u,v)
cosine=dot(u,v)/max(norm(u)*norm(v),eps);
a=acosd(max(-1,min(1,cosine)));
end
