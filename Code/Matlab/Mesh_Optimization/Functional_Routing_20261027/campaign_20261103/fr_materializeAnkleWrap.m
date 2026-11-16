function [cs,d]=fr_materializeAnkleWrap(C,kin,cs,d)
% Round the tibial control point into a finite-radius guided vessel bend.
% A 75 mm screening radius is a trial geometry assumption, not a rating.
row=cs.bpaWrapControlRow;
[P,~,ids]=fr_locForPass(C,cs,d,kin,1,true);
A=P(ids==1,:); W=P(ids==row,:); B=P(ids==cs.distalSleeveGuideRow,:);
u=W-A; v=B-W;
L=min(norm(u),norm(v));
assert(L>.015,'Wrap guide is too close to a mount');
u=u/norm(u); v=v/norm(v);
cosine=min(1,max(-1,dot(u,v)));
theta=acos(cosine);
radiusTarget=.075;
if theta<1e-4
    arc=W; radius=inf; arcLength=0;
else
    cut=min(radiusTarget*tan(theta/2),.45*L);
    radius=cut/tan(theta/2);
    normal=(v-cosine*u)/sin(theta);
    tangent=W-cut*u;
    center=tangent+radius*normal;
    angle=linspace(0,theta,7)';
    arc=center+radius*(sin(angle).*u-cos(angle).*normal);
    arcLength=radius*theta;
end
index=find([cs.movable.row]==row,1);
cs.movable(index)=[];
d.uv(index,:)=[];
remove=3*index-2:3*index;
cs.lb(remove)=[]; cs.ub(remove)=[];
cs.nX=cs.nX-3;
cs.seedLoc(row,:)=arc(1,:);
for j=2:size(arc,1)
    cs=fr_addViaGuide(cs,row+j-2,arc(j,:));
end
added=size(arc,1)-1;
cs.distalSleeveGuideRow=cs.distalSleeveGuideRow+added;
cs.shankStopRow=cs.shankStopRow+added;
cs.bpaBendArc_m=arcLength;
cs.bpaBendRadius_m=radius;
cs.bpaBendRadiusTarget_m=radiusTarget;
cs.bpaBendEndRow=row+added;
cs=rmfield(cs,'bpaWrapControlRow');
cs.xSeed=fr_packX(cs,d);
end
