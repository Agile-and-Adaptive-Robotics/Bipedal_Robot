function [p,r,branch]=fr_fittingSpheres(caps,L)
% Conservative end-fitting envelopes for each parallel branch and section.
p=zeros(0,3); r=zeros(0,1); branch=zeros(0,1);
keys=unique([[caps.branch]',[caps.section]'],'rows');
for i=1:size(keys,1)
    subset=caps([caps.branch]==keys(i,1) & [caps.section]==keys(i,2));
    a=subset(1); b=subset(end);
    u=a.p2-a.p1; v=b.p2-b.p1;
    p=[p;a.p1-L/2*u/max(norm(u),eps);b.p2+L/2*v/max(norm(v),eps)]; %#ok<AGROW>
    r=[r;hypot(a.r,L/2);hypot(b.r,L/2)]; %#ok<AGROW>
    branch=[branch;keys(i,1);keys(i,1)]; %#ok<AGROW>
end
end
