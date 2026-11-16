function out=fr_compareAnkleLayers(C,kin,poses,outer,inner,fitOuter,fitInner,checkOrder)
% "Above" = superficial to EHL, away from the tibia, at shared shank y.
% Vessel overlap is reported separately from the hard fitting checks.
if nargin<8; checkOrder=true; end
order=inf; hardware=inf; overlap=0; matches=0;
for i=1:numel(poses)
    T=fr_bodyT(kin,poses{i},'tibia_r');
    for a=outer{i}
        A=local(T,[a.p1;a.p2]);
        for b=inner{i}
            B=local(T,[b.p1;b.p2]);
            low=max(min(A(:,2)),min(B(:,2)));
            high=min(max(A(:,2)),max(B(:,2)));
            if checkOrder && high-low>1e-4
                stations=linspace(low,high,max(2,ceil((high-low)/.01)+1))';
                x=atY(A,stations); y=atY(B,stations);
                [~,da]=fr_surfaceDistance(C,'tibia_r',x);
                [~,db]=fr_surfaceDistance(C,'tibia_r',y);
                order=min(order,min(da-db));
                matches=matches+1;
            end
            overlap=max(overlap,a.r+b.r-fr_segmentDistance(A(1,:),A(2,:),B(1,:),B(2,:)));
        end
    end
    [FA,ra]=fr_fittingSpheres(outer{i},fitOuter);
    [FB,rb]=fr_fittingSpheres(inner{i},fitInner);
    hardware=min(hardware,fittingToSleeve(FA,ra,inner{i}));
    hardware=min(hardware,fittingToSleeve(FB,rb,outer{i}));
    for j=1:size(FA,1)
        hardware=min(hardware,min(vecnorm(FB-FA(j,:),2,2)-ra(j)-rb));
    end
end
out=struct('orderingMargin_m',order,'rigidMargin_m',hardware-C.clearanceEnvelope_m, ...
    'sleeveOverlap_m',overlap,'overlapPairs',matches);
end

function p=local(T,p)
p=(T(1:3,1:3)'*(p'-T(1:3,4)))';
end

function p=atY(A,y)
p=A(1,:)+(y-A(1,2))/(A(2,2)-A(1,2)).*(A(2,:)-A(1,:));
end

function d=fittingToSleeve(p,r,caps)
d=inf;
for j=1:size(p,1)
    for c=caps
        d=min(d,fr_segmentDistance(p(j,:),p(j,:),c.p1,c.p2)-r(j)-c.r);
    end
end
end
