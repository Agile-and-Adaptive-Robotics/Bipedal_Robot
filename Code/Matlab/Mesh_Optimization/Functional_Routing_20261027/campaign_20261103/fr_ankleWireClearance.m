function margin=fr_ankleWireClearance(kin,poses,a,stopA,b,stopB)
% Minimum centerline separation of two distinct distal wires, less 2 mm.
margin=inf;
for i=1:numel(poses)
    A=world(kin,poses{i},a,stopA);
    B=world(kin,poses{i},b,stopB);
    for j=1:size(A,1)-1
        for k=1:size(B,1)-1
            margin=min(margin,fr_segmentDistance(A(j,:),A(j+1,:),B(k,:),B(k+1,:))-.002);
        end
    end
end
end

function p=world(kin,pose,d,stop)
rows=stop:size(d.loc,1);
p=zeros(numel(rows),3);
for j=1:numel(rows)
    r=rows(j);
    T=fr_bodyT(kin,pose,d.frames{r});
    p(j,:)=(T(1:3,1:3)*d.loc(r,:)'+T(1:3,4))';
end
end
