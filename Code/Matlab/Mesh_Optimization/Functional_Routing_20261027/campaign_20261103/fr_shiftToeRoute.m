function [cs,shift_m]=fr_shiftToeRoute(cs,s,requested_mm)
% Project the insertion and dependent toe guides into their own patches.
shift_m=zeros(numel(s.rows),1);
if requested_mm<=0; return; end
points=s.points;
for i=1:numel(s.rows)
    p=s.patch{i};
    tangent=s.distal-dot(s.distal,p.n)*p.n;
    tangent=tangent/max(norm(tangent),eps);
    uv=requested_mm/1000/p.Rpatch_m*[dot(tangent,p.e1),dot(tangent,p.e2)];
    seed=s.points(i,:);
    standoff=1000*dot(seed-p.project(0,0,0),p.n);
    point=p.project(uv(1),uv(2),standoff);
    displacement=point-seed;
    if dot(displacement,s.distal)<=0 || norm(displacement)>.015001
        shift_m(:)=0;
        return
    end
    points(i,:)=point;
    shift_m(i)=dot(displacement,s.distal);
end
cs.seedLoc(s.rows,:)=points;
cs.allowFootLeadCorrection=false;
end
