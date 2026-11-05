function d=fr_surfaceDistance(C, body, pointsLocal)
% Signed distance to closed OpenSim meshes in their own body frame.
persistent cache
if isempty(cache);cache=containers.Map('KeyType','char','ValueType','any');end
root=fullfile(C.camp,'audit_20261105','surface_fields');key=[root '|' body];
if ~isKey(cache,key)
 f=fullfile(root,[body '.mat']);assert(isfile(f),'fr_surfaceDistance:Missing','No verified closed surface for %s',body);
 s=load(f);s.interp=griddedInterpolant({s.x(:),s.y(:),s.z(:)},double(s.distance),'linear','none');cache(key)=s;
end
s=cache(key);p=pointsLocal;
d=s.interp(p(:,1),p(:,2),p(:,3))-s.uncertainty_m;
outside=isnan(d);
if any(outside)
 bb=s.mesh_bbox;
 delta=max(max(bb(1,:)-p(outside,:),p(outside,:)-bb(2,:)),0);
 d(outside)=vecnorm(delta,2,2); % distance to containing box is a lower bound
end
end
