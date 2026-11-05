function fr_verify_repairs_20261105()
% Rotated planar fixture catches SVD-output mixups and doubled displacement.
C=fr_setup();root=fullfile(C.camp,'audit_20261105','verification_fixture');
if ~isfolder(root);mkdir(root);end
theta=.73;phi=.41;
R=[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]*[cos(phi),0,sin(phi);0,1,0;-sin(phi),0,cos(phi)];
[a,b]=meshgrid(linspace(-.06,.06,25));cloud=[a(:),b(:),zeros(numel(a),1)]*R';
writematrix(cloud,fullfile(root,'fixture_Mesh_Points.txt'),'Delimiter','space');
C.boneCloudDir=root;n=R(:,3)';seed=.01*n;
P=fr_patch(C,'fixture','fixture',seed,30);
assert(abs(dot(P.n,n))>1-1e-12,'Surface normal does not follow rotated plane');
assert(abs(norm(P.e1)-1)<1e-12&&abs(norm(P.e2)-1)<1e-12,'Non-unit surface directions');
p=P.project(.5,0,7);p0=P.project(0,0,7);
assert(abs(dot(p,n)-.007)<1e-12,'Incorrect surface standoff');
assert(abs(norm(p-p0)-.015)<1e-12,'Tangential displacement was doubled or scaled');
p=P.project(1,1,7);assert(norm(p-p0)<=.030000001,'Patch radius exceeded');
% The user-specified 20 mm calibration is exact at both endpoints.
C=fr_setup();assert(abs(fr_diameter(C,'20',0)-20)<1e-12);
assert(abs(fr_diameter(C,'20',.24)-38)<1e-12);
assert(abs(fr_segmentDistance([0,0,0],[2,0,0],[1,1,0],[3,1,0])-1)<1e-12);
assert(abs(fr_segmentDistance([0,0,0],[1,0,0],[.5,-1,1],[.5,1,1])-1)<1e-12);
assert(fr_segmentDistance([0,0,0],[1,0,0],[.5,-1,0],[.5,1,0])<1e-12);
assert(abs(fr_segmentDistance([-2,0,0],[-1,0,0],[0,0,0],[1,0,0])-1)<1e-12);
assert(abs(fr_segmentDistance([0,0,0],[0,0,0],[1,-1,0],[1,1,0])-1)<1e-12);
report=struct('rotated_surface_patch','PASS','unit_tangent_basis','PASS', ...
 'single_tangential_displacement','PASS','bounded_patch_radius','PASS','diameter20_endpoints','PASS', ...
 'segment_parallel','PASS','segment_skew','PASS','segment_crossing','PASS','segment_endpoint','PASS','segment_degenerate','PASS');
fid=fopen(fullfile(C.camp,'audit_20261105','repair_verification.json'),'w');fwrite(fid,jsonencode(report));fclose(fid);
disp(report);
end
