function export_cam_label_geometry(outputFile)
% Read the selected cam result and export each plotted active p-row in metres.
root = fileparts(mfilename('fullpath'));
while ~isfolder(fullfile(root,'Code','Matlab','Mesh_Optimization'))
    parent = fileparts(root);
    assert(~strcmp(parent,root),'Bipedal_Robot root was not found.');
    root = parent;
end
if nargin < 1
    outputFile = fullfile(fileparts(mfilename('fullpath')),'CAM_LABEL_GEOMETRY.json');
end
addpath(fullfile(root,'Code','Matlab','Functions'));
addpath(fullfile(root,'Code','Matlab','Functions','ModernRobotics'));
mat = fullfile(root,'Code','Matlab','Mesh_Optimization','Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat');
s = load(mat,'ctxCamBest','predBest');
ctx = s.ctxCamBest;
route = s.predBest.routeInfo;
angles = ctx.phiD(:);
changes = find(any(route.active(:,1:end-1) & ~route.active(:,2:end),1));
poses = unique([1,changes,numel(angles)],'stable');
if numel(poses) ~= 8
    targets = linspace(min(angles),max(angles),8);
    [~,poses] = min(abs(angles-targets),[],1);
end
assert(numel(unique(poses))==8);
panels = struct([]);
for k = 1:8
    idx = poses(k);
    points = route.raw(:,:,idx);
    for j = 6:9
        intermediate = RowVecTrans(ctx.T_ICR_t1(:,:,idx),points(j,:));
        points(j,:) = RowVecTrans(ctx.T_Pam(:,:,idx),intermediate);
    end
    panels(k).panel = char('A'+k-1);
    panels(k).frame = idx;
    panels(k).angle_deg = angles(idx);
    panels(k).point_ids = find(route.active(:,idx))';
    panels(k).femur_xy_m = points(route.active(:,idx),1:2);
end
out = struct('result','Vas_Pam_20mm_Result_cam_20261014_0416.mat', ...
    'panels',panels);
out.drum_radius_mm = 1000*route.wrap.radius(2,1);
released = find(isfinite(route.eliminatedAngleD));
out.released_rows = released(:)';
out.release_angles_deg = reshape(route.eliminatedAngleD(released),1,[]);
fid = fopen(outputFile,'w','n','UTF-8');
assert(fid >= 0,'Could not create the geometry file.');
guard = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(out,'PrettyPrint',true));
end
