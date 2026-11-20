function export_extensor_data(outputFile)
% Collect the six adopted pinned tests and the two cam endpoint poses.
sourceDir = fileparts(mfilename('fullpath'));
if nargin == 0
    outputFile = fullfile(sourceDir,'extensor_data.mat');
end
repo = fileparts(sourceDir);
while ~isfolder(fullfile(repo,'Testing_Data','2022_02_Festo'))
    higher = fileparts(repo);
    assert(~strcmp(higher,repo),'Cannot find the repository root.');
    repo = higher;
end
startDir = pwd;
back = onCleanup(@() cd(startDir)); %#ok<NASGU>
addpath(fullfile(repo,'Code','Matlab','Functions'), ...
    fullfile(repo,'Code','Matlab','Functions','ModernRobotics'), ...
    fullfile(repo,'Code','Matlab','Robot_Data'));
cd(fullfile(repo,'Testing_Data','2022_02_Festo'));
selection = load('minimizeExt10mmX3_results.mat','x_actual','results_sort_actual');
if isfield(selection,'x_actual')
    terms = selection.x_actual(32,1:4);
else
    terms = selection.results_sort_actual(32,4:7);
end
expected = [-6.378022e-3,3.997873e4,1.473375e4,1.582887e-1];
assert(all(abs(terms-expected) < 1e-6*abs(expected)), ...
    'The adopted pinned-extensor row has changed.');
[~,bpa] = minimizeExtX3(terms(1),terms(2),terms(3),terms(4),1:9,'2trans');
order = [5 6 7 1 2 8];
tests = struct([]);
for k=1:numel(order)
    trial = bpa(order(k));
    tests(k).id = order(k);
    tests(k).rest = trial.rest;
    tests(k).angle = trial.Ak(:);
    tests(k).hybrid_angle = trial.A_h(:);
    tests(k).hybrid = trial.M_h(:);
    tests(k).measured_angle = trial.Aexp(:);
    tests(k).measured = trial.Mexp(:);
    tests(k).unoptimized = trial.M(:);
    tests(k).optimized = trial.M_p(:,3);
end

record = load(fullfile(repo,'Code','Matlab','Mesh_Optimization','Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'),'ctxCamBest','predBest');
ctx = record.ctxCamBest;
route = record.predBest.routeInfo;
geometry = ctx.geo;
mesh = fullfile(repo,'Code','Matlab','Bone_Mesh_Plots','Open_Sim_Bone_Geometry');
femur = readmatrix(fullfile(mesh,'Femur_Mesh_Points.txt'));
tibia = readmatrix(fullfile(mesh,'Tibia_Mesh_Points.txt'));
femur = femur(all(isfinite(femur),2),:);
tibia = tibia(all(isfinite(tibia),2),:);
% Landmark registration from Appendix C's appendix_c_bone_background.m.
femur(:,1:2) = femur(:,1:2) + geometry.femurProfileCenter - [0.00652,-0.40800];
tibia(:,1:2) = tibia(:,1:2) + geometry.tibiaLowerCenter - [0.00600,-0.031255];
tibia = tibia(tibia(:,2)>=-0.250,:);
poses = struct([]);
frames = [1 numel(ctx.phiD)];
for k=1:2
    i = frames(k);
    transform = ctx.T_Pam(:,:,i)*ctx.T_ICR_t1(:,:,i);
    points = route.raw(:,:,i);
    points(6:9,:) = applyTransform(transform,points(6:9,:));
    poses(k).angle = ctx.phiD(i);
    poses(k).point = points;
    poses(k).active = route.active(:,i);
    poses(k).tibia_transform = transform;
    poses(k).femur = outline(femur);
    poses(k).tibia = outline(applyTransform(transform,tibia));
end
assert(all(abs([poses.angle]-[-120 10]) < 1e-6), ...
    'Cam endpoints differ from the published angles at their reporting precision.');
assert(isequal(find(poses(1).active).',1:9));
assert(isequal(find(poses(2).active).',[1 7 8 9]));
cd(startDir);
save(outputFile,'tests','terms','poses','geometry','-v7');
end

function out = applyTransform(transform,points)
out = points*transform(1:3,1:3)' + transform(1:3,4)';
end

function out = outline(points)
xy = unique(points(:,1:2),'rows');
out = xy(boundary(xy(:,1),xy(:,2),0.75),:);
end
