function fr_screenRectusCamSeed(outputDir)
% Screen the existing robot-ICR extensor cam as a rectus knee module.
% This compares leverage and stroke; it does not solve the hip connection.
C = fr_setup();
source = fullfile(C.repo,'Code','Matlab','Mesh_Optimization','Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat');
S = load(source,'predBest','ctxCamBest','xBestCam');
poses = readtable(fullfile(outputDir,'rectus_pose_grid.csv'));
knee = unique(poses.knee_deg);
hip = unique(poses.hip_deg);
ref = poses.hip_deg==hip(find(abs(hip)==min(abs(hip)),1));
human = poses.human_knee_arm_m(ref);
assert(isequal(poses.knee_deg(ref),knee),'Unexpected rectus grid ordering');
angles = S.ctxCamBest.phiD(:);
assert(min(knee)>=min(angles)-1e-5 && max(knee)<=max(angles)+1e-5);
cam = interp1(angles,S.predBest.momentArm(:), ...
    min(max(knee,min(angles)),max(angles)),'pchip');
T = table(knee,1000*human,2000*human,1000*cam,cam./human, ...
    'VariableNames',{'knee_deg','human_rectus_arm_mm','twice_human_mm', ...
    'cam_arm_about_robot_ICR_mm','cam_to_human_ratio'});
writetable(T,fullfile(outputDir,'rectus_vasti_cam_seed.csv'));
targetTravel = zeros(size(hip));
for j = 1:numel(hip)
    h = poses.hip_deg==hip(j);
    targetTravel(j) = trapz(deg2rad(knee),2*poses.human_knee_arm_m(h));
end
camTravelEquivalent = trapz(deg2rad(knee),cam);
result = struct('source',source,'saved_cam_variables',S.xBestCam, ...
    'reference_hip_deg',poses.hip_deg(find(ref,1)), ...
    'two_x_human_knee_travel_range_mm',1000*[min(targetTravel) max(targetTravel)], ...
    'cam_knee_travel_equivalent_mm',1000*camTravelEquivalent, ...
    'cam_reference_ratio_range',[min(T.cam_to_human_ratio) max(T.cam_to_human_ratio)], ...
    'reference_poses_meeting_twice_human',nnz(cam>=2*human),'reference_pose_count',numel(knee), ...
    'status','Seed only. Rectus proximal hip span, contact, equilibrium and clearance unverified.');
fid = fopen(fullfile(outputDir,'rectus_vasti_cam_seed.json'),'w');
assert(fid>=0);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(result));
disp(result);
f = figure('Visible','off','Color','w');
ax = axes(f);
plot(ax,knee,T.human_rectus_arm_mm,'--',knee,T.twice_human_mm,':', ...
    knee,T.cam_arm_about_robot_ICR_mm,'-','LineWidth',1.7);
grid(ax,'on');
xlabel(ax,'Knee angle (deg) · flexion negative');
ylabel(ax,'Extension moment arm (mm)');
title(ax,sprintf('Rectus femoris leverage | hip at %.1f deg',result.reference_hip_deg));
legend(ax,{'Stock human rectus femoris','Twice human','Existing robot extensor cam'}, ...
    'Location','best');
exportgraphics(ax,fullfile(outputDir,'rectus_vasti_cam_seed.png'),'Resolution',180);
close(f);
end

