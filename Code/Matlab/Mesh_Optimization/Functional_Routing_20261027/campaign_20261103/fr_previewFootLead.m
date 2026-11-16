function fr_previewFootLead(outputDir, chosen, branch)
% Plot a saved route against the same tibia field used for lead clearance.
C = fr_setup();
kin = fr_kinematics(C.osimStock);
file = fullfile(C.camp, 'audit_20261105', 'reoptimization', ...
    'tib_post', 'coordinated_compliant.mat');
if nargin < 2
    S = load(file, 'chosen');
    chosen = S.chosen;
end
if nargin < 3; branch = 'compliant'; end
cs = chosen.caseDef;
cs.allowFootLeadCorrection = true;
design = chosen.design;
[P, frames, ids] = fr_locForPass(C, cs, design, kin, 1);
poses = cs.passes{1}.qStructs;
owning = cs.rowOwning(ids);
q = poses{ceil(numel(poses)/2)};
lead = fr_tibialLeadMargin(C, kin, P, frames, owning, poses, 0);
fprintf('Tibialis posterior nominal tibial lead margin: %.3f mm\n', lead*1000);

field = load(fullfile(C.camp, 'audit_20261105', 'surface_fields', 'tibia_r.mat'));
[X,Y,Z] = meshgrid(field.x, field.y, field.z);
bone = isosurface(X,Y,Z,permute(double(field.distance), [2 1 3]), 0);
world = zeros(size(P));
for j = 1:size(P,1)
    T = fr_chainT(kin, q, frames{j}, 'tibia_r');
    world(j,:) = (T(1:3,1:3)*P(j,:)' + T(1:3,4))';
end
f = figure('Visible', 'off', 'Color', 'w');
cleanup = onCleanup(@() close(f));
ax = axes(f);
patch(ax, 'Faces', bone.faces, 'Vertices', 1000*bone.vertices, ...
    'FaceColor', [.68 .68 .68], 'EdgeColor', 'none', 'FaceAlpha', .5);
hold(ax, 'on');
routeLine = plot3(ax, 1000*world(:,1), 1000*world(:,2), 1000*world(:,3), ...
    '-o', 'Color', [0 .35 .8], 'LineWidth', 2, 'MarkerFaceColor', [0 .35 .8]);
poseIdx = ceil(numel(poses)/2);
env = fr_envelope(C, cs, design, kin, poses, branch);
T = fr_bodyT(kin, q, 'tibia_r');
sleeveLine = gobjects(0);
for sleeve = env.capsules{poseIdx}
    ends = 1000*(T(1:3,1:3)'*([sleeve.p1;sleeve.p2]'-T(1:3,4)))';
    sleeveLine = plot3(ax, ends(:,1), ends(:,2), ends(:,3), '-', ...
        'Color',[.9 .4 .12], 'LineWidth',6);
end
hardware = env.spheres{poseIdx};
stop = strcmp(hardware.what,'shank_tendon_stop');
stopMarker = gobjects(0);
if any(stop)
    at = 1000*(T(1:3,1:3)'*(hardware.p(stop,:)'-T(1:3,4)))';
    stopMarker = plot3(ax,at(:,1),at(:,2),at(:,3),'ks','MarkerFaceColor','k','MarkerSize',6);
end
legend(ax,[routeLine sleeveLine stopMarker], ...
    [{'Cable route'},repmat({'Solved sleeve axis'},size(sleeveLine)), ...
    repmat({'Shank stop'},size(stopMarker))],'Location','northeast');
axis(ax, 'equal'); grid(ax, 'on'); view(ax, 38, 17);
xlabel(ax, 'Tibia x (mm)'); ylabel(ax, 'Tibia y (mm)'); zlabel(ax, 'Tibia z (mm)');
title(ax, sprintf('Tibialis posterior | %d x %s mm BPA | %gx travel', ...
    cs.seed.nBPA, cs.diaStr, max(cs.txConfig.gain)));
if ~isfolder(outputDir)
    mkdir(outputDir);
end
exportgraphics(ax, fullfile(outputDir, 'tibialis_posterior_lead.png'), 'Resolution', 180);
% Export the evaluated points in their native owning frames for OpenSim.
native = zeros(size(P));
for j = 1:size(P,1)
    T = fr_chainT(kin,q,frames{j},owning{j});
    native(j,:) = (T(1:3,1:3)*P(j,:)' + T(1:3,4))';
end
record = struct('points_m',native,'bodies',{owning},'reference',C.osimRobot, ...
    'gain',max(cs.txConfig.gain),'branch',branch);
fid = fopen(fullfile(outputDir,'tibialis_posterior_route.json'),'w');
assert(fid>=0);
fileCleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(record));
end
