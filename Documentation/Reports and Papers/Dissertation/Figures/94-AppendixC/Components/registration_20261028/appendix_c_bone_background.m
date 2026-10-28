function appendix_c_bone_background(ax,ctx,ii)
% Presentation-only registration; saved route and clearance data are untouched.
persistent fem tib
if isempty(fem)
 mesh='D:/Github/Bipedal_Robot/Code/Matlab/Bone_Mesh_Plots/Open_Sim_Bone_Geometry';
 fem=readmatrix(fullfile(mesh,'Femur_Mesh_Points.txt'));tib=readmatrix(fullfile(mesh,'Tibia_Mesh_Points.txt'));
 fem=fem(all(isfinite(fem),2),:);tib=tib(all(isfinite(tib),2),:);
end
xl=xlim(ax);yl=ylim(ax);
% Condyle center and proximal tibial groove registered to the robot envelopes.
ff=fem;ff(:,1:2)=ff(:,1:2)+ctx.geo.femurProfileCenter-[0.00652,-0.40800];
tt=tib;tt(:,1:2)=tt(:,1:2)+ctx.geo.tibiaLowerCenter-[0.00600,-0.031255];
% Only proximal tibia/shaft is shown; a transverse crop is not a foot.
tt=tt(tt(:,2)>=-0.250,:);
T=ctx.T_Pam(:,:,ii)*ctx.T_ICR_t1(:,:,ii);tp=tt*T(1:3,1:3)'+T(1:3,4)';
clouds={ff,tp};
for j=1:2
 q=unique(clouds{j}(:,1:2),'rows');k=boundary(q(:,1),q(:,2),.75);
 h=patch(ax,q(k,1),q(k,2),[.89 .89 .89],'EdgeColor',[.65 .65 .65],'LineWidth',.6,'HandleVisibility','off');uistack(h,'bottom');
end
xlim(ax,xl);ylim(ax,yl);
end
