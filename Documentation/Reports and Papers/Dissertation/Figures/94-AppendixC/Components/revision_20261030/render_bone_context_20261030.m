function render_bone_context_20261030()
root='D:/Github/Bipedal_Robot';
addpath(fullfile(root,'Code/Matlab'));
addpath(fullfile(root,'Code/Matlab/Robot_Data'));
addpath(fullfile(root,'Code/Matlab/Human_Data'));
addpath(fullfile(root,'Code/Matlab/Functions'));
addpath(fullfile(root,'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root,'Code/Matlab/Bone_Mesh_Plots'));
addpath(fullfile(root,'Code/Matlab/Bone_Mesh_Plots/Open_Sim_Bone_Geometry'));
S=load(fullfile(root,'Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result_20260925.mat'),'ctx','predBest');
ctx=S.ctx;predBest=S.predBest;phiD=ctx.phiD;phi=ctx.phi;positions=ctx.N;T=ctx.T;

cdeg = pi/180;
rect_fem_x = [0.0156367;0.0179948;0.0274274;0.029683;0.0306;0.0366;0.0422;0.0451;0.0484;0.0533;0.0617;0.0634;0.067;0.0733];
rect_fem_xD = cdeg*[-120.118;-114.871;-90.068;-83.532;-80;-60;-40;-30;-20;-10;0;1.6;5;10];
rect_fem_y = [0.0234;0.0238;0.0251;0.0253;0.025284;0.0249;0.0243;0.0239;0.0234;0.0228;0.0210;0.0206;0.0192;0.0160];
rect_fem_yD = cdeg*[-120;-114.6;-90;-83.5;-80.01;-60;-40;-30;-20;-10;0;1.6;5;10];

fcn5 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.009811),'smoothingspline');
fcn6 = fit(rect_fem_yD,rect_fem_y-(0.02346-0.02242),'smoothingspline');
fcn7 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.008285),'smoothingspline');
fcn8 = fit(rect_fem_yD,rect_fem_y+(0.0256239-0.02346),'smoothingspline');
fcn9 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.0142881),'smoothingspline');
fcn10 = fit(rect_fem_yD,rect_fem_y-(0.02346-0.0215281),'smoothingspline');

Name = 'Vastus Medialis';
MIF = 1294;
OFL = 0.089;
TSL = 0.126;
Pennation = 0.08726646;
HumanLocation = zeros(5,3,positions);

for i = 1:positions
    if phiD(i) < -101
        HumanLocation(:,:,i) = [0.014,-0.210,0.019; ...
                                0.036,-0.277,0.001; ...
                                0.037,-0.405,-0.013; ...
                                0.027,-0.425,-0.013; ...
                                fcn5(phi(i)),fcn6(phi(i)),-0.0146];
    elseif phiD(i) < -69 && phiD(i) >= -101
        HumanLocation(:,:,i) = [0.014,-0.210,0.019; ...
                                0.036,-0.277,0.001; ...
                                0.037,-0.405,-0.013; ...
                                0.037,-0.405,-0.013; ...
                                fcn5(phi(i)),fcn6(phi(i)),-0.0146];
    else
        HumanLocation(:,:,i) = [0.014,-0.210,0.019; ...
                                0.036,-0.277,0.001; ...
                                0.036,-0.277,0.001; ...
                                0.036,-0.277,0.001; ...
                                fcn5(phi(i)),fcn6(phi(i)),-0.0146];
    end
end

Vas_Med = MonoMuscleData(Name,HumanLocation,5,MIF,TSL,Pennation,OFL,T);

Name = 'Vastus Intermedius';
MIF = 1365;
OFL = 0.087;
TSL = 0.136;
Pennation = 0.05235988;
HumanLocation = zeros(4,3,positions);

for i = 1:positions
    if phiD(i) < -80
        HumanLocation(:,:,i) = [0.029,-0.192,0.031; ...
                                0.034,-0.208,0.029; ...
                                0.034,-0.403,0.005; ...
                                fcn7(phi(i)),fcn8(phi(i)),0.0018];
    else
        HumanLocation(:,:,i) = [0.029,-0.192,0.031; ...
                                0.034,-0.208,0.029; ...
                                0.034,-0.208,0.029; ...
                                fcn7(phi(i)),fcn8(phi(i)),0.0018];
    end
end

Vas_Int = MonoMuscleData(Name,HumanLocation,4,MIF,TSL,Pennation,OFL,T);

Name = 'Vastus Lateralis';
MIF = 1871;
OFL = 0.084;
TSL = 0.157;
Pennation = 0.08726646;
HumanLocation = zeros(5,3,positions);

for i = 1:positions
    if phiD(i) < -110
        HumanLocation(:,:,i) = [0.005,-0.185,0.035; ...
                                0.027,-0.259,0.041; ...
                                0.036,-0.403,0.021; ...
                                0.025,-0.424,0.018; ...
                                fcn9(phi(i)),fcn10(phi(i)),0.0165];
    elseif phiD(i) < -69 && phiD(i) >= -110
        HumanLocation(:,:,i) = [0.005,-0.185,0.035; ...
                                0.027,-0.259,0.041; ...
                                0.036,-0.403,0.021; ...
                                0.036,-0.403,0.021; ...
                                fcn9(phi(i)),fcn10(phi(i)),0.0165];
    else
        HumanLocation(:,:,i) = [0.005,-0.185,0.035; ...
                                0.027,-0.259,0.041; ...
                                0.027,-0.259,0.041; ...
                                0.027,-0.259,0.041; ...
                                fcn9(phi(i)),fcn10(phi(i)),0.0165];
    end
end

Vas_Lat = MonoMuscleData(Name,HumanLocation,5,MIF,TSL,Pennation,OFL,T);


AnimateKneeBoneMuscle(ctx.T,ctx.T_ICR_t1,ctx.phi,ctx.pos,predBest.p1,predBest.pEnd, ...
 {Vas_Int,Vas_Lat,Vas_Med},'Static',true,'StaticFullSkeleton',true, ...
 'DisplayRotation',eye(3),'DisplayAxisMap',eye(3), ...
 'Location',predBest.Location,'CrossPoint',ctx.CrossPoint,'T_Pam',ctx.T_Pam, ...
 'HumanLabels',{'Vastus intermedius','Vastus lateralis','Vastus medialis'}, ...
 'XLim',[],'YLim',[],'ZLim',[],'StaticPadding',0.01);
fig=gcf;ax=gca; Colors; palette=d;
set(fig,'Units','inches','Position',[1 1 6.5 6.1],'Color','w');
lg=findobj(fig,'Type','legend');legHandles=lg.PlotChildren;legLabels=cellstr(lg.String);delete(lg);
tags={'FemurMesh','TibiaMesh','PelvisMesh','SacrumMesh','SpineMesh','FootMesh'};
labels={'Femur','Tibia','Pelvis / hip','Sacrum','Lumbar spine','Foot'};ci=[2 3 6 4 5 7];
for j=1:6
 h=findall(ax,'Tag',tags{j});set(h,'Color',palette(ci(j),:),'MarkerFaceColor',palette(ci(j),:),'MarkerEdgeColor','none');
 li=find(strcmpi(legLabels,labels{j}));if ~isempty(li),set(legHandles(li),'MarkerFaceColor',palette(ci(j),:),'Color',palette(ci(j),:));end
end
set(findobj(ax,'Tag','RobotRoute'),'Color',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',2.3);
humanNames={'Vastus intermedius','Vastus lateralis','Vastus medialis'};
for j=1:3
 h=findobj(ax,'Tag',sprintf('HumanRoute%d',j));styles={'--',':','-.'};set(h,'Color',[.36 .36 .36],'LineStyle',styles{j},'LineWidth',1.6,'Marker','none');
 li=find(contains(legLabels,humanNames{j}));if ~isempty(li),set(legHandles(li),'Color',[.36 .36 .36],'LineStyle',styles{j},'Marker','none');end
end
set(ax,'FontName','Arial','FontSize',12,'FontWeight','normal','LineWidth',1,'XMinorTick','off','YMinorTick','off','ZMinorTick','off');
clouds=findall(ax,'-property','XData');XYZ=[];
for j=1:numel(clouds)
 h=clouds(j);if isprop(h,'ZData'),q=[h.XData(:),h.YData(:),h.ZData(:)];q=q(all(isfinite(q),2),:);XYZ=[XYZ;q];end
end
lims=[min(XYZ);max(XYZ)];pad=max(.01,.025*diff(lims));xlim(ax,lims(:,1)'+[-pad(1),pad(1)]);ylim(ax,lims(:,2)'+[-pad(2),pad(2)]);zlim(ax,lims(:,3)'+[-pad(3),pad(3)]);
set(ax,'Position',[.02 .31 .45 .56]);title(ax,'(A) Lower-body context','FontSize',13,'FontWeight','bold');ax.Title.Units='normalized';ax.Title.Position(2)=1.10;
ax.CameraUpVector=[0 1 0];view(ax,145,12);ax.CameraUpVector=[0 1 0];
xlabel(ax,'x (m)');ylabel(ax,'y (m)');zlabel(ax,'z (m)');set(ax,'XTick',[],'ZTick',[]);
ax2=copyobj(ax,fig);set(ax2,'Position',[.55 .31 .41 .56]);view(ax2,0,90);ax2.CameraUpVector=[0 1 0];
xlim(ax2,[-.17 .08]);ylim(ax2,[-.64 -.25]);set(ax2,'XTick',[-.15 -.05 .05]);title(ax2,'(B) Knee close-up','FontSize',13,'FontWeight','bold');ax2.Title.Units='normalized';ax2.Title.Position(2)=1.10;
zlabel(ax2,'');xlabel(ax2,'x (m)');ylabel(ax2,'y (m)');
lg=legend(ax,legHandles,legLabels,'NumColumns',3,'FontName','Arial','FontSize',10.5,'Box','off','AutoUpdate','off');lg.Units='normalized';lg.Position=[.03 .055 .94 .19];
F=fullfile(root,'Documentation/Reports and Papers/Dissertation/Figures/94-AppendixC');
% Point-cloud image is exported as a high-resolution raster; route arrays are unchanged.
set(fig,'PaperUnits','inches','PaperPosition',[0 0 6.5 6.1],'PaperSize',[6.5 6.1]);
print(fig,fullfile(F,'bone_context_extensor20mm_20261030.png'),'-dpng','-r300');

close(fig);
end
