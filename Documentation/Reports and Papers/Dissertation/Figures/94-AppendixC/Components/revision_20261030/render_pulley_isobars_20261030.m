function render_pulley_isobars_20261030()
root='D:/Github/Bipedal_Robot';F=fullfile(root,'Documentation/Reports and Papers/Dissertation/Figures/94-AppendixC');
addpath(fullfile(root,'Code/Matlab'));Colors;palette=d;
S=load(fullfile(root,'Code/Matlab/Mesh_Optimization/Results/PulleyDoctrine_figdata_20261021.mat'),'allData');
names={'Soleus','Tibialis anterior','BFSH','Medial gastrocnemius','Rectus femoris'};
fig=figure('Color','w','Units','inches','Position',[1 1 8.4 4.8]);t=tiledlayout(2,3,'TileSpacing','loose','Padding','loose');
ci=[7 2 6 4 5];styles={'-',':','--','-.','-'};marks={'o','s','^','d','v'};pressureLabels={};
for k=1:5
 ax=nexttile(t,k);dd=S.allData(k).T_2BPA_G2;theta=dd.th1m(:);p=dd.pressSurf(:);Z=abs(dd.tauSurfP);req=abs(dd.tauReqM(:));
 assert(size(Z,1)==numel(theta)&&size(Z,2)==numel(p),'Saved torque grid dimensions do not match axes.');
 ix=unique(round(linspace(1,numel(p),5)));assert(numel(ix)==5,'Five saved pressure levels required.');
 hold(ax,'on');hh=gobjects(6,1);pl=cell(1,5);
 for j=1:5
  pl{j}=sprintf('%g kPa',p(ix(j)));
  hh(j)=plot(ax,theta,Z(:,ix(j)),'Color',palette(ci(j),:),'LineStyle',styles{j},'Marker',marks{j},'MarkerIndices',1:max(1,round(numel(theta)/6)):numel(theta),'MarkerSize',3.4,'LineWidth',1.4,'MarkerFaceColor',palette(ci(j),:),'MarkerEdgeColor','none');
 end
 hh(6)=plot(ax,theta,req,'k--','LineWidth',2.0);
 if k==1,pressureLabels=pl;legendHandles=hh;else,assert(isequal(pressureLabels,pl),'Pressure grids differ between panels.');end
 title(ax,sprintf('(%c) %s',64+k,names{k}),'FontSize',11.7,'FontWeight','bold');
 xlabel(ax,'Joint angle (deg)','FontSize',11);ylabel(ax,'Torque magnitude (N m)','FontSize',11);
 set(ax,'FontName','Arial','FontSize',10.5,'LineWidth',1,'Box','off','XMinorTick','off','YMinorTick','off','TickDir','out');
 xlim(ax,[min(theta) max(theta)]);ylim(ax,[0 1.08*max([Z(:);req(:)])]);grid(ax,'off');ax.Title.Units='normalized';ax.Title.Position(2)=1.04;
 fprintf('%s: saved pressure indices %s; kPa %s\n',names{k},mat2str(ix),mat2str(p(ix)'));
end
lg=legend(legendHandles,[pressureLabels,{'Human reference'}],'FontName','Arial','FontSize',10.5,'Box','off','NumColumns',1);lg.Layout.Tile=6;
exportgraphics(fig,fullfile(F,'pulley_torque_isobars_20261030.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(F,'pulley_torque_isobars_20261030.png'),'Resolution',220);close(fig);
end
