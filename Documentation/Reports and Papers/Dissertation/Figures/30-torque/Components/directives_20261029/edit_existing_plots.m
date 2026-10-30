r='D:/Github/Bipedal_Robot';d=fullfile(r,'Documentation/Reports and Papers/Dissertation');o=fullfile(r,'tmp/chapter45_directives_20261029');
checkcode(fullfile(o,'edit_existing_plots.m'));
f=openfig(fullfile(d,'Figures/30-torque/Components/reconciliation_20261028/FlxGrp_current.fig'),'invisible');aa=findall(f,'Type','axes');a=[];b=[];
for j=1:numel(aa),if aa(j).Layout.Tile==1,a=aa(j);elseif aa(j).Layout.Tile==3,b=aa(j);end,end
assert(~isempty(a)&&~isempty(b));old=[a b];g=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 24 10]);
for j=1:2
 ax=copyobj(old(j),g);ax.Units='normalized';ax.Position=[.075+(j-1)*.50 .19 .40 .72];ax.FontSize=16;ax.Title.FontSize=18;
 xlabel(ax,'\theta_{\bf k}, ^\circ');ylabel(ax,'Torque, N m');text(ax,.02,.95,['(' char(64+j) ')'],'Units','normalized','FontWeight','bold','FontSize',18);
 sc=findall(ax,'Type','scatter');set(sc,'MarkerFaceAlpha',.75,'MarkerEdgeColor','none');
 before=findall(old(j),'-property','XData');after=findall(ax,'-property','XData');assert(numel(before)==numel(after));for k=1:numel(before),assert(isequaln(before(k).XData,after(k).XData)&&isequaln(before(k).YData,after(k).YData));end
 if j==2
 hold(ax,'on');hh=gobjects(5,1);hh(1)=plot(ax,nan,nan,'--','Color',[.980 .529 .459],'LineWidth',2);hh(2)=plot(ax,nan,nan,'-','Color',[.804 .204 .710],'LineWidth',2.5);hh(3)=scatter(ax,nan,nan,40,'b','filled','MarkerFaceAlpha',.75,'MarkerEdgeColor','none');hh(4)=plot(ax,nan,nan,':k','LineWidth',3);hh(5)=scatter(ax,nan,nan,40,[.918 .373 .580],'filled','MarkerFaceAlpha',.75,'MarkerEdgeColor','none');
 legend(ax,hh,{'Unoptimized','Optimized','Measured','Human model (max)','Measured at reduced pressure'},'Location','southwest','FontSize',15,'Box','off');
 end
end
exportgraphics(g,fullfile(d,'Figures/30-torque/FlxGrp.pdf'),'ContentType','vector');exportgraphics(g,fullfile(o,'FlxGrp_AB.png'),'Resolution',180);savefig(g,fullfile(o,'FlxGrp_AB.fig'));close([f g]);
f=openfig(fullfile(d,'Figures/30-torque/Components/publication_20261028/Redesigned_flexor_pressure_publication.fig'),'invisible');ax=findall(f,'Type','axes');sc=findall(ax,'Type','scatter');before=arrayfun(@(x)struct('x',x.XData,'y',x.YData),sc);
set(sc,'MarkerFaceAlpha',.75,'MarkerEdgeColor','none');m=[];
for j=1:numel(sc)
 if numel(sc(j).XData)==15 && ~contains(string(sc(j).DisplayName),'Prediction'),m=sc(j);end
end
assert(~isempty(m),'Measured series not found');parent=m;while ~isequal(parent.Parent,ax),parent=parent.Parent;end
uistack(parent,'top');
for j=1:numel(sc),assert(isequaln(before(j).x,sc(j).XData)&&isequaln(before(j).y,sc(j).YData));end
exportgraphics(f,fullfile(d,'Figures/30-torque/Redesigned_flexor_pressure_publication.pdf'),'ContentType','vector');exportgraphics(f,fullfile(o,'flexor_measured_last.png'),'Resolution',180);savefig(f,fullfile(o,'Redesigned_flexor_measured_last.fig'));close(f);
fid=fopen(fullfile(o,'plot_verification.txt'),'w');fprintf(fid,'Figure 4.11 retains only source tiles 1 and 3. All retained XData and YData bit-identical. Legend southwest in panel B. Figure 4.14 measured data group stacked last/top. Every scatter MarkerFaceAlpha=0.75 and MarkerEdgeColor=none. No model evaluation or optimization.\n');fclose(fid);
