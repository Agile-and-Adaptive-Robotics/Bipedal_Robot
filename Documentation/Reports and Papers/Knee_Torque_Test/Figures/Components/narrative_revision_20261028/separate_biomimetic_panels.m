root='D:/Github/Bipedal_Robot';
out=fullfile(root,'Documentation/Reports and Papers/Dissertation/Archive/reconciliation_20261028/narrative_revision');
src=openfig(fullfile(root,'tmp/reconciliation_figures_20261028/FlxGrp_current.fig'),'invisible');
aa=findall(src,'Type','axes');
for j=1:numel(aa)
    tile=aa(j).Layout.Tile;
    if tile==1,a=aa(j);elseif tile==3,b=aa(j);elseif tile==4,c=aa(j);end
end
% Copy existing plotted objects; no torque or measurement calculation is run.
f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 19 10]);
for j=1:2
    old=[a b]; ax=copyobj(old(j),f);ax.Units='normalized';
    ax.Position=[.10+(j-1)*.48,.25,.38,.65];
    ax.FontSize=11;ax.Title.FontSize=12;
    xlabel(ax,'\theta_{\bf k}, ^\circ');ylabel(ax,'M_{\bf k}, N\cdotm');
    text(ax,-.19,1.04,['(' char(64+j) ')'],'Units','normalized','FontWeight','bold','FontSize',12);
end
la=axes(f,'Position',[.1 .01 .8 .12],'Visible','off');hold(la,'on');
plot(la,nan,nan,'--','Color',[.980 .529 .459],'LineWidth',2,'DisplayName','Improved BPA');
plot(la,nan,nan,'-','Color',[.804 .204 .710],'LineWidth',2.5,'DisplayName','Optimized model');
scatter(la,nan,nan,40,'b','filled','DisplayName','Measured');
plot(la,nan,nan,':k','LineWidth',3,'DisplayName','Human model (max)');
scatter(la,nan,nan,40,[.918 .373 .580],'filled','DisplayName','Reduced pressure');
legend(la,'Location','south','Orientation','horizontal','NumColumns',3,'FontSize',10,'Box','off');
exportgraphics(f,fullfile(out,'Biomimetic_flexor_AB.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(out,'Biomimetic_flexor_AB.png'),'Resolution',170);savefig(f,fullfile(out,'Biomimetic_flexor_AB.fig'));
g=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 16 10]);
ax=copyobj(c,g);ax.Units='normalized';ax.Position=[.14 .18 .82 .75];ax.YTickLabelMode='auto';ax.FontSize=12;ax.Title.String='';
xlabel(ax,'\theta_{\bf k}, ^\circ');ylabel(ax,'M_{\bf k}, N\cdotm');
hold(ax,'on');
plot(ax,nan,nan,':k','LineWidth',3,'DisplayName','Human model (BFSH)');
plot(ax,nan,nan,'-','Color',[.804 .204 .710],'LineWidth',2.5,'DisplayName','Optimized model, 620 kPa');
scatter(ax,nan,nan,40,'b','filled','DisplayName','Measured');
legend(ax,'Location','southwest','FontSize',10,'Box','off');
exportgraphics(g,fullfile(out,'Redesigned_flexor_validation.pdf'),'ContentType','vector');
exportgraphics(g,fullfile(out,'Redesigned_flexor_validation.png'),'Resolution',170);savefig(g,fullfile(out,'Redesigned_flexor_validation.fig'));
close([src f g]);disp('EXISTING_SERIES_COPIED_NO_RECALCULATION');
