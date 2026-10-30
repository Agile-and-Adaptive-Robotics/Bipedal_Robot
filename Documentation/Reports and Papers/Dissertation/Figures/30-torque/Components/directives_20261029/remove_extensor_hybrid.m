root='D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation';
src=fullfile(root,'Figures/30-torque/Components/publication_20261028/Ext10mm_52cm_current_publication.fig');
out=fullfile(root,'Figures/30-torque/Components/directives_20261029');
f=openfig(src,'invisible');
axesAll=findall(f,'Type','axes');
ax=axesAll(arrayfun(@(a) any(arrayfun(@(s) strcmp(s.DisplayName,'Approximate hybrid'),findall(a,'Type','scatter'))),axesAll));
assert(numel(ax)==1);
hybrid=findall(ax,'Type','scatter','DisplayName','Approximate hybrid');
assert(numel(hybrid)==1);
keepers=findall(ax,'-property','XData');keepers=keepers(keepers~=hybrid);
before=arrayfun(@(h) struct('x',h.XData,'y',h.YData),keepers);
delete(hybrid);
lg=findall(f,'Type','legend');delete(lg);
objects=findall(ax,'-property','DisplayName');
names={'Measured','Improved BPA','Optimized model'};handles=gobjects(1,3);
for k=1:3
    match=objects(arrayfun(@(h) strcmp(h.DisplayName,names{k}),objects));assert(numel(match)==1);handles(k)=match;
end
lg=legend(ax,handles,names,'Orientation','horizontal','NumColumns',3,'FontName','Arial','FontSize',16);
lg.Layout.Tile='south';
for h=findall(ax,'Type','scatter')'
    h.MarkerFaceAlpha=.75;h.MarkerEdgeColor='none';
end
set(ax,'FontName','Arial','FontSize',16);
set(findall(f,'Type','textboxshape'),'FontSize',16);
after=arrayfun(@(h) struct('x',h.XData,'y',h.YData),keepers);
assert(isequaln(before,after),'Retained plot data changed');
savefig(f,fullfile(out,'Ext10mm_52cm_no_hybrid.fig'));
exportgraphics(f,fullfile(root,'Figures/30-torque/Ext10mm_52cm_current_publication.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(out,'Ext10mm_52cm_no_hybrid_preview.png'),'Resolution',160);
fid=fopen(fullfile(out,'extensor_hybrid_removal_verification.txt'),'w');
fprintf(fid,'Removed one Approximate hybrid scatter; retained all other XData/YData exactly.\nPhoto panel A retained. Legend outside axes in south layout tile. Scatter alpha 0.75; edges none.\n');fclose(fid);
close(f);disp('HYBRID_REMOVAL_VERIFIED');
