folder='D:/Github/Bipedal_Robot/tmp/reconciliation_figures_20261028';
fig=openfig(fullfile(folder,'FlxGrp_current.fig'),'invisible');
a=findall(fig,'Type','textboxshape');
for j=1:numel(a)
 if contains(a(j).String,'(C)')
  pos=a(j).Position;pos(1)=pos(1)-0.035;a(j).Position=pos;
 end
end
leg=findall(fig,'Type','legend');set(leg,'FontSize',12);
savefig(fig,fullfile(folder,'FlxGrp_current.fig'));
exportgraphics(fig,fullfile(folder,'FlxGrp_20261023.eps'),'ContentType','vector');
exportgraphics(fig,fullfile(folder,'FlxGrp_2x2_20261026.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(folder,'FlxGrp_2x2_20261026.png'),'Resolution',300);
close(fig);disp('LAYOUT_POLISHED_NO_DATA_CHANGE');
