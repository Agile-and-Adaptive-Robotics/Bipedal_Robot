% Change only the azimuth of the two saved comparison surfaces.
disp(version);
root=fileparts(mfilename('fullpath'));
asset=fullfile(root,'assets','mechanics','historical_bflh_common_scales');
audit=fullfile(root,'build','plot_rotation_180');
if ~isfolder(audit),mkdir(audit);end
for ext={'.fig','.png','.pdf'}
    prior=fullfile(audit,['before',ext{1}]);
    if ~isfile(prior),copyfile([asset,ext{1}],prior);end
end
f=openfig(fullfile(audit,'before.fig'),'invisible');
axesHandles=findobj(f,'Type','axes');
assert(numel(axesHandles)==2);
records=struct([]);
for i=1:numel(axesHandles)
    ax=axesHandles(i);s=findobj(ax,'Type','surface');assert(numel(s)==1);
    before={s.XData,s.YData,s.ZData,s.CData,ax.XLim,ax.YLim,ax.ZLim,ax.CLim,ax.Colormap,ax.XDir,ax.YDir,ax.ZDir};
    [az,el]=view(ax);assert(abs(az+37.5)<1e-9 && abs(el-30)<1e-9);
    view(ax,[az+180,el]);
    after={s.XData,s.YData,s.ZData,s.CData,ax.XLim,ax.YLim,ax.ZLim,ax.CLim,ax.Colormap,ax.XDir,ax.YDir,ax.ZDir};
    assert(isequaln(before,after),'A data, scale, color, or axis-direction property changed');
    records(i).before=[az,el];records(i).after=[az+180,el];
    records(i).data_scales_colors_unchanged=true;
end
drawnow;
exportgraphics(f,[asset,'.png'],'Resolution',180);
exportgraphics(f,[asset,'.pdf'],'ContentType','vector');
savefig(f,[asset,'.fig']);close(f);
fid=fopen(fullfile(audit,'matlab_check.json'),'w');
fprintf(fid,'%s',jsonencode(records));fclose(fid);
fprintf('Two surface views rotated 180 degrees. Data, scales, and colors unchanged.\n');
