function info = Polish_KneeTorque_Figures(mode, fig, axesAB)
% Re-export retained plot series without evaluating a force or torque model.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
diss=fullfile(root,'Documentation','Reports and Papers','Dissertation');
paper=fullfile(root,'Documentation','Reports and Papers','Knee_Torque_Test');
review=fullfile(diss,'Notes','figure_polish_20261030');
if nargin<2
    if strcmp(mode,'AB')
        source=fullfile(diss,'Figures','30-torque','Components','directives_20261029','FlxGrp_AB.fig');
    else
        source=fullfile(diss,'Figures','30-torque','Components','approved_20261030','Paper_flexor_pressure.fig');
    end
    fig=openfig(source,'new','invisible');
    axesAB=findall(fig,'Type','axes');
    keep=arrayfun(@(a)strcmp(a.Visible,'on'),axesAB);axesAB=axesAB(keep);
    x=arrayfun(@(a)a.Position(1),axesAB);[~,ix]=sort(x);axesAB=axesAB(ix);
else
    source='Caller-supplied regenerated plot';
end
objects=findall(fig,'-property','XData');
objects=objects(arrayfun(@(h)~isempty(h.XData)&&any(isfinite(h.XData(:))),objects));
before=cell(size(objects));
for k=1:numel(objects)
    before{k}=struct('x',objects(k).XData,'y',objects(k).YData,'size',[]);
    if isprop(objects(k),'SizeData'),before{k}.size=objects(k).SizeData;end
end

if strcmp(mode,'AB')
    assert(numel(axesAB)==2,'A/B figure must have two data axes');
    axesAB(1).XLim=[-120 15];axesAB(2).XLim=[-125 10];
    oldtext=findall(fig,'Type','text');
    for k=1:numel(oldtext)
        if any(strcmp(strtrim(string(oldtext(k).String)),["(A)","(B)","(C)"]))
            delete(oldtext(k));
        end
    end
    delete(findall(fig,'Type','legend'));
    ax=axesAB(2); real=findall(ax,'-property','XData');
    real=real(arrayfun(@(h)~isempty(h.XData)&&any(isfinite(h.XData(:))),real));
    order=[byColor(real,'line',[.980 .529 .459]), ...
           byColor(real,'line',[.804 .204 .710]), ...
           byColor(real,'line',[0 0 0]), ...
           byColor(real,'scatter',[0 0 1]), ...
           byColor(real,'scatter',[.918 .373 .580])];
    lg=legend(ax,order,{'Unoptimized','Optimized','Human model (max)',...
        'Measured','Measured at reduced pressure'},'Location','southwest',...
        'FontSize',15,'FontName','Arial','Box','off','AutoUpdate','off');
    lg.Orientation='vertical';lg.NumColumns=1;
    fig.Units='points';
    letters=gobjects(1,2);
    for k=1:2
        ax=axesAB(k);ax.Units='points';
        ax.FontSize=16;ax.FontName='Arial';ax.FontWeight='bold';
        ax.XLabel.FontSize=16;ax.YLabel.FontSize=16;ax.Title.FontSize=18;
        ax.Title.Units='points';ax.Title.VerticalAlignment='bottom';
        ax.YLabel.Units='points';
        letters(k)=text(ax,0,0,['(' char(64+k) ')'],...
            'FontName',ax.Title.FontName,'FontSize',ax.Title.FontSize,...
            'FontWeight','bold','Clipping','off','Tag','TorquePanelLetter',...
            'VerticalAlignment','bottom','HorizontalAlignment','left');
        letters(k).Units='points';
    end
    drawnow;
    left=zeros(1,2);
    for k=1:2
        ax=axesAB(k); e=ax.YLabel.Extent;w=letters(k).Extent(3);
        letters(k).Position=[e(1)-8-w,ax.Title.Position(2),0];
        left(k)=ax.Position(1)+letters(k).Position(1);
    end
    padding=max(0,4-min(left));
    p=fig.Position;p(3)=p(3)+padding;fig.Position=p;
    for k=1:2
        p=axesAB(k).Position;p(1)=p(1)+padding;axesAB(k).Position=p;
    end
    drawnow;
    info.panelLetters=struct([]);
    for k=1:2
        ax=axesAB(k); e=ax.YLabel.Extent;
        le=letters(k).Extent;te=ax.Title.Extent;
        assert(le(1)+le(3)<e(1),'Panel letter must be left of y-axis label');
        assert(abs(letters(k).Position(2)-ax.Title.Position(2))<1e-10);
        row=struct('text',letters(k).String,'letterBounds',le+ [ax.Position(1:2),0,0],...
            'ylabelBounds',e+[ax.Position(1:2),0,0],'titleBounds',te+[ax.Position(1:2),0,0]);
        if k==1,info.panelLetters=row;else,info.panelLetters(k)=row;end
    end
    name='FlxGrp'; destinations={fullfile(diss,'Figures','30-torque','FlxGrp.pdf'),...
        fullfile(paper,'Figures','Biomimetic_flexor_AB_publication.pdf')};
else
    assert(strcmp(mode,'physical')&&isscalar(axesAB));
    ax=axesAB(1);ax.FontUnits='points';ax.FontSize=11;
    alltext=findall(fig,'Type','text');
    for k=1:numel(alltext)
        alltext(k).FontUnits='points';
        if alltext(k).FontSize>12,alltext(k).FontSize=12;end
    end
    ax.XLabel.FontSize=12;ax.YLabel.FontSize=12;ax.Title.FontSize=12;
    sc=findall(ax,'Type','scatter');
    set(sc,'MarkerFaceAlpha',.75,'MarkerEdgeColor','none');
    prediction=byColor(sc,'scatter',[1 .843 0]);
    measured=byColor(sc,'scatter',[0 0 1]);
    uistack(prediction,'bottom');uistack(measured,'top');
    uistack(plotParent(prediction,ax),'bottom');
    uistack(plotParent(measured,ax),'top');
    lg=findall(fig,'Type','legend');lg=lg(1);
    lg.FontUnits='points';lg.FontSize=10.5;lg.NumColumns=1;lg.Orientation='vertical';
    lg.Location='southwest';lg.AutoUpdate='off';
    % Retain the current paper's uniform curve weight alongside its extensor.
    lines=findall(ax,'Type','line');set(lines,'LineWidth',1.2);
    drawnow;
    ax.Units='points';lg.Units='points';p=lg.Position;
    p(1)=max(p(1),ax.Position(1)+14);lg.Position=p;
    drawnow;
    info.legendBounds=lg.Position;
    info.axesBounds=ax.Position;
    ordered=orderedSeries(ax);
    info.measuredDrawnLast=find(ordered==measured,1)<find(ordered==prediction,1);
    assert(info.measuredDrawnLast);
    name='Redesigned_flexor_pressure_publication';
    destinations={fullfile(diss,'Figures','30-torque',[name '.pdf']),...
        fullfile(paper,'Figures',[name '.pdf'])};
end

sc=findall(fig,'Type','scatter');
set(sc,'MarkerFaceAlpha',.75,'MarkerEdgeColor','none');
for k=1:numel(objects)
    assert(isequaln(before{k}.x,objects(k).XData)&&isequaln(before{k}.y,objects(k).YData),'Plot series changed');
    if isprop(objects(k),'SizeData'),assert(isequaln(before{k}.size,objects(k).SizeData),'Marker sizes changed');end
end
info.mode=mode;info.source=source;info.seriesCount=numel(objects);info.dataConserved=true;
info.legendStrings=lg.String;info.legendFont=lg.FontSize;
info.scatterAlpha=[sc.MarkerFaceAlpha];
info.scatterEdges={sc.MarkerEdgeColor};
drawnow; % Complete native legend/annotation geometry BEFORE vector export.
candidate=getenv('BPA_FIGURE_REVIEW_DIR');
if isempty(candidate)
    candidate=review;
end
if ~isfolder(candidate),mkdir(candidate);end
if strcmp(mode,'physical')
    % Font sizes compensate the existing manuscript inclusion scales. Native
    % 10.5 pt becomes about 10.14 pt at .97*468/470 in the dissertation;
    % 12 pt becomes 10.2 pt at the paper's unchanged .85 scale.
    profiles={'dissertation','paper'};legendSizes=[10.5 12];tickSizes=[11 12];
    rows=cell(1,2);
    for j=1:2
        ax.FontSize=tickSizes(j);lg.FontSize=legendSizes(j);
        ax.XLabel.FontSize=12;ax.YLabel.FontSize=12;ax.Title.FontSize=12;
        lg.Location='southwest';drawnow;
        p=lg.Position;p(1)=max(p(1),ax.Position(1)+14);lg.Position=p;
        drawnow;
        for k=1:numel(objects)
            assert(isequaln(before{k}.x,objects(k).XData)&&isequaln(before{k}.y,objects(k).YData),'Plot series changed');
            if isprop(objects(k),'SizeData'),assert(isequaln(before{k}.size,objects(k).SizeData),'Marker sizes changed');end
        end
        output=fullfile(candidate,'physical_profiles',profiles{j});
        if ~isfolder(output),mkdir(output);end
        pdf=fullfile(output,[name '.pdf']);
        exportgraphics(fig,pdf,'ContentType','vector');
        exportgraphics(fig,fullfile(output,[name '_preview.png']),'Resolution',180);
        savefig(fig,fullfile(output,[name '.fig']));
        row=info;row.profile=profiles{j};row.legendFont=lg.FontSize;
        row.tickFont=ax.FontSize;row.legendBounds=lg.Position;row.output=pdf;
        fid=fopen(fullfile(output,[name '_verification.json']),'w');fprintf(fid,'%s',jsonencode(row));fclose(fid);
        rows{j}=row;
        if ~strcmp(getenv('BPA_FIGURE_STAGE_ONLY'),'1')
            staged=[destinations{j} '.pending'];copyfile(pdf,staged,'f');movefile(staged,destinations{j},'f');
            if j==1,home=fullfile(diss,'Figures','30-torque');else,home=fullfile(paper,'Figures');end
            component=fullfile(home,'Components','figure_polish_20261030');
            if ~isfolder(component),mkdir(component);end
            copyfile(fullfile(output,[name '.fig']),fullfile(component,[name '.fig']),'f');
        end
    end
    info.profiles=rows;close(fig);return
end
pdf=fullfile(candidate,[name '.pdf']);
exportgraphics(fig,pdf,'ContentType','vector');
exportgraphics(fig,fullfile(candidate,[name '_preview.png']),'Resolution',180);
savefig(fig,fullfile(candidate,[name '.fig']));
fid=fopen(fullfile(candidate,[name '_verification.json']),'w');fprintf(fid,'%s',jsonencode(info));fclose(fid);
if ~strcmp(getenv('BPA_FIGURE_STAGE_ONLY'),'1')
    for k=1:numel(destinations),copyfile(pdf,destinations{k},'f');end
    component=fullfile(diss,'Figures','30-torque','Components','figure_polish_20261030');
    if ~isfolder(component),mkdir(component);end
    copyfile(fullfile(candidate,[name '.fig']),fullfile(component,[name '.fig']),'f');
end
close(fig);
end

function h=byColor(objects,type,color)
hit=false(size(objects));
for k=1:numel(objects)
    o=objects(k);
    if ~strcmp(o.Type,type)||isempty(o.XData)||~any(isfinite(o.XData(:))),continue;end
    if strcmp(type,'line'),c=o.Color;else,c=o.MarkerFaceColor;if ischar(c)||isstring(c),c=o.CData;end,end
    hit(k)=isnumeric(c)&&numel(c)==3&&max(abs(c(:)-color(:)))<1e-8;
end
h=objects(hit);assert(isscalar(h),'Expected exactly one retained series of this style');
end
function p=plotParent(p,ax)
while ~isequal(p.Parent,ax),p=p.Parent;end
end
function out=orderedSeries(parent)
out=gobjects(0);
children=allchild(parent);
for k=1:numel(children)
    h=children(k);
    if isprop(h,'XData')&&~isempty(h.XData)&&any(isfinite(h.XData(:)))
        out(end+1)=h; %#ok<AGROW>
    elseif isprop(h,'Children')&&~isempty(h.Children)
        out=[out,orderedSeries(h)]; %#ok<AGROW>
    end
end
end
