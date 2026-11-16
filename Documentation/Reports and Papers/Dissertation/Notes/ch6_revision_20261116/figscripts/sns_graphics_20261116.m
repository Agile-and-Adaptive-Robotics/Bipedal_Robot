function sns_graphics_20261116()
% Publish current SNS glyphs and recorded reflex tests without changing models.
root = fileparts(mfilename('fullpath'));
while ~endsWith(root,'Bipedal_Robot'), root = fileparts(root); end
src = fullfile(root,'Code','Matlab','SNS_Simscape');
out = fullfile(root,'Documentation','Reports and Papers','Dissertation','Figures','60-sns');
here = fileparts(mfilename('fullpath'));
addpath(src);
set(groot,'DefaultFigureVisible','off','DefaultAxesFontName','Arial', ...
    'DefaultAxesFontSize',10.5,'DefaultTextFontName','Arial','DefaultTextFontSize',10.5, ...
    'DefaultAxesLineWidth',1,'DefaultLineLineWidth',1.5);
fprintf('MATLAB %s\n',version);
inventory = jsondecode(fileread(fullfile(here,'sns_graphics_sources.json')));
assert(numel(inventory.icons)==12);
geometryCheck();
libraryFigure(inventory,out);
knee = load(fullfile(src,'results','knee_transient_reflex_20261028','replay_summary.mat'));
kneeFigure(knee,out);
cup = load(fullfile(src,'results','sns_beer_results.mat'));
cupFigure(cup.runs,out);
fid=fopen(fullfile(here,'graphics_metrics.json'),'w');
metrics=struct('matlab',version,'knee',knee.summary, ...
    'cup_final_on_deg',rad2deg(cup.runs(1).th.Data(end)), ...
    'cup_final_off_deg',rad2deg(cup.runs(2).th.Data(end)), ...
    'cup_final_load_kg',cup.runs(1).mCup.Data(end));
fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true)); fclose(fid);
fprintf('All three figures exported. No model, result, or parameter file saved.\n');
end

function geometryCheck()
g=sns_icon_geom('BPA');
code=sns_icon_muscle('BPA',true);
p=regexp(code,'plot\(\[(.*?)\], \[(.*?)\]\);','tokens');
assert(numel(p)==size(g.braid,1));
for k=1:numel(p)
    x=sscanf(p{k}{1},'%f'); y=sscanf(p{k}{2},'%f');
    assert(max(abs([x(:)' y(:)']-g.braid(k,[1 3 2 4])))<6e-5);
end
q=regexp(code,'patch\(\[.*?\], \[(.*?)\], \[0 0 0\]\);','tokens','once');
y=sscanf(q{1},'%f');
assert(max(abs(y(:)'-(g.bodyY+g.edge*sign(g.bodyY))))<6e-5);
fprintf('Icon geometry PASS: %d braid segments and outer body match shared vertices.\n',numel(p));
end

function libraryFigure(inv,out)
names={'NonSpikingNeuron','NonSpikingSynapse','SynSum','SpikingLIFNeuron', ...
    'IaMuscleSpindle','IbGolgiTendon','MuscleActivation','BioMuscle','BPAForce', ...
    'BPA_10mm','BPA_20mm','BPA_40mm'};
labels={'Nonspiking neuron','Graded synapse','Synaptic summation', ...
    'Spiking LIF neuron','Ia spindle encoder','Ib force encoder', ...
    'Activation dynamics','Biological muscle','Simple force element', ...
    '10 mm BPA','20 mm BPA','40 mm BPA'};
sub={'Membrane voltage and drive','Excitation or inhibition','Conductance pairs', ...
    'Threshold and reset','Length and stretch velocity','Tendon-force signal', ...
    'First-order filter','Thelen-style formulation','Reduced-order force law', ...
    'Pressure and length','Pressure and length','Pressure and length'};
f=newFig(14,15.2);
for k=1:12
    row=floor((k-1)/3); col=mod(k-1,3);
    cx=.04+col*.326; top=.985-row*.247;
    annotation(f,'textbox',[cx top-.232 .305 .23], 'String','', ...
        'EdgeColor',[.8 .8 .8],'LineWidth',.6,'BackgroundColor','w');
    ax=axes(f,'Position',[cx+.071 top-.126 .165 .115]);
    axis(ax,'equal'); axis(ax,[-1.15 1.15 -1.12 1.12]); axis(ax,'off'); hold(ax,'on');
    idx=find(strcmp({inv.icons.name},names{k}),1);
    code=inv.icons(idx).display;
    if contains(names{k},'BPA')
        txt='BPA'; if contains(names{k},'_'),txt=strrep(names{k},'_',' '); end
        code=sns_icon_muscle(txt,true);
        set(ax,'Position',[cx+.035 top-.108 .235 .079]);
        axis(ax,[-1.05 1.05 -.41 .41]);
    elseif strcmp(names{k},'BioMuscle')
        code=sns_icon_muscle('BIO',true);
        set(ax,'Position',[cx+.035 top-.108 .235 .079]);
        axis(ax,[-1.05 1.05 -.41 .41]);
    end
    renderMask(ax,code);
    textBox(f,[cx+.006 top-.171 .293 .034],labels{k},10.5,'bold');
    textBox(f,[cx+.005 top-.207 .295 .033],sub{k},10,'normal');
end
saveFigure(f,out,'sns_library_20261116');
altwrite(out,'sns_library_20261116', ...
    'Twelve SNS-Simscape library glyphs matching the block selection in defense slide 69. Circles distinguish graded and spiking waveforms; a compact synapse has an open excitatory terminal; spindle and tendon-force encoders have distinct capsules; activation is a pentagon; biological muscle is pennate; BPAs have braided sleeves and clamp fittings. The simple force element is distinguished from diameter-specific pressure-length models.');
end

function renderMask(ax,code)
axes(ax); hold(ax,'on');
Esyn=0; iconColour=[0 0 0]; %#ok<NASGU>
code=regexprep(code,'color\(''black''\);','iconColour=[0 0 0];');
code=regexprep(code,'color\(''white''\);','iconColour=[1 1 1];');
code=regexprep(code,'plot\(([^;]+)\);','line($1,''Color'',iconColour,''LineWidth'',1.2);');
code=regexprep(code,'disp\(([^;]+)\);','text(0,0,$1,''HorizontalAlignment'',''center'',''FontSize'',10.5,''FontWeight'',''bold'');');
code=regexprep(code,'port_label\([^;]+\);','');
eval(code);
set(findall(ax,'Type','patch'),'EdgeColor','none');
end

function kneeFigure(s,out)
f=newFig(14,18.7);
ax=axes(f,'Position',[.025 .645 .95 .335]);
reflexCircuit(ax,{'Extensor','Flexor'},'A  Implemented antagonist feedback');
ax=axes(f,'Position',[.15 .51 .79 .105]);
sel=s.t>=29.85&s.t<=30.25; x=s.t(sel)-30;
plot(ax,x,s.P(sel,9),'Color',[.2 .2 .2]); xlim(ax,[-.15 .25]); ylim(ax,[-.1 2.3]);
ylabel(ax,{'Ia current','(nA)'});
title(ax,'B  Sensory-current stimulus','FontSize',11,'HorizontalAlignment','left'); styleAxes(ax);
ax=axes(f,'Position',[.15 .295 .79 .14]);
sel=s.t>=29&s.t<=40; x=s.t(sel)-30;
plot(ax,x,-rad2deg(s.P(sel,1)),'Color',blue); hold(ax,'on');
plot(ax,x,-rad2deg(s.C(sel,1)),'--','Color',[.3 .3 .3]); xlim(ax,[-1 10]);
xline(ax,0,':','Color',[.5 .5 .5]); ylabel(ax,{'Knee angle','(deg)'});
title(ax,'C  Extension and recovery','FontSize',11,'HorizontalAlignment','left');
legend(ax,{'Stimulus','No stimulus'},'Location','southeast','Box','off','FontSize',10.5);
text(ax,.99,.92,sprintf('Peak +%.2f deg at %d ms',s.summary.paper_extension_excursion_deg,round(1000*s.summary.peak_latency_s)), ...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',10);
styleAxes(ax);
ax=axes(f,'Position',[.15 .095 .79 .12]);
plot(ax,x,s.P(sel,3),'Color',blue); hold(ax,'on');
plot(ax,x,s.P(sel,4),'--','Color',orange); xlim(ax,[-1 10]);
xline(ax,0,':','Color',[.5 .5 .5]); ylabel(ax,{'MN voltage','(mV)'}); xlabel(ax,'Time from pulse onset (s)');
title(ax,'D  Motor-neuron response','FontSize',11,'HorizontalAlignment','left');
legend(ax,{'Extensor','Flexor'},'Location','northeast','Box','off','FontSize',10.5); styleAxes(ax);
saveFigure(f,out,'sns_knee_reflex_20261116');
altwrite(out,'sns_knee_reflex_20261116', ...
    'Panel A shows four sensory-neuron channels and both motor neurons, with all six synapses: each Ia channel excites its own motor neuron and inhibits its antagonist, and each Ib channel inhibits its own motor neuron. Panel B shows the 2 nA, 30 ms sensory-current pulse. Panel C shows the knee extending 1.09 degrees at 111 ms and returning toward the no-stimulus condition. Panel D shows extensor and flexor motor-neuron voltages. The saved graded-response configuration uses simple force elements, and flexor drive is already at its floor.');
end

function cupFigure(runs,out)
f=newFig(14,19.1);
ax=axes(f,'Position',[.04 .705 .92 .27]);
reflexCircuit(ax,{'Biceps','Triceps'},'A  Implemented loaded-arm feedback');
ax=axes(f,'Position',[.07 .445 .34 .22]);
armMechanism(ax,rad2deg(runs(1).th.Data(end)),rad2deg(runs(2).th.Data(end)));
title(ax,'B  Loaded-arm geometry','FontSize',11);
ax=axes(f,'Position',[.57 .48 .36 .15]);
plot(ax,runs(1).mCup.Time,runs(1).mCup.Data,'Color',[.2 .2 .2]);
xlim(ax,[0 10]); ylim(ax,[0 .66]); ylabel(ax,'Added mass (kg)'); xlabel(ax,'Time (s)');
title(ax,'C  Applied load','FontSize',11); styleAxes(ax);
ax=axes(f,'Position',[.14 .255 .8 .135]);
plot(ax,runs(1).th.Time,-rad2deg(runs(1).th.Data),'Color',blue); hold(ax,'on');
plot(ax,runs(2).th.Time,-rad2deg(runs(2).th.Data),'--','Color',orange);
xlim(ax,[0 10]); ylim(ax,[-9 .5]); ylabel(ax,{'Forearm elevation','(deg; 0 = level)'});
title(ax,'D  Feedback limits arm sag','FontSize',11,'HorizontalAlignment','left');
legend(ax,{'Feedback on','Feedback off'},'Location','southwest','Box','off','FontSize',10.5); styleAxes(ax);
ax=axes(f,'Position',[.14 .065 .8 .125]);
plot(ax,runs(1).A.Time,runs(1).A.Data,'Color',blue); hold(ax,'on');
plot(ax,runs(1).At.Time,runs(1).At.Data,'--','Color',orange);
xlim(ax,[0 10]); ylim(ax,[0 .55]); ylabel(ax,{'Normalized','pressure command'}); xlabel(ax,'Time (s)');
title(ax,'E  Feedback-on antagonist commands','FontSize',11,'HorizontalAlignment','left');
legend(ax,{'Biceps','Triceps'},'Location','southwest','Box','off','FontSize',10.5); styleAxes(ax);
saveFigure(f,out,'sns_cup_reflex_20261116');
altwrite(out,'sns_cup_reflex_20261116', ...
    'Panel A contains all six actual loaded-arm feedback synapses, including separate biceps and triceps Ib autogenic inhibition. Panel B shows the elbow and cup, side-view braided BPA glyphs, and the small recorded feedback-on/off final arm inclinations. Panel C shows the ramp to 0.625 kg over 10 s. Panel D plots elevation: feedback-on arm droops about 2 degrees versus about 7.85 degrees with feedback off. Panel E shows biceps and triceps normalized pressure commands in the feedback-on recording.');
end

function reflexCircuit(ax,names,heading)
cla(ax); hold(ax,'on'); axis(ax,[0 14 0 5.9]); axis(ax,'off');
text(ax,.1,5.65,heading,'FontWeight','bold','FontSize',11);
% Source neurons are drawn separately from the sensory transducers.
sx=[2.0 .8 10.0 12.8]; sy=4.5; my=1.45; mx=[3.0 11.0];
% Six disjoint semantic routes. Crossed segments have no junction markers.
edge(ax,[sx(1) sy-.38;mx(1)-.22 my+.45],true,blue);
edge(ax,[sx(2)+.38 2.8;mx(1)-.46 my],false,orange);
edge(ax,[sx(3) sy-.38;mx(2)-.22 my+.45],true,blue);
edge(ax,[sx(4)-.38 2.8;mx(2)+.46 my],false,orange);
edge(ax,[sx(1)+.38 sy;6.6 sy;6.6 my;mx(2)-.46 my],false,blue);
edge(ax,[sx(3)-.38 sy;7.5 sy;7.5 my-.8;mx(1) my-.8;mx(1) my-.46],false,blue);
for k=1:4
    label='Ia'; if ismember(k,[2 4]), label='Ib'; end
    y=sy; if ismember(k,[2 4]),y=2.8;end
    neuron(ax,[sx(k) y],.38,label,[.94 .97 .94]);
    side=1+(k>=3);
    text(ax,sx(k),y+.6,names{side},'HorizontalAlignment','center','FontSize',10.5);
end
for k=1:2
    neuron(ax,[mx(k) my],.46,'MN',[1 1 1]);
    labelx=mx(k)+(-1.3+2.6*(k-1));
    text(ax,labelx,my-.55,names{k},'HorizontalAlignment','center','FontSize',10.5);
end
% Legend uses shape and line style in addition to color.
edge(ax,[.3 .25;1.0 .25],true,[.2 .2 .2]);
text(ax,1.25,.25,'excitation','FontSize',10);
edge(ax,[4.5 .25;5.2 .25],false,[.2 .2 .2]);
text(ax,5.45,.25,'inhibition','FontSize',10);
text(ax,9.4,.25,'Ia: stretch   Ib: force','FontSize',10);
end

function neuron(ax,p,r,label,col)
t=linspace(0,2*pi,65);
patch(ax,p(1)+r*cos(t),p(2)+r*sin(t),col,'EdgeColor','k','LineWidth',1.3);
text(ax,p(1),p(2),label,'HorizontalAlignment','center','FontSize',10.5,'FontWeight','bold');
end

function edge(ax,p,exc,col)
style='-'; if ~exc, style='--'; end
plot(ax,p(:,1),p(:,2),style,'Color',col,'LineWidth',1.3);
u=p(end,:)-p(end-1,:); u=u/norm(u); n=[-u(2) u(1)]; q=p(end,:);
if exc
    v=[q-u*.17;q-n*.12;q+n*.12];
    patch(ax,v(:,1),v(:,2),'w','EdgeColor','k','LineWidth',1);
else
    t=linspace(0,2*pi,25); patch(ax,q(1)+.095*cos(t),q(2)+.095*sin(t),'k','EdgeColor','k');
end
end

function armMechanism(ax,on,off)
hold(ax,'on'); axis(ax,'equal'); axis(ax,[-.17 .44 -.44 .05]); axis(ax,'off');
sh=[-.06 0]; el=[-.06 -.27];
plot(ax,[sh(1) el(1)],[sh(2) el(2)],'-','Color',[.6 .6 .6],'LineWidth',12);
plot(ax,[el(1) el(1)+.3*cosd(on)],[el(2) el(2)-.3*sind(on)],'-','Color',blue,'LineWidth',8);
plot(ax,[el(1) el(1)+.3*cosd(off)],[el(2) el(2)-.3*sind(off)],'--','Color',orange,'LineWidth',1.5);
plot(ax,el(1),el(2),'ko','MarkerFaceColor','w','MarkerSize',7);
bpaGlyph(ax,[-.025 -.045],[.005 -.255],.028);
bpaGlyph(ax,[-.095 -.045],[-.09 -.28],.022);
hand=el+[.3*cosd(on) -.3*sind(on)];
patch(ax,hand(1)+[-.035 .035 .04 -.04],hand(2)+[0 0 .1 .1],[.93 .96 1], ...
    'EdgeColor',[.25 .25 .25],'LineWidth',1);
patch(ax,hand(1)+[-.035 .035 .038 -.038],hand(2)+[.005 .005 .075 .075],[.65 .85 .95], ...
    'EdgeColor','none');
text(ax,.06,-.07,'20 mm BPAs','FontSize',10);
text(ax,.09,-.395,sprintf('ON %.2f deg\nOFF %.2f deg',on,off),'FontSize',10,'HorizontalAlignment','center');
end

function bpaGlyph(ax,p1,p2,w)
g=sns_icon_geom('BPA'); u=(p2-p1)/norm(p2-p1); n=[-u(2) u(1)]; c=(p1+p2)/2;
trans=@(x,y)c+((x(:)/1.8)*norm(p2-p1)).*u+(y(:)*w/.6).*n;
p=trans(g.bodyX,g.bodyY); patch(ax,p(:,1),p(:,2),g.fill,'EdgeColor','k','LineWidth',.8);
for k=1:size(g.braid,1)
 p=trans(g.braid(k,[1 3]),g.braid(k,[2 4])); plot(ax,p(:,1),p(:,2),'k-','LineWidth',.5);
end
for q={g.capL,g.capR}
 p=trans(q{1}(:,1),q{1}(:,2)); patch(ax,p(:,1),p(:,2),g.cap,'EdgeColor','k','LineWidth',.5);
end
end

function f=newFig(w,h)
f=figure('Units','centimeters','Position',[1 1 w h],'Color','w','Visible','off');
end
function styleAxes(ax)
set(ax,'Box','off','TickDir','out','FontName','Arial','FontSize',10.5,'Layer','top');
grid(ax,'on'); ax.GridAlpha=.13;
end
function textBox(f,pos,txt,fs,weight)
annotation(f,'textbox',pos,'String',txt,'EdgeColor','none','FontName','Arial', ...
    'FontSize',fs,'FontWeight',weight,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','none','Margin',0);
end
function c=blue(), c=[68 119 170]/255; end
function c=orange(), c=[238 102 119]/255; end
function saveFigure(f,out,stem)
exportgraphics(f,fullfile(out,[stem '.pdf']),'ContentType','vector','BackgroundColor','white');
exportgraphics(f,fullfile(out,[stem '.png']),'Resolution',200,'BackgroundColor','white');
close(f); fprintf('Wrote %s\n',stem);
end
function altwrite(out,stem,txt)
fid=fopen(fullfile(out,[stem '.alt.txt']),'w'); fprintf(fid,'%s\n',txt); fclose(fid);
end
