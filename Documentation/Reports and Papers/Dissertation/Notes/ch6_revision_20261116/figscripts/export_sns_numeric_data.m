function export_sns_numeric_data()
% Export recorded arrays and exact icon vertices without graphics or Simulink.
root=fileparts(mfilename('fullpath'));
while ~endsWith(root,'Bipedal_Robot'),root=fileparts(root);end
src=fullfile(root,'Code','Matlab','SNS_Simscape'); addpath(src);
knee=load(fullfile(src,'results','knee_transient_reflex_20261028','replay_summary.mat'));
ix=knee.t>=29 & knee.t<=40;
knee.t=knee.t(ix); knee.C=knee.C(ix,:); knee.P=knee.P(ix,:);
knee=rmfield(knee,'F');
cup=load(fullfile(src,'results','sns_beer_results.mat'));
d=struct('matlab',version,'knee',knee,'bpa',sns_icon_geom('BPA'),'bio',sns_icon_geom('BIO'));
inv=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'sns_graphics_sources.json')));
for k=1:numel(inv.icons)
    code=inv.icons(k).display; name=inv.icons(k).name;
    if contains(name,'BPA')
        tag='BPA'; if contains(name,'_'), tag=strrep(name,'_',' '); end
        code=sns_icon_muscle(tag,true);
    elseif strcmp(name,'BioMuscle')
        code=sns_icon_muscle('BIO',true);
    end
    d.icons(k)=maskData(name,code);
end
for k=1:numel(cup.runs)
    d.cup(k).label=cup.runs(k).label;
    for field={'th','mCup','F','A','Ft','At'}
        s=cup.runs(k).(field{1});
        d.cup(k).(field{1})=struct('t',s.Time(:),'y',s.Data(:));
    end
end
fid=fopen(fullfile(root,'tmp','ch6_20261116','sns_numeric_data.json'),'w');
fprintf(fid,'%s',jsonencode(d)); fclose(fid);
fprintf('Exported saved arrays and shared geometry in MATLAB %s. No graphics or model calls.\n',version);
end

function data=maskData(name,code)
global SNS_CAPTURE_ITEMS SNS_CAPTURE_COLOUR
SNS_CAPTURE_ITEMS=struct('kind',{},'x',{},'y',{},'colour',{},'text',{});
SNS_CAPTURE_COLOUR=[0 0 0];
Esyn=0; %#ok<NASGU>
code=strrep(code,'patch(','capturePatch(');
code=strrep(code,'plot(','capturePlot(');
code=strrep(code,'disp(','captureText(');
code=strrep(code,'color(','captureColour(');
eval(code);
data=struct('name',name,'items',SNS_CAPTURE_ITEMS);
assert(~isempty(data.items),'Empty icon: %s',name);
end
function capturePatch(x,y,c)
global SNS_CAPTURE_ITEMS
SNS_CAPTURE_ITEMS(end+1)=struct('kind','patch','x',x(:)','y',y(:)','colour',c,'text','');
end
function capturePlot(x,y)
global SNS_CAPTURE_ITEMS SNS_CAPTURE_COLOUR
SNS_CAPTURE_ITEMS(end+1)=struct('kind','line','x',x(:)','y',y(:)','colour',SNS_CAPTURE_COLOUR,'text','');
end
function captureText(s)
global SNS_CAPTURE_ITEMS
SNS_CAPTURE_ITEMS(end+1)=struct('kind','text','x',0,'y',0,'colour',[0 0 0],'text',s);
end
function captureColour(c)
global SNS_CAPTURE_COLOUR
if strcmp(c,'black'),SNS_CAPTURE_COLOUR=[0 0 0];else,SNS_CAPTURE_COLOUR=[1 1 1];end
end
