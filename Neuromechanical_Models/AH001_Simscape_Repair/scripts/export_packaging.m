% Export saved survey fields only. No constructor, optimization, or simulation.
here=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(fileparts(here)));
addpath(fullfile(repo,'Code','Matlab','Robot_Data'));
addpath(fullfile(repo,'Code','Matlab','Mesh_Optimization'));
pkg=fullfile(repo,'Code','Matlab','Mesh_Optimization','Results','KompareXXVII');
out=fullfile(fileparts(here),'packaging','build'); if ~isfolder(out), mkdir(out); end
c=readtable(fullfile(pkg,'route_configurations.csv'),'TextType','string');
c=c(c.side=="selected",:); routes=cell(height(c),1);
for k=1:height(c)
    [~,n,e]=fileparts(c.mat_file(k));
    f=fullfile(repo,'Code','Matlab','Mesh_Optimization','Results','BiPulley_surveys',n+e);
    s=load(f); o=s.g.obj;
    r=struct('site',c.site(k),'spec',s.spec,'T',s.T,'saved_strain',s.g.strain);
    names={'RestingL','Kmax','TendonL','FittingLength','Diameter','L_p','sContraction','PulleySlack','PulleyInfeasible','PulleyGain','PulleyExitIndex','ActiveCrossing','FsingleBPA'};
    for j=1:numel(names), r.(names{j})=o.(names{j}); end
    r.source=char(f); routes{k}=r;
end
fid=fopen(fullfile(out,'saved_routes.json'),'w'); fwrite(fid,jsonencode(routes),'char'); fclose(fid);
fprintf('Exported %d saved routes; MATLAB %s\n',numel(routes),version);
