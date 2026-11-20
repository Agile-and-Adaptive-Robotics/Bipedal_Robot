function read_plot_data(destination)
% Reproduce the plotted flexor trials and read the marker properties of the
% physical validation figure. The optimizers are not run.
here = fileparts(mfilename('fullpath'));
repo = here;
while ~isfolder(fullfile(repo,'Testing_Data','2022_02_Festo'))
    parent = fileparts(repo);
    assert(~strcmp(parent,repo),'Cannot locate Testing_Data.');
    repo = parent;
end
if nargin == 0
    destination = fullfile(here,'plot_data.mat');
end
initial = pwd;
guard = onCleanup(@() cd(initial));
addpath(fullfile(repo,'Code','Matlab','Functions'), ...
        fullfile(repo,'Code','Matlab','Functions','ModernRobotics'), ...
        fullfile(repo,'Code','Matlab','Robot_Data'));
cd(fullfile(repo,'Testing_Data','2022_02_Festo'));
pick = load('minimizeFlxPin10_results.mat','filtered_results','xCols');
terms = pick.filtered_results(77,pick.xCols);
assert(max(abs(terms(1:3)./[0.003937370 39978.73 14733.75]-1)) < 1e-6);
[~,rigid] = minimizeFlxPin(0,Inf,Inf,1:5);
[~,adopted] = minimizeFlxPin(terms(1),terms(2),terms(3),1:5);
order = [1 5 2 4];
for k=1:numel(order)
    source = rigid(order(k));
    estimate = adopted(order(k));
    strain = source.strain(:);
    evalc('force = HuntEq(strain,0.171,source.P,0);');
    hunt = -source.mA(:).*max(force(:,2),0);
    hunt(strain>0.171) = 0;
    hunt(strain<0) = NaN;
    flex(k) = struct('angle',source.Ak(:),'measured_angle',source.Aexp(:), ...
        'measured',source.Mexp(:),'hybrid_angle',source.A_h(:),'hybrid',source.M_h(:), ...
        'unoptimized',source.M(:),'hunt',hunt,'optimized',estimate.M_p(:,3)); %#ok<AGROW>
end
paperFigs = fileparts(fileparts(here));
native = fullfile(paperFigs,'Components','figure_polish_20261030','Redesigned_flexor_pressure_publication.fig');
f = openfig(native,'invisible');
closeFig = onCleanup(@() close(f));
ax = findall(f,'Type','axes');
sc = findall(ax,'Type','scatter');
props = struct([]);
for k = 1:numel(sc)
    h = sc(k);
    if numel(h.XData)~=15
        continue
    end
    entry = struct('name',h.DisplayName,'face',h.MarkerFaceColor, ...
        'edge',h.MarkerEdgeColor,'alpha',h.MarkerFaceAlpha,'edgeAlpha',h.MarkerEdgeAlpha, ...
        'size',h.SizeData,'lineWidth',h.LineWidth);
    props = [props,entry]; %#ok<AGROW>
end
assert(numel(props)==2);
save(destination,'flex','props','-v7');
end
