%inspect_figs.m — dump structure of the exemplar .fig files + BPASet mats (read-only)
%Run: matlab -batch "run('D:/GitHub/Bipedal_Robot/tmp/paper_bridge/figscripts/inspect_figs.m')"
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');

%% ---- helper dump function ----
function dump_fig(fdir, fname)
    fprintf('\n################ %s ################\n', fname);
    try
        h = openfig(fullfile(fdir, fname), 'reuse');
    catch ME
        fprintf('  OPEN FAILED: %s\n', ME.message); return;
    end
    axList = findall(h, 'Type', 'axes');
    for k = 1:numel(axList)
        ax = axList(k);
        kids = ax.Children;
        if isempty(kids), continue; end
        fprintf('--- AXES #%d (pos=[%s]) tag=%s title=[%s] xl=[%s] yl=[%s]\n', ...
            k, num2str(round(ax.Position*100)/100), ax.Tag, ...
            ax.Title.String, ax.XLabel.String, ax.YLabel.String);
        fprintf('    XLim=[%s] YLim=[%s] FontSize=%s FontName=%s LW=%s\n', ...
            num2str(ax.XLim), num2str(ax.YLim), num2str(ax.FontSize), ax.FontName, num2str(ax.LineWidth));
        for c = 1:numel(kids)
            o = kids(c);
            try
                dn = o.DisplayName; if isempty(dn), dn='(none)'; end
            catch, dn='?'; end
            try, col = num2str(round(o.Color*1000)/1000); catch, col='?'; end
            try, ls = o.LineStyle; catch, ls=''; end
            try, lw = num2str(o.LineWidth); catch, lw='?'; end
            try, mk = o.Marker; catch, mk=''; end
            try, mfc = num2str(round(o.MarkerFaceColor*1000)/1000); catch, mfc=''; end
            try
                xd = o.XData; yd = o.YData;
                if isnumeric(xd)
                    fprintf('    [%s] "%s" col=%s ls=%s lw=%s mk=%s mfc=%s n=%d x[%.3g..%.3g] y[%.4g..%.4g]\n', ...
                        class(o), dn, col, ls, lw, mk, mfc, numel(xd), min(xd), max(xd), min(yd(:)), max(yd(:)));
                    if numel(xd) <= 30
                        fprintf('        X='); fprintf('%.4g,', xd); fprintf('\n        Y='); fprintf('%.5g,', yd); fprintf('\n');
                    else
                        fprintf('        Xfirst3='); fprintf('%.4g,', xd(1:3));
                        fprintf(' Xlast3='); fprintf('%.4g,', xd(end-2:end));
                        fprintf('\n        Yfirst3='); fprintf('%.5g,', yd(1:3));
                        fprintf(' Ylast3='); fprintf('%.5g,', yd(end-2:end)); fprintf('\n');
                    end
                end
            catch, fprintf('    [%s] "%s" (no XY data)\n', class(o), dn); end
        end
    end
    lg = findall(h, 'Type', 'legend');
    for k = 1:numel(lg)
        fprintf('--- LEGEND loc=%s fs=%s strings: ', lg(k).Location, num2str(lg(k).FontSize));
        disp(lg(k).String);
    end
    close(h);
end

dump_fig(festodir, 'FlxPin_group.fig');
dump_fig(festodir, 'FlxGrp.fig');
dump_fig(festodir, 'Ext10mm_52cm.fig');

%% ---- Plot data mats ----
T = load('Plot_KneeFlx_20mm_42cm.mat', 'Angle', 'Torque', 'TorqueHand');
if isfield(T,'Angle'), fprintf('\nPlot_KneeFlx_20mm Angle(%d): ', numel(T.Angle)); fprintf('%.4g,', T.Angle); fprintf('\n'); end
if isfield(T,'Torque'), fprintf('Plot_KneeFlx_20mm Torque(%d): ', numel(T.Torque)); fprintf('%.5g,', T.Torque); fprintf('\n'); end
if isfield(T,'TorqueHand'), fprintf('Plot_KneeFlx_20mm TorqueHand(%d): ', numel(T.TorqueHand)); fprintf('%.5g,', T.TorqueHand); fprintf('\n'); end

fprintf('\nDONE_INSPECT\n');
