% Patch the live SNS_Library.slx muscle-family icons to the side-view
% symbols (Ben, 2026-10-21): McKibben BPA with braided sleeve + clamp
% fittings; pennate biological muscle. Display strings come from
% sns_icon_muscle.m, the same source the builders now use.
here = 'C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape';
addpath(here);
load_system(fullfile(here, 'SNS_Library.slx'));
set_param('SNS_Library', 'Lock', 'off');
targets = {'BPAForce', 'BPA'; 'BPA_10mm', 'BPA 10'; 'BPA_20mm', 'BPA 20'; ...
           'BPA_40mm', 'BPA 40'; 'BioMuscle', 'BIO'};
for k = 1:size(targets, 1)
    blk = ['SNS_Library/' targets{k, 1}];
    m = Simulink.Mask.get(blk);
    m.Display = sns_icon_muscle(targets{k, 2});
    % wide default aspect so the side-view drawing does not squash
    g = sns_icon_geom(targets{k, 2});
    p = get_param(blk, 'Position');
    h = round((p(3) - p(1)) / g.blockAspect);
    set_param(blk, 'Position', [p(1) p(2) p(3) p(2) + h]);
    fprintf('patched %s\n', targets{k, 1});
end
save_system('SNS_Library');
set_param('SNS_Library', 'Lock', 'on');
% proof render of the canvas with the new icons
open_system('SNS_Library');
drawnow; pause(1.0);
print('-sSNS_Library', '-dpng', '-r200', ...
      fullfile(here, 'results', 'pictures', 'library_icons_proof.png'));
close_system('SNS_Library', 0);
fprintf('PATCH DONE\n');
