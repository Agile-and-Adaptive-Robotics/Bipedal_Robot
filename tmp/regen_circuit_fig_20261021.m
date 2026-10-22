% Regenerate the journal circuit figure with the new side-view muscle
% symbol, then refresh the dissertation chapter copies.
here = 'C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape';
figDir = 'C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\60-sns';
addpath(here);
cd(here);
sns_draw_circuit;   % -> figures\KneeReflex_circuit.{png,pdf,svg}
copyfile(fullfile(here, 'figures', 'KneeReflex_circuit.pdf'), ...
         fullfile(figDir, 'KneeReflex_circuit.pdf'));
copyfile(fullfile(here, 'figures', 'KneeReflex_circuit.png'), ...
         fullfile(figDir, 'components', 'KneeReflex_circuit.png'));
fprintf('CIRCUIT FIG regenerated\n');
