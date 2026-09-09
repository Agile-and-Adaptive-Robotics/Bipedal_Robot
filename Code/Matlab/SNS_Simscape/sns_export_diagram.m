function sns_export_diagram(outDir, models, arrange)
% SNS_EXPORT_DIAGRAM  Print Simulink models to publication-support images.
%
%   sns_export_diagram()                       — KneeReflexDemo only, no rearrange
%   sns_export_diagram(outDir, models, arrange)
%
% Writes <model>_simulink.png (300 dpi) and <model>_simulink.pdf per model.
% The auto-arrange option (Simulink.BlockDiagram.arrangeSystem) rebuilds the
% layout algorithmically — off by default because the hand layout matches the
% paper circuit; use it only for quick internal figures.
% NOTE: the primary publication figure is sns_draw_circuit.m (vector redraw),
% not a Simulink screenshot.

if nargin < 1 || isempty(outDir), outDir = fullfile(fileparts(mfilename('fullpath')), 'figures'); end
if nargin < 2 || isempty(models), models = {'KneeReflexDemo'}; end
if nargin < 3, arrange = false; end
if ~exist(outDir, 'dir'), mkdir(outDir); end

for k = 1:numel(models)
    mdl = models{k};
    if ~bdIsLoaded(mdl), load_system(mdl); end
    if arrange
        Simulink.BlockDiagram.arrangeSystem(mdl);
    end
    pngPath = fullfile(outDir, [mdl '_simulink.png']);
    pdfPath = fullfile(outDir, [mdl '_simulink.pdf']);
    print(['-s' mdl], '-dpng', '-r300', pngPath);
    print(['-s' mdl], '-dpdf', '-bestfit', pdfPath);
    fprintf('Wrote %s and %s\n', pngPath, pdfPath);
end
end
