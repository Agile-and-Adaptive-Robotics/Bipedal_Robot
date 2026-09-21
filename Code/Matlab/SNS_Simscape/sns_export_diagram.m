function sns_export_diagram(outDir, models, arrange)
% SNS_EXPORT_DIAGRAM  Print Simulink models as dissertation-ready figures.
%
%   sns_export_diagram()                        — KneeReflexDemo only
%   sns_export_diagram(outDir, models, arrange)
%   sns_export_diagram(figures, {'KneeReflexDemo','SNS_Deng_RG'}, false)
%
% Per model writes <model>_simulink.{eps,pdf,png} into outDir (default
% figures\):
%   .eps — color VECTOR (painters renderer). LaTeX/Overleaf include it with
%          plain \includegraphics{KneeReflexDemo_simulink.eps} (Overleaf
%          auto-converts eps->pdf at compile; a local pdflatex needs
%          \usepackage{epstopdf} or shell-escape, so the twin .pdf is the
%          zero-friction option — identical vector content, tight page).
%   .pdf — color vector, best-fit page.
%   .png — 300 dpi raster for slides / quick checks.
%
% Dissertation styling (Ben 2026-09-20: 10–12 pt, Rybak-style) applied to
% the in-memory diagram before printing; the .slx on disk is NEVER modified
% (nothing is saved, and original param values are restored after):
%   * model default font size -> 12 pt (block names, port labels)
%   * annotation font sizes normalized into 10–12 pt
%   * synapse-block names hidden — in the Rybak/Szczecinski diagram language
%     a synapse is a CONNECTION (marker carries the meaning), not a labelled
%     component
%   * canvas forced white
%
% The auto-arrange option (Simulink.BlockDiagram.arrangeSystem) rebuilds the
% layout algorithmatically — off by default because the hand layout matches
% the paper circuit; use it only for quick internal figures.
% NOTE: the primary hand-drawn publication circuit remains sns_draw_circuit.m;
% these exports are for showing the actual Simulink models.

if nargin < 1 || isempty(outDir), outDir = fullfile(fileparts(mfilename('fullpath')), 'figures'); end
if nargin < 2 || isempty(models), models = {'KneeReflexDemo'}; end
if nargin < 3, arrange = false; end
if ~exist(outDir, 'dir'), mkdir(outDir); end

synMaskTypes = {'SNS NonSpiking Synapse', 'SNS Deng Synapse'};  % hide-name set

for k = 1:numel(models)
    mdl = models{k};
    loadedBefore = bdIsLoaded(mdl);
    if ~loadedBefore, load_system(mdl); end

    % ---- styling (all restored below; nothing saved) ------------------------
    % R2025b probe (2026-09-20): the MODEL has no font parameter; each BLOCK
    % carries FontSize (-1 = inherit the ~10 pt diagram default). Set every
    % block to 12 pt and normalize annotations into 10-12 pt.
    set_param(mdl, 'ScreenColor', 'white');

    blks = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
    origBlkFont = nan(numel(blks), 1);
    for b = 1:numel(blks)
        try
            origBlkFont(b) = str2double(get_param(blks{b}, 'FontSize'));
            set_param(blks{b}, 'FontSize', '12');
        catch
        end
    end

    anns = find_system(mdl, 'SearchDepth', 1, 'Type', 'annotation');
    origAnn = nan(numel(anns), 1);
    for a = 1:numel(anns)
        try
            origAnn(a) = str2double(get_param(anns{a}, 'FontSize'));
            fs = origAnn(a);
            if isnan(fs) || fs < 10, fs = 11; end   % normalize into 10-12 pt
            fs = min(fs, 12);
            set_param(anns{a}, 'FontSize', num2str(fs));
        catch
        end
    end

    syns = find_system(mdl, 'LookUnderMasks', 'all', 'Type', 'Block');
    hideNames = {};
    for b = 1:numel(syns)
        try
            mt = get_param(syns{b}, 'MaskType');
        catch
            mt = '';
        end
        if any(strcmp(mt, synMaskTypes)) && strcmp(get_param(syns{b}, 'ShowName'), 'on')
            set_param(syns{b}, 'ShowName', 'off');
            hideNames{end+1} = syns{b}; %#ok<SAGROW>
        end
    end

    if arrange
        Simulink.BlockDiagram.arrangeSystem(mdl);
    end

    % ---- print ---------------------------------------------------------------
    % NOTE (R2025b, 2026-09-20): print('-s<mdl>','-depsc') is REFUSED for
    % Simulink systems ("'epsc' format is not supported with Simulink or
    % Stateflow printing"). The EPS is therefore produced by converting the
    % vector PDF with MATLAB's bundled ghostscript (eps2write device).
    epsPath = fullfile(outDir, [mdl '_simulink.eps']);
    pdfPath = fullfile(outDir, [mdl '_simulink.pdf']);
    pngPath = fullfile(outDir, [mdl '_simulink.png']);
    print(['-s' mdl], '-dpdf', '-bestfit', pdfPath);
    print(['-s' mdl], '-dpng', '-r300', pngPath);
    [gs, pdfcrop] = findEpsTools();
    if ~isempty(gs) && ~isempty(pdfcrop)
        % temp crop file must live in a WRITABLE folder — pdfcrop cannot move
        % its temp file to C:\ root (tempname('C:') failed 2026-09-20)
        tmpCrop = fullfile(outDir, [mdl '_crop_tmp.pdf']);
        [st1] = system(['"' pdfcrop '" "' pdfPath '" "' tmpCrop '"']);
        [st2, so] = system(['"' gs '" -dNOPAUSE -dBATCH -dSAFER -sDEVICE=eps2write' ...
            ' -sOutputFile="' epsPath '" "' tmpCrop '"']);
        delete(tmpCrop);
        if st1 == 0 && st2 == 0
            fprintf('Wrote %s.{eps,pdf,png} (tight vector EPS via pdfcrop+ghostscript)\n', ...
                fullfile(outDir, [mdl '_simulink']));
        else
            fprintf('Wrote %s.{pdf,png} — EPS conversion FAILED (pdfcrop st=%d, gs st=%d):\n%s\n', ...
                fullfile(outDir, [mdl '_simulink']), st1, st2, so);
        end
    else
        fprintf('Wrote %s.{pdf,png} — pdfcrop/ghostscript not found (MiKTeX bin), EPS skipped\n', ...
            fullfile(outDir, [mdl '_simulink']));
    end
    % ---- restore -------------------------------------------------------------
    for b = 1:numel(blks)
        if ~isnan(origBlkFont(b))
            set_param(blks{b}, 'FontSize', num2str(origBlkFont(b)));
        end
    end
    for b = 1:numel(hideNames)
        set_param(hideNames{b}, 'ShowName', 'on');
    end
    for a = 1:numel(anns)
        if ~isnan(origAnn(a))
            set_param(anns{a}, 'FontSize', num2str(origAnn(a)));
        end
    end
    if ~loadedBefore, close_system(mdl, 0); end
end
end

function [gs, pdfcrop] = findEpsTools()
% R2025b ships NO ghostscript and print(-depsc) refuses Simulink systems;
% MiKTeX's mgs.exe (real ghostscript) + pdfcrop.exe are the working EPS
% pipeline on this machine (2026-09-20).
gs = ''; pdfcrop = '';
miktexBin = fullfile(getenv('LOCALAPPDATA'), 'Programs', 'MiKTeX', 'miktex', 'bin', 'x64');
gsCands = { ...
    fullfile(matlabroot, 'bin', 'win64', 'gswin64c.exe'), ...
    fullfile(matlabroot, 'sys', 'ghostscript', 'bin', 'gswin64c.exe'), ...
    fullfile(miktexBin, 'mgs.exe')};
for k = 1:numel(gsCands)
    if isfile(gsCands{k}), gs = gsCands{k}; break; end
end
pc = fullfile(miktexBin, 'pdfcrop.exe');
if isfile(pc), pdfcrop = pc; end
end
