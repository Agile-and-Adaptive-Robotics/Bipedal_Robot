%% probe_font_params_20260920.m — find the real font-size parameter names
%% for models, blocks, and annotations on THIS release (R2025b).

mdl = 'fnt_probe';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
add_block('simulink/Sources/Constant', [mdl '/C'], 'Position', [30 30 60 60]);
add_block('simulink/Sinks/Out1', [mdl '/O'], 'Position', [90 30 120 44]);
add_line(mdl, 'C/1', 'O/1');
add_block('built-in/Note', [mdl '/note'], 'Position', [30 90 60 120], 'Text', 'hi');

mp = get_param(mdl, 'ObjectParameters');
fprintf('MODEL params containing "Font":\n');
fn = fieldnames(mp);
for k = 1:numel(fn)
    if contains(fn{k}, 'Font', 'IgnoreCase', true)
        try
            v = get_param(mdl, fn{k});
            fprintf('  %-24s = %s\n', fn{k}, mat2str(v));
        catch
            fprintf('  %-24s = <unreadable>\n', fn{k});
        end
    end
end

fprintf('BLOCK params containing "Font" (Constant C):\n');
bp = get_param([mdl '/C'], 'ObjectParameters');
fn = fieldnames(bp);
for k = 1:numel(fn)
    if contains(fn{k}, 'Font', 'IgnoreCase', true)
        try
            v = get_param([mdl '/C'], fn{k});
            fprintf('  %-24s = %s\n', fn{k}, mat2str(v));
        catch
            fprintf('  %-24s = <unreadable>\n', fn{k});
        end
    end
end

fprintf('ANNOTATION params containing "Font" (note):\n');
ap = get_param([mdl '/note'], 'ObjectParameters');
fn = fieldnames(ap);
for k = 1:numel(fn)
    if contains(fn{k}, 'Font', 'IgnoreCase', true)
        try
            v = get_param([mdl '/note'], fn{k});
            fprintf('  %-24s = %s\n', fn{k}, mat2str(v));
        catch
            fprintf('  %-24s = <unreadable>\n', fn{k});
        end
    end
end

% try the likely model-level setters
for cand = {'DefaultBlockFontSize', 'BlockFontSize', 'FontSize'}
    try
        set_param(mdl, cand{1}, '12');
        fprintf('SET OK: model %s\n', cand{1});
    catch ME
        fprintf('set %s failed: %s\n', cand{1}, ME.message);
    end
end
% block-level
try
    set_param([mdl '/C'], 'FontSize', '12');
    fprintf('SET OK: block FontSize (C now %s)\n', get_param([mdl '/C'], 'FontSize'));
catch ME
    fprintf('block FontSize set failed: %s\n', ME.message);
end
close_system(mdl, 0);
fprintf('PROBE DONE\n');
