% Dig_FlxExcels_20260913.m
% Corrected #3 audit: label-search the 2026_06 pinned excel (7 Flx sheets) for
% resting length / max contraction / tendon length / min+max measured muscle
% length, then compare against the Plot-script constants and the fitted data
% (min Lm_h per test from the kf build mats).
here = fileparts(mfilename('fullpath'));
root = fileparts(here);                                  %2022_02_Festo
F = fullfile(fileparts(root), '2026_06_Festo', 'Results_table_10mm_pinned.xlsx');

fprintf('=== 2026_06_Festo\\Results_table_10mm_pinned.xlsx, Flx sheets ===\n');
for s = 1:7
    sh = sprintf('FlxTest10mm_%d', s);
    C = readcell(F, 'Sheet', sh, 'Range', 'A1:P30');
    rest = NaN; kmax = NaN; ten = NaN; mn = NaN; mx = NaN;
    for rr = 1:size(C,1)
        for cc = 1:size(C,2)
            v = C{rr,cc};
            if ischar(v) || isstring(v)
                t = lower(char(v));
                if contains(t, 'resting')
                    rest = firstnum(C, rr, cc, 'right');
                elseif contains(t, 'max contraction')
                    kmax = firstnum(C, rr, cc, 'right');
                elseif contains(t, 'tendon')
                    ten = firstnum(C, rr, cc, 'right');
                elseif contains(t, 'min length')
                    mn = firstnum(C, rr, cc, 'both');
                elseif contains(t, 'max length')
                    mx = firstnum(C, rr, cc, 'both');
                end
            end
        end
    end
    fprintf('%s : rest=%s kmax=%s tendon=%s | minLm=%s maxLm=%s\n', sh, ...
        vs(rest), vs(kmax), vs(ten), vs(mn), vs(mx));
end

fprintf('\n=== script constants + min/max of the ACTUAL fitted data (Lm_h) ===\n');
mats = {'Plot_KneeFlxPin10mm_48cm.mat','Plot_KneeFlxPin10mm_46cm.mat', ...
        'Plot_KneeFlxPin10mm_47cm.mat','Plot_KneeFlxPin10mm_40cm.mat', ...
        'Plot_KneeFlxPin10mm_42cm.mat'};
srest = [485 457 479 406 410];
skmax = [398 387 403 338 335];
sten  = [0 0 0 35 0];
for k = 1:5
    S = load(fullfile(root, mats{k}), 'InflatedLength');
    Lh = S.InflatedLength;
    if iscell(Lh), Lh = cell2mat(Lh); end
    Lh = Lh(:) * 1000;   %mm
    fprintf('test %d (%s): script rest=%d kmax=%d tendon=%d | Lm_h min=%.1f max=%.1f | min-kmax = %+.1f mm %s\n', ...
        k, mats{k}(end-6:end), srest(k), skmax(k), sten(k), min(Lh), max(Lh), ...
        min(Lh)-skmax(k), iff(min(Lh) < skmax(k), '<<< CONCERN (min < kmax)', '(ok)'));
end

function s = vs(v)
if isnumeric(v) && isscalar(v) && ~isnan(v), s = num2str(v); else, s = '--'; end
end

function v = firstnum(C, rr, cc, dir)
%first numeric cell scanning right (and down for 'both') from the label
v = NaN;
if strcmpi(dir, 'right') || strcmpi(dir, 'both')
    for kk = 1:4
        if cc+kk > size(C,2), break; end
        w = C{rr, cc+kk};
        if isnumeric(w) && isscalar(w) && ~isnan(w), v = w; return; end
    end
end
if strcmpi(dir, 'both')
    for kk = 1:3
        if rr+kk > size(C,1), break; end
        w = C{rr+kk, cc};
        if isnumeric(w) && isscalar(w) && ~isnan(w), v = w; return; end
    end
end
end

function s = iff(c, a, b)
if c, s = a; else, s = b; end
end
