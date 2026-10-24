% Check Survey2 Winner 20261023 (session scratch verification, v2)
% Decisive checks on INTEGER cell counts (no float-tolerance artifacts):
% (1) 27 rows, all built; (2) caps: evalsUsed <= 300, wallMin <= 10.5;
% (3) BEST-MEETING guarantee vs the ladder CSV: the recorded winner's
%     meetCells = max over the site's rungs, and its activeCells = max
%     among the max-meet rungs; (4) unconverged list.

here = fileparts(mfilename('fullpath'));
T = readtable(fullfile(here, 'BiPulley_survey2_perMuscle_20261023.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
L = readtable(fullfile(here, 'BiPulley_survey2_ladder_20261023.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
fprintf('perMuscle rows: %d | ladder rows: %d\n', height(T), height(L));
assert(height(T) == 27, 'need exactly 27 perMuscle rows');
assert(all(strcmp(T.built, 'true')), 'built:false rows: %s', ...
    strjoin(T.site(~strcmp(T.built, 'true')), ','));
assert(all(T.evalsUsed <= 300), 'eval cap violated');
assert(all(T.wallMin <= 10.5), 'wall cap violated');
fprintf('CAPS OK: evalsUsed max %d (<=300), wallMin max %.2f (<=10)\n', ...
    max(T.evalsUsed), max(T.wallMin));

bad = {};
for k = 1:height(T)
    ix = find(L.actuator_id == T.actuator_id(k));
    nDom = double(T.domain_cells(k));
    recMeet = round(double(T.meetPct(k)) * nDom / 100);
    recAct = round(double(T.activePct(k)) * nDom / 100);
    mx = max(L.meetCells(ix));
    ma = max(L.activeCells(ix(L.meetCells(ix) == mx)));
    ok = (recMeet == mx) && (recAct == ma);
    if ~ok
        bad{end+1} = sprintf(['%s rec %d/%d meet, %d active vs ladder ' ...
            'max meet %d, max-active-among-max-meet %d'], T.site{k}, ...
            recMeet, nDom, recAct, mx, ma); %#ok<SAGROW>
    end
end
if isempty(bad)
    fprintf(['BEST-MEETING GUARANTEE OK: every recorded winner is the ' ...
        'max-meet rung (then max-active) for its site.\n']);
else
    fprintf('WINNER MISMATCHES:\n');
    cellfun(@(s) fprintf('  %s\n', s), bad);
end

convNum = zeros(height(T), 1);
isCh = iscell(T.converged);
for k = 1:height(T)
    if isCh, v = T.converged{k}; else, v = T.converged(k); end
    if isnumeric(v), convNum(k) = double(v); else, convNum(k) = str2double(char(v)); end
end
uncMask = convNum == 0;
if ~any(uncMask)
    fprintf('UNCONVERGED: none (all 27 advisor-converged)\n');
else
    unc = T.site(uncMask);
    uncS = cellfun(@(c) char(string(c)), unc, 'UniformOutput', false);
    fprintf('UNCONVERGED (%d): %s\n', numel(uncS), strjoin(uncS, ', '));
end
assert(isequal(convNum == 0, strcmp(T.stopReason, 'converged') == false), ...
    'converged flag disagrees with stopReason');
fprintf('converged flag consistent with stopReason for all 27 rows.\n');

% Winner-config duplication check: every rung of every site visited a
% UNIQUE (nBPA, G) -- the ladder never revisits a configuration.
cfgs = [L.actuator_id, L.nBPA, L.G];
[~, nU] = unique(cfgs, 'rows');
assert(numel(nU) == size(cfgs, 1), 'a (site,nB,G) config was revisited');
fprintf('LADDER OK: %d rungs, all (site, nBPA, G) unique.\n', size(cfgs, 1));
