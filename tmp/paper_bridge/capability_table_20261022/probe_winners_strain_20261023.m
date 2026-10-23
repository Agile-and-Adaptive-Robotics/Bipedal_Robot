% Probe Winners Strain 20261023 (scratch): for every winner with
% activePct == 0, report the post-solve strain range over the graded
% domain + taut/slack/infeasible counts from its winner mat (grounded
% notes for the result rows).
here = fileparts(mfilename('fullpath'));
T = readtable(fullfile(here, 'BiPulley_survey2_perMuscle_20261023.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
for k = 1:height(T)
    act = double(T.activePct(k));
    if ~(abs(act) < 1e-12), continue; end
    mf = char(string(T.matFile(k)));
    S = load(mf);
    g = S.g; d = S.dom;
    fprintf(['%-12s meet %5.1f%% active 0%% | strain_p over domain [%+.3f, %+.3f] | ' ...
        'taut %d slack %d infeas %d | |tauB| [%6.2f, %6.2f] | tau_h peak %6.2f\n'], ...
        char(string(T.site(k))), double(T.meetPct(k)), min(g.strain(d)), ...
        max(g.strain(d)), sum(g.taut(d)), sum(g.slack(d)), sum(g.infeas(d)), ...
        min(abs(g.tauB(d))), max(abs(g.tauB(d))), max(abs(S.tauHuman(d))));
end
