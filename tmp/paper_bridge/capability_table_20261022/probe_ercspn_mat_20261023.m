% Probe Ercspn Mat 20261023 (scratch): characterize the winner mat grids.
S = load('D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\BiPulley_survey2_ercspn_r_nB1_G1_20261023_051731.mat');
fprintf('mat vars: %s\n', strjoin(fieldnames(S), ', '));
g = S.g; d = S.dom;
fprintf('size(dom) = %s, class(dom) = %s, sum(dom) = %d\n', ...
    mat2str(size(d)), class(d), sum(d(:)));
fprintf('size(g.meet) = %s, class = %s\n', mat2str(size(g.meet)), class(g.meet));
fprintf('size(g.active) = %s\n', mat2str(size(g.active)));
fprintf('size(g.strain) = %s\n', mat2str(size(g.strain)));
fprintf('meet: %d / %d\n', sum(g.meet(d)), sum(d(:)));
fprintf('active: %d / %d\n', sum(g.active(d)), sum(d(:)));
fprintf('strain range over domain: [%.4f, %.4f]\n', ...
    min(g.strain(d)), max(g.strain(d)));
fprintf('taut: %d, slack: %d, infeas: %d\n', sum(g.taut(d)), ...
    sum(g.slack(d)), sum(g.infeas(d)));
fprintf('|tauB| range: [%.2f, %.2f]; |tau_h| peak: %.2f\n', ...
    min(abs(g.tauB(d))), max(abs(g.tauB(d))), max(abs(S.tauHuman(d))));
fprintf('surveyMeta: meetPct %.2f activePct %.2f\n', ...
    S.surveyMeta.meetPct, S.surveyMeta.activePct);
