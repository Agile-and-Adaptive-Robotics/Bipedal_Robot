%Dig_rehydrate_fronts.m — one-shot: rebuild the full run state for the existing
%extensor front captures (ctx, humanBest, predBest, cBest, routeInfo) and re-save
%each capture as a FULL workspace, so the display sections of Opt_run_Ext.m run
%directly after loading them (Ben, 2026-09-10).
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization');

pairs = { ...
 'minimizeExt10mmX3_results_20260910_noT3.mat',          '2trans'; ...
 'minimizeExt10mmX3_results_20260909_1trans_Lbolt.mat',  '1trans'};

for c = 1:size(pairs,1)
    fn = pairs{c,1}; conv = pairs{c,2};
    S = load(fn, 'sol_actual', 'g', 'filtered_results', 'xCols', 'results_cv', 'a0', 'a1', 'a2');
    sa = S.sol_actual;

    ctx = buildKneeExtContext20mm();                    %deterministic rebuild
    ctx.Xi0 = sa(1); ctx.Xi1 = sa(2); ctx.Xi2 = sa(3); ctx.Xi3 = sa(4);

    predBest = predictKneeExt20mm(sa(1:3), ctx);        %same call the run used
    cBest     = nonlconExt20mm(sa(1:3), ctx);
    humanBest = interp1(ctx.humanAngleD, ctx.humanTorque, ctx.phiD, 'pchip');

    xBest = sa(1:3); xSeed = sa(1:3); fBest = NaN; fRefined = NaN;
    f0 = objective_KneeExt20mm(xSeed, ctx);

    assignin('base', 'ctx', ctx); %#ok<NASGU>
    display('Rehydrated: '); disp(fn);

    % Re-save the FULL rehydrated workspace under the same convention name.
    evalin('base', sprintf('save(''%s'')', fn));
end
fprintf('All captures rehydrated and re-saved (full workspaces).\n');
