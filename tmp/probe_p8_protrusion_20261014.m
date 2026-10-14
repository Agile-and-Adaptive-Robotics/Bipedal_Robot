%% Probe (2026-10-14): p8 slight protrusion past the pEnd -> p7 chord
%% for the run-8 design. Sweeps perpendicular protrusion eps, reports
%% whether the release logic contacts p8, the added path length, and
%% the no-stretch (cRestLength) impact + the tendon trim that offsets
%% it. No constraints -- just the CAD numbers Ben asked about (1-2 mm
%% of path length).
root = 'D:/GitHub/Bipedal_Robot';
moDir = fullfile(root,'Code','Matlab','Mesh_Optimization');
addpath(genpath(fullfile(root,'Code','Matlab')));
addpath(moDir);
addpath(fullfile(root,'Testing_Data','2022_02_Festo'),'-end');

ctx = buildKneeExtContext20mm();
S = load(fullfile(moDir,'Results','Vas_Pam_20mm_Result.mat'),'xBest');
xBest = reshape(S.xBest,1,[]);
p1 = xBest(1:3); pEnd = xBest(4:6); rest = xBest(7);

% run-8 winner geometry
tendonR8 = 0.090440;
p6R8 = [0.116323, 0.015040];
p7R8 = [0.131323, -0.018974];

% chord pEnd -> p7 (t1 frame; all t1-fixed so this is pose-independent)
u = (p7R8 - pEnd(1:2)); u = u/norm(u);
p8stock = [0.06836, -0.01599];
sFoot = dot(p8stock - pEnd(1:2), u);
sFoot = min(max(sFoot, 0), norm(p7R8 - pEnd(1:2)));
foot = pEnd(1:2) + sFoot*u;
n = [u(2), -u(1)];   % unit normal giving offChord = -eps (stock side)

fprintf('run-8 chord: pEnd->p7 length %.1f mm; foot of stock p8 at [%.1f %.1f] t1 mm;\n', ...
    1000*norm(p7R8 - pEnd(1:2)), 1000*foot);
fprintf('stock p8 sits %.1f mm from that chord (inside, tibia side)\n\n', ...
    1000*norm(p8stock - foot));

epsList = [0, 2, 5, 10, 14];   % mm of protrusion beyond the chord
baseMaxLmt = NaN;
fprintf('eps_mm  p8_contacted  release_deg  extra_path_mm  cRest_mm  cRest_trimmed_mm\n');
for kE = 1:numel(epsList)
    p8e = foot + (epsList(kE)/1000)*n;
    ctxAdj = ctx;
    ctxAdj.geo.seedP6XY = p6R8;
    ctxAdj.geo.seedP7XY = p7R8;
    ctxAdj.geo.seedP8XY = p8e;
    ctxAdj.geo.tibiaWrapCircumcircle = true;

    pred = predictKneeExt20mm([p1, pEnd, rest, tendonR8], ctxAdj);

    if kE == 1
        baseMaxLmt = pred.maxPathLength0;
    end
    extra = 1000*(pred.maxPathLength0 - baseMaxLmt);
    cRest = 1000*pred.cRestLength;
    % tendon trimmed by the added length restores the eps=0 budget
    predTrim = predictKneeExt20mm([p1, pEnd, rest, tendonR8 - extra/1000], ctxAdj);
    cRestTrim = 1000*predTrim.cRestLength;

    if pred.routeInfo.active(8,1)
        rel = -120.0;   % contacted from the start
    else
        rel = pred.routeInfo.eliminatedAngleD(8);
        if isnan(rel), rel = 999; end   % never contacted at all
    end

    fprintf('%5.1f  %8d  %12.2f  %12.2f  %8.2f  %8.2f\n', ...
        epsList(kE), pred.routeInfo.active(8,1), rel, extra, cRest, cRestTrim);
end

disp('PROBE DONE');
