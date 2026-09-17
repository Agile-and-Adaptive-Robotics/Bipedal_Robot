function e0_multibody_license()
% E0: can THIS machine build + RUN hand-built Simscape Multibody models
% (primitive blocks, no CAD import)? Builds World Frame -> Revolute Joint ->
% Cylindrical Solid, sims 0.1 s. Verdict printed + any license error shown.
% Also reports whether a generic Force Element exists (needed to later
% implement a muscle pulling between two frames - Ben's elbow idea).

try
    load_system('sm_lib');
catch e
    fprintf(['E0 VERDICT: FAIL - sm_lib will not even load: %s\n'], e.message);
    return
end
% FollowLinks on: sm_lib/Utilities is a LINKED sublibrary, and the
% Mechanism Configuration block lives behind it (verified 2026-09-16 from
% an smimport-produced model: ReferenceBlock = sm_lib/Utilities/...
% "Mechanism Configuration"). Simscape MULTIBODY models take NO Solver
% Configuration block (that is foundation-Simscape only) - the imported
% model carries none either.
% NOTE: search-criterion option pairs MUST precede the 'Type','Block'
% parameter-value pair in find_system, or the search silently returns 0.
blks = find_system('sm_lib', 'LookUnderMasks', 'all', 'FollowLinks', 'on', ...
    'SearchDepth', 5, 'Type', 'Block');

hits = @(pat) blks(~cellfun(@isempty, regexp(blks, pat, 'once')));
wf  = hits('World Frame');
rj  = hits('Joints[/\\]Revolute Joint');
cyl = hits('Cylindrical Solid');
fe  = hits('Force Element');
% Simscape Multibody DOES require a Solver Configuration block (E0 sim error
% 2026-09-16: "Each physical network must be connected to exactly one Solver
% Configuration block") even though smimport's own output carries none; it
% lives in nesl_utility and resolves by direct path.
load_system('nesl_utility');
scPath = 'nesl_utility/Solver Configuration';
sc = {scPath};
% Mechanism Configuration resolves by DIRECT PATH but find_system cannot
% enumerate it even with FollowLinks (the Utilities subsystem is protected).
% Verified 2026-09-16 from an smimport model's ReferenceBlock.
mcPath = 'sm_lib/Utilities/Mechanism Configuration';
mc = {mcPath};
try
    get_param(mcPath, 'MaskType');
catch
    mc = {};
end
if ~isempty(fe)
    fprintf('library has Force Element: %s\n', fe{1});
else
    fprintf('no generic Force Element found in top-5 library levels\n');
end
if isempty(wf) || isempty(rj) || isempty(cyl) || isempty(mc)
    fprintf('E0 VERDICT: INCONCLUSIVE - block paths not found:\n');
    fprintf('  world=%d revolute=%d cylinder=%d mechcfg=%d\n', ...
            numel(wf), numel(rj), numel(cyl), numel(mc));
    bdclose('sm_lib');
    return
end

mdl = 'sm_license_probe';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
add_block(wf{1}, [mdl '/WF'], 'Position', [100 100 130 130]);
add_block(rj{1}, [mdl '/RJ'], 'Position', [200 100 250 150]);
add_block(cyl{1}, [mdl '/Body'], 'Position', [320 100 390 160]);
add_block(mc{1}, [mdl '/MC'], 'Position', [220 220 250 250]);
add_block(sc{1}, [mdl '/SC'], 'Position', [340 220 370 250]);
% conserving connections: Simscape frame ports are LConn/RConn port handles,
% NOT Simulink 'blk/1' ports. World -> Joint base(B), Joint follower(F) -> Solid.
% Stage-labeled errors so a failure pinpoints the exact step.
ok = false; err = ''; stage = 'port handles';
try
    wf = get_param([mdl '/WF'], 'PortHandles');
    mc = get_param([mdl '/MC'], 'PortHandles');
    scP = get_param([mdl '/SC'], 'PortHandles');
    rj = get_param([mdl '/RJ'], 'PortHandles');
    bd = get_param([mdl '/Body'], 'PortHandles');
    fprintf('  conn ports: WF %d, MC %d, RJ %d, Body %d\n', ...
        numel([wf.RConn wf.LConn]), numel([mc.RConn mc.LConn]), ...
        numel([rj.RConn rj.LConn]), numel([bd.RConn bd.LConn]));
    pWorld = pickConn(wf);
    pMC    = pickConn(mc);
    pSC    = pickConn(scP);
    pSolid = pickConn(bd);
    cj = [rj.LConn rj.RConn];               % joint has TWO frame ports (B, F)
    assert(numel(cj) >= 2, 'revolute joint has fewer than 2 frame ports');
    pB = cj(1);
    pF = cj(end);

    stage = 'world->joint line';
    add_line(mdl, pWorld, pB, 'autorouting', 'on');
    stage = 'mech-config branch';
    add_line(mdl, pMC, pWorld, 'autorouting', 'on');   % branch off the net
    stage = 'solver-config branch';
    add_line(mdl, pSC, pWorld, 'autorouting', 'on');   % branch off the net
    stage = 'joint->solid line';
    add_line(mdl, pF, pSolid, 'autorouting', 'on');
    try
        % gravity param is 'GravityVector' on this block (verified from an
        % smimport model); default is already [0 0 -9.80665]
        set_param([mdl '/MC'], 'GravityVector', '[0 0 -9.80665]');
    catch
    end
    set_param(mdl, 'StopTime', '0.1');
    stage = 'sim';
    sim(mdl);
    ok = true;
catch e
    err = sprintf('at stage "%s": %s', stage, e.message);
end
if ok
    fprintf(['E0 VERDICT: PASS - hand-built Simscape Multibody RUNS on this ' ...
             'machine (World+Revolute+Cylinder, 0.1 s clean)\n']);
else
    fprintf('E0 VERDICT: FAIL - sim error:\n%s\n', err);
end
try
    close_system(mdl, 0);
catch
end
bdclose('sm_lib');
bdclose('nesl_utility');
end

function p = pickConn(ph)
% first available conserving/frame port handle (RConn preferred, then LConn)
c = [ph.RConn ph.LConn];
assert(~isempty(c), 'block has no conserving ports');
p = c(1);
end
