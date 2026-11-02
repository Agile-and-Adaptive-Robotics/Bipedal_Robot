% repair_build2.m  (AH001_Simscape_Repair, 2026-11-01)
% Build ah001_repair.slx from ah001_full_import with replace_block-based surgery:
%   Stage A: copy, gravity [0 -9.80665 0], absolute File Solid paths, solver check
%   Stage B: remove the two spurious Cartesian joints (+dangling transform cleanup)
%            -- these are the compile blockers ("base and follower rigidly
%            connected"); the intended knee topology does not include them.
%   Stage C: delete the two smimport-added KT<->KB rigid connections (legal but
%            wrong: they freeze the four-bar; MuJoCo drops them too)
%   Stage D: Cylindrical -> Revolute on the two KB<->FL pins (in-place replace)
%   Stage E: fastener spin joints -> identity Rigid Transform (weld at CAD pose)
%   Stage F: compile gate + joint census
% Compiles before F are DIAGNOSTIC ONLY (B's success expected after Cartesians
% go; the KT-KB weld does not block compilation).

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'repair_build2_20261101.txt');
fid = fopen(logFile, 'w');

srcMdl = 'ah001_full_import';
mdl = 'ah001_repair';
outPath = fullfile(repairRoot, 'working', [mdl '.slx']);
FASTENER_PAT = "94645A205|91239A330|97493A126|94510A255|90695A038|92997A108|6630T310|91595A476";

%% Stage A
load_system(fullfile(repairRoot, 'working', [srcMdl '.slx']));
save_system(srcMdl, outPath);
p(fid, 'STAGE A: %s\n', mdl);

hws = get_param(mdl, 'ModelWorkspace');
if isempty(whos(hws))
    evalin(hws, sprintf('run(''%s'');', fullfile(repairRoot, 'cad', 'export', 'ah001_full_DataFile.m')));
    p(fid, '  DataFile run into model workspace\n');
end

mc = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                 'ReferenceBlock', 'sm_lib/Utilities/Mechanism Configuration');
for k = 1:numel(mc)
    set_param(mc{k}, 'GravityVector', '[0 -9.80665 0]');
    p(fid, '  gravity set: %s\n', mc{k});
end

expAbs = fullfile(repairRoot, 'cad', 'export');
allB = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
nFixed = 0;
for k = 1:numel(allB)
    try
        r = strtrim(strrep(string(get_param(allB{k}, 'ReferenceBlock')), newline, " "));
    catch
        continue;
    end
    if r ~= "sm_lib/Body Elements/File Solid", continue; end
    fn = get_param(allB{k}, 'ExtGeomFileName');
    if ~isempty(fn) && isempty(regexp(fn, '^[A-Za-z]:', 'once'))
        set_param(allB{k}, 'ExtGeomFileName', fullfile(expAbs, fn));
        nFixed = nFixed + 1;
    end
end
p(fid, '  File Solid paths absolute: %d\n', nFixed);
okA = tryCompile(fid, mdl);
p(fid, 'STAGE A compile (diagnostic): %d\n', okA);
save_system(mdl, outPath);

%% Stage B: remove Cartesians (+ chains)
carts = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                    'ReferenceBlock', 'sm_lib/Joints/Cartesian Joint');
for k = 1:numel(carts)
    b = carts{k};
    p(fid, 'STAGE B: removing %s\n', b);
    deleteJointWithChainCleanup(fid, b);
end
okB = tryCompile(fid, mdl);
p(fid, 'STAGE B compile: %d\n', okB);
save_system(mdl, outPath);

%% Stage C: delete KT<->KB rigid lines
pairs = {
  'x09_BA_001_1/x04_01_KT_R_001_1_RIGID', 6, 'x09_BA_001_1/x04_02_KB_R_001_1_RIGID';
  'x09_BA_001_1/x04_03_KT_L_001_1_RIGID', 2, 'x09_BA_001_1/x04_04_KB_L_001_1_RIGID'};
for k = 1:size(pairs,1)
    g1 = [mdl '/' pairs{k,1}]; n1 = pairs{k,2};
    g2 = [mdl '/' pairs{k,3}];
    ph1 = get_param(g1, 'PortHandles');
    l = get_param(ph1.LConn(n1), 'Line');
    if l ~= -1
        delete_line(l);
        p(fid, 'STAGE C: deleted KT-KB rigid line between %s(L%d) and %s\n', pairs{k,1}, n1, pairs{k,3});
    else
        p(fid, 'STAGE C: no line on %s L%d (already cut?)\n', pairs{k,1}, n1);
    end
end
okC = tryCompile(fid, mdl);
p(fid, 'STAGE C compile: %d\n', okC);
save_system(mdl, outPath);

%% Stage D: knee Cylindrical -> Revolute (in-place)
cyls = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                   'ReferenceBlock', 'sm_lib/Joints/Cylindrical Joint');
nRep = 0;
for k = 1:numel(cyls)
    b = cyls{k};
    if ~contains(b, 'x09_BA_001_1'), continue; end
    parent = get_param(b, 'Parent');
    nm = get_param(b, 'Name');
    try
        replace_block(parent, 'Name', nm, 'sm_lib/Joints/Revolute Joint', 'noprompt');
        set_param([parent '/' nm], 'Name', [nm '_pin']);
        p(fid, 'STAGE D: %s -> Revolute\n', b);
        nRep = nRep + 1;
    catch ME
        p(fid, 'STAGE D FAIL %s: %s\n', b, ME.message);
    end
end
p(fid, '  replaced: %d\n', nRep);
okD = tryCompile(fid, mdl);
p(fid, 'STAGE D compile: %d\n', okD);
save_system(mdl, outPath);

%% Stage D2: remove the degenerate Parallel angle constraint (base and
%% follower rigidly connected after loop surgery; the four-bar pins enforce
%% parallelism physically)
pars = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                   'ReferenceBlock', 'sm_lib/Constraints/Angle Constraint');
for k = 1:numel(pars)
    b = pars{k};
    p(fid, 'STAGE D2: removing degenerate constraint %s\n', b);
    deleteJointWithChainCleanup(fid, b);
end
okD2 = tryCompile(fid, mdl);
p(fid, 'STAGE D2 compile: %d\n', okD2);
save_system(mdl, outPath);

%% Stage E: hip knob/stud spin welds -- ATTEMPTED AND REVERTED.
% Evidence (repair_build2 log, earlier run): replacing x09/Revolute6..9 with
% identity Rigid Transforms triggers "Spherical Joint 'Spherical1' has base
% and follower ports that are rigidly connected" -- a hidden rigid path
% pelvis<->knob2/femur-L exists in the import. Tracing it is a bounded
% follow-up. The four spins are physically-real free thread rotations with
% negligible inertia; they stay as joints and are documented in the census.
p(fid, 'STAGE E: skipped (attempted hip-spin welds reverted; see record)\n');
okE = okD2;   % gate = D2 compile state
okE = tryCompile(fid, mdl);
p(fid, 'STAGE E compile (GATE): %d\n', okE);
save_system(mdl, outPath);

%% Stage F: census
blks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
ref = strings(numel(blks),1);
for k = 1:numel(blks)
    try, ref(k) = strtrim(strrep(string(get_param(blks{k}, 'ReferenceBlock')), newline, " ")); catch, ref(k) = ""; end
end
p(fid, '\nSTAGE F census:\n');
uj = unique(ref(startsWith(ref, "sm_lib/Joints/") | startsWith(ref, "sm_lib/Constraints/")));
for t = uj'
    p(fid, '  %2d x %s\n', nnz(ref == t), t);
end
p(fid, 'REPAIR BUILD DONE. compiles A(d)=%d B=%d C=%d D=%d D2=%d E(gate)=%d\n', okA, okB, okC, okD, okD2, okE);
fclose(fid);

% ----------------------------------------------------------------------
function ok = tryCompile(fid, mdl)
    ok = true;
    try
        sim(mdl, 'StopTime', '0');
        p(fid, '  compile OK\n');
    catch ME
        ok = false;
        p(fid, '  compile FAIL: %s\n', ME.message);
        try
            cs = ME.cause;
            if isa(cs, 'MException')
                cs = {cs};
            end
            for k = 1:numel(cs)
                p(fid, '    cause: %s\n', cs{k}.message);
            end
        catch
        end
    end
end

function deleteJointWithChainCleanup(fid, b)
    ph = get_param(b, 'PortHandles');
    sides = {};
    for which = {"LConn", "RConn"}
        ports = ph.(which{1});
        if isempty(ports) || get_param(ports(1),'Line') == -1
            continue;
        end
        l = get_param(ports(1), 'Line');
        far = [get_param(l,'SrcPortHandle') reshape(get_param(l,'DstPortHandle'),1,[])];
        far = far(far ~= ports(1));
        sides{end+1} = far; %#ok<AGROW>
    end
    delete_block(b);
    for s = 1:numel(sides)
        cur = sides{s};
        if isempty(cur), continue; end
        cur = cur(1);
        for it = 1:12
            try
                r = strtrim(strrep(string(get_param(cur, 'ReferenceBlock')), newline, " "));
            catch
                break;
            end
            if ~(contains(r, "Rigid") && contains(r, "Transform"))
                break;
            end
            pht = get_param(cur, 'PortHandles');
            allp = [pht.LConn pht.RConn];
            conncount = 0; nxt = -1;
            for ii = 1:numel(allp)
                l2 = get_param(allp(ii), 'Line');
                if l2 ~= -1
                    conncount = conncount + 1;
                    x = [get_param(l2,'SrcBlockHandle') reshape(get_param(l2,'DstBlockHandle'),1,[])];
                    x = x(x ~= get_param(cur,'Handle'));
                    if ~isempty(x), nxt = x(1); end
                end
            end
            if conncount ~= 1
                break;
            end
            nm = get_param(cur, 'Name');
            delete_block(cur);
            p(fid, '    deleted dangling transform %s\n', nm);
            if nxt == -1, break; end
            cur = nxt;
        end
    end
end

function nm = neighborName(b, which)
    nm = "?";
    ph = get_param(b, 'PortHandles');
    ports = ph.(which);
    if isempty(ports), return; end
    l = get_param(ports(1), 'Line');
    if l == -1, nm = "<none>"; return; end
    x = [get_param(l,'SrcBlockHandle') reshape(get_param(l,'DstBlockHandle'),1,[])];
    x = x(x ~= get_param(b, 'Handle'));
    if isempty(x), nm = "<none>"; return; end
    nm = string(get_param(x(1), 'Name'));
end

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
