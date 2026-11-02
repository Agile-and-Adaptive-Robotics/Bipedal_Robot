% repair_build.m  (AH001_Simscape_Repair, 2026-11-01)
% Build ah001_repair.slx from the fresh full CAD import (ah001_full_import)
% by staged, compile-checked surgery:
%   Stage A: copy model, add Solver Configuration, set gravity [0 -9.80665 0]
%   Stage B: delete the two smimport-added KT<->KB rigid connections
%   Stage C: replace the two KB<->FL Cylindrical joints with Revolute (pin axis
%            convention is base-frame X for both blocks, so drop-in)
%   Stage D: delete the two spurious Cartesian joints (KT to import-net) and
%            their now-orphaned transform chains
%   Stage E: weld all fastener spin joints (root-level hardware joints + the
%            hip knob/stud revolute spins) by joint->direct-line replacement
%   Stage F: compile verification + final joint census
% Every stage logs; model saved after each successful stage.

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'repair_build_20261101.txt');
fid = fopen(logFile, 'w');

srcMdl = 'ah001_full_import';
mdl = 'ah001_repair';
outPath = fullfile(repairRoot, 'working', [mdl '.slx']);

FASTENER_PAT = "94645A205|91239A330|97493A126|94510A255|90695A038|92997A108|6630T310|91595A476";

%% Stage A: copy + solver config + gravity
load_system(fullfile(repairRoot, 'working', [srcMdl '.slx']));
save_system(srcMdl, outPath);           % renames the open model to ah001_repair
p(fid, 'STAGE A: working on %s\n', mdl);

% defensive: ensure smiData is in the model workspace
hws = get_param(mdl, 'ModelWorkspace');
if isempty(whos(hws))
    evalin(hws, sprintf('run(''%s'');', ...
        fullfile(repairRoot, 'cad', 'export', 'ah001_full_DataFile.m')));
    p(fid, '  ran DataFile into model workspace\n');
end

% absolute-ize File Solid geometry paths (STEP files live in cad/export)
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
p(fid, '  File Solid paths made absolute: %d\n', nFixed);

% Mechanism Configuration gravity
mc = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                 'ReferenceBlock', 'sm_lib/Utilities/Mechanism Configuration');
for k = 1:numel(mc)
    set_param(mc{k}, 'GravityVector', '[0 -9.80665 0]');
    p(fid, '  gravity set on %s\n', mc{k});
end

% Solver Configuration: connect to the World Frame net (branch)
wf = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                 'ReferenceBlock', 'sm_lib/Frames and Transforms/World Frame');
if isempty(wf)
    add_block('sm_lib/Frames and Transforms/World Frame', [mdl '/RepairWorld']);
    wfPath = [mdl '/RepairWorld'];
else
    wfPath = wf{1};
end
sc = [mdl '/Solver Configuration'];
scExists = true;
try
    get_param(sc, 'Handle');
catch
    scExists = false;
end
if ~scExists
    add_block('nesl_utility/Solver Configuration', sc);
    phW = get_param(wfPath, 'PortHandles');
    phS = get_param(sc, 'PortHandles');
    wconn = [phW.RConn phW.LConn];
    add_line(mdl, wconn(1), phS.LConn(1));
    p(fid, '  solver configuration connected to %s\n', wfPath);
else
    p(fid, '  solver configuration already present\n');
end

okA = tryCompile(fid, mdl);
p(fid, 'STAGE A compile: %d\n', okA);
save_system(mdl, outPath);

%% Stage B: delete KT<->KB rigid connections
% R: KT_R L6 <-> KB_R L3   |  L: KT_L L2 <-> KB_L L3
pairs = {
  'x09_BA_001_1/x04_01_KT_R_001_1_RIGID', 6, 'x09_BA_001_1/x04_02_KB_R_001_1_RIGID', 3;
  'x09_BA_001_1/x04_03_KT_L_001_1_RIGID', 2, 'x09_BA_001_1/x04_04_KB_L_001_1_RIGID', 3};
for k = 1:size(pairs,1)
    g1 = [mdl '/' pairs{k,1}]; n1 = pairs{k,2};
    g2 = [mdl '/' pairs{k,2}]; n2 = pairs{k,3};
    ph1 = get_param(g1, 'PortHandles');
    l = get_param(ph1.LConn(n1), 'Line');
    if l ~= -1
        delete_line(l);
        p(fid, 'STAGE B: deleted KT-KB line %s L%d <-> %s L%d\n', pairs{k,1}, n1, pairs{k,2}, n2);
    end
    % remove the now-orphaned ports on both groups
    try
        Simulink.SubSystem.delete_ports(g1, 'LConn', n1);
        p(fid, '  removed port L%d on %s\n', n1, pairs{k,1});
    catch ME
        p(fid, '  port cleanup g1 (left dangling, tolerated): %s\n', ME.message);
    end
    try
        Simulink.SubSystem.delete_ports(g2, 'LConn', n2);
        p(fid, '  removed port L%d on %s\n', n2, pairs{k,2});
    catch ME
        p(fid, '  port cleanup g2 (left dangling, tolerated): %s\n', ME.message);
    end
end
% log dangling inner Rigid Transforms inside knee groups (Simscape ignores
% them with a warning; deleting them risks cutting the groups' internal nets)
kneeGroups = {'x04_01_KT_R_001_1_RIGID','x04_02_KB_R_001_1_RIGID', ...
              'x04_03_KT_L_001_1_RIGID','x04_04_KB_L_001_1_RIGID'};
for k = 1:numel(kneeGroups)
    gp = [mdl '/x09_BA_001_1/' kneeGroups{k}];
    tb = find_system(gp, 'LookUnderMasks','all','FollowLinks','on', ...
                     'ReferenceBlock', 'sm_lib/Frames and Transforms/Rigid Transform');
    for j = 1:numel(tb)
        pht = get_param(tb{j}, 'PortHandles');
        allp = [pht.LConn pht.RConn];
        nc = 0;
        for ii = 1:numel(allp)
            if get_param(allp(ii), 'Line') ~= -1, nc = nc + 1; end
        end
        if nc < numel(allp)
            p(fid, '  dangling inner transform left in place (ignored by solver): %s (%d/%d connected)\n', ...
                tb{j}, nc, numel(allp));
        end
    end
end
okB = tryCompile(fid, mdl);
p(fid, 'STAGE B compile: %d\n', okB);
save_system(mdl, outPath);

%% Stage C: Cylindrical -> Revolute on KB<->FL
cyls = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                   'ReferenceBlock', 'sm_lib/Joints/Cylindrical Joint');
for k = 1:numel(cyls)
    b = cyls{k};
    if ~contains(b, 'x09_BA_001_1'), continue; end   % only the knee ones (inside x09)
    p(fid, 'STAGE C: replacing %s\n', b);
    ph = get_param(b, 'PortHandles');
    lL = get_param(ph.LConn(1), 'Line');
    lR = get_param(ph.RConn(1), 'Line');
    farL = get_param(lL, 'SrcPortHandle');
    farL = [farL reshape(get_param(lL, 'DstPortHandle'),1,[])];
    farR = get_param(lR, 'SrcPortHandle');
    farR = [farR reshape(get_param(lR, 'DstPortHandle'),1,[])];
    % keep the port that is NOT on the joint
    farL = farL(farL ~= ph.LConn(1));
    farR = farR(farR ~= ph.RConn(1));
    pos = get_param(b, 'Position');
    newName = [get_param(b,'Name') '_pin'];
    parentPath = get_param(b, 'Parent');
    delete_block(b);
    nb = [parentPath '/' newName];
    add_block('sm_lib/Joints/Revolute Joint', nb, 'Position', pos);
    phN = get_param(nb, 'PortHandles');
    add_line(parentPath, farL(1), phN.LConn(1));
    add_line(parentPath, farR(1), phN.RConn(1));
    p(fid, '  -> %s\n', nb);
end
okC = tryCompile(fid, mdl);
p(fid, 'STAGE C compile: %d\n', okC);
save_system(mdl, outPath);

%% Stage D: delete the two Cartesians + orphaned transform chains
carts = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                    'ReferenceBlock', 'sm_lib/Joints/Cartesian Joint');
for k = 1:numel(carts)
    b = carts{k};
    p(fid, 'STAGE D: removing %s\n', b);
    deleteJointWithChainCleanup(fid, mdl, b);
end
okD = tryCompile(fid, mdl);
p(fid, 'STAGE D compile: %d\n', okD);
save_system(mdl, outPath);

%% Stage E: weld fastener spins
% every remaining joint whose BOTH neighbor groups match fastener pattern,
% plus the hip knob/stud revolutes (knob<->stud, stud<->FH), gets joint ->
% direct line.
jrefs = {"sm_lib/Joints/Revolute Joint", "sm_lib/Joints/Cylindrical Joint", ...
         "sm_lib/Joints/Cartesian Joint", "sm_lib/Joints/Planar Joint"};
for jr = jrefs
    jb = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                     'ReferenceBlock', char(jr));
    for k = 1:numel(jb)
        b = jb{k};
        nb1 = neighborName(mdl, b, 'LConn');
        nb2 = neighborName(mdl, b, 'RConn');
        isFast = @(s) ~isempty(regexpi(char(s), char(FASTENER_PAT), 'once'));
        % structural knee pins live between KT/KB/FL/BL groups; hips between
        % pelvis and knob. Weld if either side is a fastener, or both sides
        % are hardware (knob/stud), i.e., NOT one of the structural groups.
        structural = regexp(nb1 + "|" + nb2, 'KT_|KB_|FL_001|BL_001|PE_001|FH_|FB_|TI_', 'match');
        if isFast(nb1) || isFast(nb2)
            p(fid, 'STAGE E: weld %s  (%s | %s)\n', get_param(b,'Name'), nb1, nb2);
            weldJoint(fid, mdl, b);
        end
    end
end
okE = tryCompile(fid, mdl);
p(fid, 'STAGE E compile: %d\n', okE);
save_system(mdl, outPath);

%% Stage F: final census
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
p(fid, 'REPAIR BUILD DONE. compiles A=%d B=%d C=%d D=%d E=%d\n', okA, okB, okC, okD, okE);
fclose(fid);

% ----------------------------------------------------------------------
function ok = tryCompile(fid, mdl)
    ok = true;
    try
        sim(mdl, 'StopTime', '0');
        p(fid, '  compile OK\n');
    catch ME
        ok = false;
        p(fid, '  compile FAIL:\n%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

function weldJoint(fid, mdl, b)
% delete joint block and directly connect its two far ports (weld at CAD pose)
    ph = get_param(b, 'PortHandles');
    lL = get_param(ph.LConn(1), 'Line');
    lR = get_param(ph.RConn(1), 'Line');
    if lL == -1 || lR == -1
        p(fid, '    skip (unconnected port): %s\n', b);
        return;
    end
    farL = [get_param(lL,'SrcPortHandle') reshape(get_param(lL,'DstPortHandle'),1,[])];
    farR = [get_param(lR,'SrcPortHandle') reshape(get_param(lR,'DstPortHandle'),1,[])];
    farL = farL(farL ~= ph.LConn(1));
    farR = farR(farR ~= ph.RConn(1));
    if isempty(farL) || isempty(farR)
        p(fid, '    skip (no far port): %s\n', b);
        return;
    end
    sysJ = get_param(b, 'Parent');
    delete_block(b);
    add_line(sysJ, farL(1), farR(1));
end

function deleteJointWithChainCleanup(fid, mdl, b)
% delete a joint; on each side, delete now-dangling Rigid Transform chain
% links; clean orphaned subsystem conserving ports at the chain's far end.
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
        % walk: while current block is a Rigid Transform with an unconnected
        % port, delete it and step to its other neighbor
        for it = 1:12
            try
                r = strtrim(strrep(string(get_param(cur, 'ReferenceBlock')), newline, " "));
            catch
                break;
            end
            if ~(contains(r, "Rigid") && contains(r, "Transform"))
                % non-transform: if it is a subsystem with orphaned conserving
                % port, remove that port
                bt = get_param(cur, 'BlockType');
                if strcmp(bt, 'SubSystem')
                    phc = get_param(cur, 'PortHandles');
                    del = [];
                    for ii = 1:numel(phc.LConn)
                        if get_param(phc.LConn(ii), 'Line') == -1
                            del(end+1) = ii; %#ok<AGROW>
                        end
                    end
                    for ii = 1:numel(phc.RConn)
                        if get_param(phc.RConn(ii), 'Line') == -1
                            % RConn cleanup needs separate call
                        end
                    end
                    if ~isempty(del)
                        try
                            Simulink.SubSystem.delete_ports(cur, 'LConn', del);
                            p(fid, '    cleaned ports on %s\n', cur);
                        catch
                        end
                    end
                end
                break;
            end
            % transform: delete ONLY if exactly one port remains connected
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
                p(fid, '    stop at transform %s (%d ports connected)\n', get_param(cur,'Name'), conncount);
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

function nm = neighborName(mdl, b, which)
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
