% repair_foot.m  (AH001_Simscape_Repair, 2026-11-01, Ben's direction)
% Turn the raw Foot Ankle import into Ben's articulated foot spec, twice
% (foot_L, foot_R), each with EXACTLY 3 DOF:
%   subtalar-calcaneus 2-DOF: Revolute (Ankle Interior <-> Bottom Joint) [native]
%                            + Planar2 (Interior <-> Exterior) -> Revolute
%   MTP 1-DOF:                Revolute2 (main plate <-> foot_axle_1) [native]
%                            + Cylindrical (front foot <-> foot_axle_1) -> Revolute
%   welded: mid breaks (Planar, Planar1, Planar3, Cylindrical1, Cylindrical2,
%           Revolute1, Revolute3), all five toe revolutes (Revolute4..8)
%   deleted: the two World anchor welds (front foot F1, mid foot F2), the
%            redundant Parallel constraint (Exterior <-> Bottom)
% Then expose the ankle-housing datum port (branch Ankle_Exterior group).

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'repair_foot_20261101.txt');
fid = fopen(logFile, 'w');
expDir = fullfile(repairRoot, 'cad', 'export_foot');

WELD_LIST  = {'Planar', 'Planar1', 'Planar3', 'Cylindrical1', 'Cylindrical2', ...
              'Revolute1', 'Revolute3', 'Revolute4', 'Revolute5', ...
              'Revolute6', 'Revolute7', 'Revolute8'};
TO_REV     = {'Planar2', 'Cylindrical'};    % 2nd ankle axis, MTP second half
TO_DELETE  = {'Parallel'};                  % redundant constraint
GROUNDPORTS = {'New_Front_foot_1_RIGID', 1, 'Left_mid_foot_1_RIGID', 1};  % group, RConn idx (probed live)

for mdlBase = ["foot_L", "foot_R"]
    mdl = char(mdlBase);
    % fresh import each time from the same XML
    cd(expDir);
    smimport('foot_ankle.xml', 'ModelName', mdl);
    p(fid, '\n=== %s ===\n', mdl);

    load_system(mdl); %#ok<LLOAD>  (smimport leaves it loaded)
    % absolute STEP paths
    allB = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
    for k = 1:numel(allB)
        try
            r = strtrim(strrep(string(get_param(allB{k}, 'ReferenceBlock')), newline, " "));
        catch
            continue;
        end
        if r ~= "sm_lib/Body Elements/File Solid", continue; end
        fn = get_param(allB{k}, 'ExtGeomFileName');
        if ~isempty(fn) && isempty(regexp(fn, '^[A-Za-z]:', 'once'))
            set_param(allB{k}, 'ExtGeomFileName', fullfile(expDir, fn));
        end
    end

    % delete World anchor welds: find lines from World to the two groups
    wf = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                     'ReferenceBlock', 'sm_lib/Frames and Transforms/World Frame');
    if ~isempty(wf)
        phW = get_param(wf{1}, 'PortHandles');
        for pp = [phW.LConn phW.RConn]
            l = get_param(pp, 'Line');
            if l == -1, continue; end
            x = [get_param(l,'SrcBlockHandle') reshape(get_param(l,'DstBlockHandle'),1,[])];
            for xx = x
                if xx == -1, continue; end
                nm = get_param(xx, 'Name');
                if strcmp(nm, 'New_Front_foot_1_RIGID') || strcmp(nm, 'Left_mid_foot_1_RIGID')
                    delete_line(l);
                    p(fid, '  deleted World anchor line to %s\n', nm);
                    break;
                end
            end
        end
    end

    % joint type swaps
    load_system('sm_lib');
    RT_PATH = '';
    tb0 = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                      'ReferenceBlock', 'sm_lib/Frames and Transforms/Rigid Transform');
    if ~isempty(tb0), RT_PATH = get_param(tb0{1}, 'ReferenceBlock'); end

    for nm = WELD_LIST
        b = [mdl '/' nm];
        try
            replace_block(mdl, 'Name', nm{1}, RT_PATH, 'noprompt');
            set_param(b, 'Name', [nm{1} '_weld']);
            p(fid, '  welded %s\n', nm{1});
        catch ME
            p(fid, '  WELD FAIL %s: %s\n', nm{1}, ME.message);
        end
    end
    for nm = TO_REV
        try
            replace_block(mdl, 'Name', nm{1}, 'sm_lib/Joints/Revolute Joint', 'noprompt');
            set_param([mdl '/' nm{1}], 'Name', [nm{1} '_rev']);
            p(fid, '  -> revolute %s\n', nm{1});
        catch ME
            p(fid, '  REV FAIL %s: %s\n', nm{1}, ME.message);
        end
    end
    for nm = TO_DELETE
        try
            delete_block([mdl '/' nm{1}]);
            p(fid, '  deleted %s\n', nm{1});
        catch ME
            p(fid, '  DEL FAIL %s: %s\n', nm{1}, ME.message);
        end
    end

    % expose the ankle-housing datum: branch Ankle_Exterior group's Solid base net
    extGrp = [mdl '/Ankle_Exterior_1_RIGID'];
    sol = find_system(extGrp, 'LookUnderMasks','all','FollowLinks','on', ...
                      'ReferenceBlock', 'sm_lib/Body Elements/File Solid');
    if ~isempty(sol)
        phS = get_param(sol{1}, 'PortHandles');
        cands = [phS.LConn phS.RConn];
        for ii = 1:numel(cands)
            if cands(ii) ~= -1 && get_param(cands(ii), 'Line') ~= -1
                cpPath = [extGrp '/housing_datum_port'];
                add_block('nesl_utility/Connection Port', cpPath);
                for pv = {'PortSide','right'; 'ConnectionType','physical'}
                    try, set_param(cpPath, pv{1}, pv{2}); catch, end
                end
                phC = get_param(cpPath, 'PortHandles');
                allc = [phC.LConn phC.RConn];
                add_line(extGrp, cands(ii), allc(1));
                p(fid, '  housing datum port branched\n');
                break;
            end
        end
    end

    % ground the housing to World for a standalone compile (cut when mounted)
    try
        wf2 = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on', ...
                          'ReferenceBlock', 'sm_lib/Frames and Transforms/World Frame');
        phG = get_param(extGrp, 'PortHandles');
        newPort = -1;
        for pp = [phG.LConn phG.RConn]
            if pp ~= -1 && get_param(pp, 'Line') == -1
                newPort = pp; break;
            end
        end
        phW2 = get_param(wf2{1}, 'PortHandles');
        wc = [phW2.LConn phW2.RConn];
        add_line(mdl, newPort, wc(1));
        p(fid, '  housing grounded to World for standalone compile (cut when mounted)\n');
    catch ME
        p(fid, '  world ground note: %s\n', ME.message);
    end

    ok = true;
    try
        sim(mdl, 'StopTime', '0');
        p(fid, '  compile OK\n');
    catch ME
        ok = false;
        p(fid, '  compile FAIL: %s\n', ME.message);
        try
            cs = ME.cause; if isa(cs,'MException'), cs = {cs}; end
            for k2 = 1:numel(cs), p(fid, '    cause: %s\n', cs{k2}.message); end
        catch
        end
    end

    save_system(mdl, fullfile(repairRoot, 'working', [mdl '.slx']));
    close_system(mdl, 0);
    p(fid, '  %s saved (compile=%d)\n', mdl, ok);
end

cd(repairRoot);
p(fid, 'FOOT REPAIR DONE\n');
fclose(fid);

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
