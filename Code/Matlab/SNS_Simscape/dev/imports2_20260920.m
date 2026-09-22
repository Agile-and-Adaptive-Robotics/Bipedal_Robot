%% imports2_20260920.m — (A) REAL 09_BA_003 Multibody-Link XML import + joint
%% inventory (the existing mdl_knee_rig_import_tmp_imported.slx turned out to
%% be the 1-link sw2urdf stub: 12 blocks, 0 joints); (B) patch the imported
%% OpenSim model's 19 Visual File Solids to absolute mesh paths, then compile,
%% 0.01 s sim, joint inventory, save. Rig import runs FIRST; the osim compile
%% LAST (a model whose update fails leaves a GUI-tree state that crashes the
%% process at teardown — 2026-09-20; everything before it must already be
%% saved/logged).

sns = fileparts(fileparts(mfilename('fullpath')));   % dev -> SNS_Simscape
addpath(sns);
fprintf('=== MATLAB %s ===\n', version);

%% --- A. 09_BA_003.xml -> smimport -> mdl_knee_rig_xml ----------------------
try
    xml = fullfile(sns, '..', '..', '..', 'Solid_Models', ...
        'Biomimetics_2022-Knee_Test', 'Knee assembly', '09_BA_003.xml');
    xml = fullfile(xml);  % normalize
    assert(isfile(xml), 'XML not found: %s', xml);
    srcDir = fileparts(xml);
    tmp = fullfile(srcDir, 'mdl_knee_rig_xml.xml');
    if bdIsLoaded('mdl_knee_rig_xml'), close_system('mdl_knee_rig_xml', 0); end
    copyfile(xml, tmp, 'f');
    clob = onCleanup(@() delete(tmp));
    smimport(tmp);
    mdl = 'mdl_knee_rig_xml';
    assert(bdIsLoaded(mdl), 'smimport did not create mdl_knee_rig_xml');
    blks = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
    fprintf('\nRIG XML IMPORT: %d blocks\n', numel(blks));
    nJ = 0; nSol = 0;
    for k = 1:numel(blks)
        rb = '';
        try, rb = get_param(blks{k}, 'ReferenceBlock'); if ~ischar(rb), rb = ''; end, catch, end
        if ~isempty(strfind(rb, 'sm_lib/Joints/'))
            nJ = nJ + 1;
            fprintf('  JOINT  %-44s <- %s\n', strrep(blks{k}, [mdl '/'], ''), rb);
        elseif ~isempty(strfind(rb, 'sm_lib/Body/Elements/'))
            nSol = nSol + 1;
        end
    end
    fprintf('rig XML: %d joint blocks, %d solid elements\n', nJ, nSol);
    % gravity: CAD Z is up
    try
        mc = find_system(mdl, 'LookUnderMasks', 'all', 'MaskType', 'Mechanism Configuration');
        if ~isempty(mc)
            set_param(mc{1}, 'GravityVector', '[0 0 -9.80665]');
            fprintf('gravity -Z set\n');
        end
    catch ME
        fprintf('gravity tweak skipped: %s\n', ME.message);
    end
    save_system(mdl, fullfile(sns, 'mdl_knee_rig_xml_imported.slx'));
    close_system(mdl, 0);
    fprintf('RIG XML IMPORT SAVED\n');
catch ME
    fprintf('RIG XML IMPORT FAILED: %s\n', ME.message);
end

%% --- B. osim model: absolute mesh paths + compile + sim --------------------
try
    mdl = 'Gait2392_simbody_simscape';
    load_system(fullfile(sns, 'osim_import', [mdl '.slx']));
    geoDir = fullfile(sns, 'osim_import', 'Geometry');
    vis = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', ...
        'MaskType', 'File Solid');
    fprintf('\nOSIM PATCH: %d File Solids\n', numel(vis));
    nPatched = 0;
    for k = 1:numel(vis)
        fn = '';
        try, fn = get_param(vis{k}, 'FileName'); catch, end
        [~, nm, xt] = fileparts(fn);
        absF = fullfile(geoDir, [nm xt]);
        if isfile(absF)
            set_param(vis{k}, 'FileName', absF);
            nPatched = nPatched + 1;
        else
            fprintf('  could not resolve %s (param was "%s")\n', vis{k}, fn);
        end
    end
    fprintf('patched %d/%d File Solids to %s\n', nPatched, numel(vis), geoDir);
    set_param(mdl, 'StopTime', '0.01');
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('OSIM COMPILE OK\n');
    nJ = 0;
    jt = {};
    blks = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
    for k = 1:numel(blks)
        rb = '';
        try, rb = get_param(blks{k}, 'ReferenceBlock'); if ~ischar(rb), rb = ''; end, catch, end
        if ~isempty(strfind(rb, 'sm_lib/Joints/'))
            nJ = nJ + 1; jt{end+1} = rb; %#ok<SAGROW>
        end
    end
    [uJ, ~, iu] = unique(jt);
    for u = 1:numel(uJ), fprintf('  %-42s x%d\n', uJ{u}, sum(iu == u)); end
    fprintf('osim model: %d blocks, %d joint refs\n', numel(blks), nJ);
    sim(mdl);
    fprintf('OSIM 0.01 s SIM OK\n');
    save_system(mdl);
    close_system(mdl, 0);
    fprintf('OSIM MODEL SAVED (mesh paths absolute)\n');
catch ME
    fprintf('OSIM PATCH/COMPILE FAILED: %s\n', ME.message);
    for c = 1:numel(ME.cause)
        m = ME.cause{c}.message;
        fprintf('CAUSE %d: %s\n', c, m(1:min(end, 140)));
    end
    % deliberately no further model interaction (teardown-crash hazard)
end
fprintf('=== ALL DONE ===\n');
