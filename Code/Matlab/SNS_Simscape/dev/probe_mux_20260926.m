function probe_mux_20260926()
sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));
fprintf('EXT_Mux Inputs = %s\n', get_param([mdl '/EXT_Mux'], 'Inputs'));
ph = get_param([mdl '/EXT_Mux'], 'PortHandles');
fprintf('EXT_Mux in=%d out=%d\n', numel(ph.Inport), numel(ph.Outport));
dp = fieldnames(get_param([mdl '/EXT_sps'], 'DialogParameters'));
fprintf('SPS params: %s\n', strjoin(dp', ' | '));
for k = 1:numel(dp)
    v = '';
    try, v = get_param([mdl '/EXT_sps'], dp{k}); catch, end
    if ischar(v) && ~isempty(v), fprintf('  %s = %s\n', dp{k}, v); end
end
% what feeds each mux input (block names)?
pc = get_param([mdl '/EXT_Mux'], 'PortConnectivity');
for k = 1:numel(pc)
    src = '<none>';
    try
        h = pc(k).SrcBlockHandle;
        if isscalar(h) && h > 0, src = get_param(h, 'Name'); end
    catch, end
    fprintf('mux in%d <- %s\n', k, src);
end
close_system(mdl, 0);
end
