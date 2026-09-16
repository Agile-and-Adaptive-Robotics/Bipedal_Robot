function probe_deng_hc()
% Inspect HCNeuron mask + integrator ICs, then try numeric ICs directly.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
load_system(fullfile(here, 'SNS_Deng_RG.slx'));
blk = 'SNS_Deng_RG/pair1_HC_ext';
fprintf('mask params:\n');
for p = {'Cm','Gm','Vrest','GNa','ENa','Sm','Em','Km','Sh','Eh','Kh','tauH','h0'}
    try
        fprintf('  %-6s = %s\n', p{1}, get_param(blk, p{1}));
    catch
        fprintf('  %-6s = <MISSING>\n', p{1});
    end
end
fprintf('Vint IC = "%s"  hint IC = "%s"\n', ...
    get_param([blk '/Vint'], 'InitialCondition'), ...
    get_param([blk '/hint'], 'InitialCondition'));
close_system('SNS_Deng_RG', 0);
end
