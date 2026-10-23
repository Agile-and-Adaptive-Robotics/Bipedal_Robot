% Probe Survey2 Sizes Dofs 20261023
% Session probe (scratch) for the survey2 driver: (1) per-muscle OLD-CODE
% BPA diameter from PAM_Prop.mat MIF via Size.m's own rule (Fc = 0.9
% thresholds on the Festo theoretical maxima); (2) per-actuator dof-pair
% mapping check: surfaces CSV (primary, secondary) vs the masters
% builder's kin.dofs slots -- which slot models the primary, whether the
% secondary is modeled; (3) domain exclusions per actuator from the CSV.

scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot'), break; end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(scriptDir);

%% 1. Old-code diameters (Size.m rule, PAM_Prop data of record)
P = load(fullfile(root, 'tmp', 'old_code_copy_20261022', 'PAM_Prop.mat'));
F1 = 630; F2 = 1500; F3 = 6000; Fc = 0.9;   % Size.m constants verbatim
fprintf('\n=== OLD-CODE DIAMETERS (Size.m rule on PAM_Prop MIF) ===\n');
names = biPulleySpecsFromMasters_20261023();   % map order, ids 1..27
diaOf = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:numel(names)
    nm = names(k).name;
    mif = names(k).fmax;   % = PAM_Prop.(mifField), the builder's own load
    if mif < Fc*F1
        dia = 10;
    elseif mif >= Fc*F1 && mif < Fc*F2
        dia = 20;
    elseif mif >= Fc*F2 && mif <= Fc*F3
        dia = 40;
    else
        dia = NaN;   % 'Use additional muscles'
    end
    diaOf(nm) = dia;
    fprintf('%2d %-12s mif %8.1f -> %s\n', names(k).actuatorId, nm, mif, ...
        iif(isnan(dia), 'ERROR(Use additional muscles)', sprintf('%d mm', dia)));
end

%% 2. Surfaces CSV dof pairs + domain mask vs builder slots
Sraw = readtable(fullfile(scriptDir, 'human_surfaces_5x5_grading_20261023.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char', ...
    'CommentStyle', '#');
fprintf('\n=== DOF MAPPING (CSV pair vs builder slots) ===\n');
fprintf('%-3s %-12s %-18s %-18s %-22s %-18s %s\n', 'id', 'spec', 'csv primary', ...
    'csv secondary', 'builder dofs', 'secondary modeled?', 'excluded cells');
for k = 1:numel(names)
    sp = names(k);
    ix = Sraw.actuator_id == sp.actuatorId;
    pri = strtrim(Sraw.primary_dof{find(ix, 1)});
    sec = strtrim(Sraw.secondary_dof{find(ix, 1)});
    nEx = sum(Sraw.in_graded_domain(ix) ~= 1);
    priSlot = find(strcmp(sp.kin.dofs, pri), 1);
    secSlot = find(strcmp(sp.kin.dofs, sec), 1);
    if isempty(priSlot)
        pmap = 'PRIMARY NOT MODELED';
    elseif priSlot == sp.kin.primarySlot
        pmap = sprintf('slot %d (== primarySlot)', priSlot);
    else
        pmap = sprintf('slot %d (=/= primarySlot %d)', priSlot, sp.kin.primarySlot);
    end
    if isempty(secSlot)
        smap = 'NO (robot constant along it)';
    else
        smap = sprintf('YES slot %d', secSlot);
    end
    fprintf('%-3d %-12s %-18s %-18s [%-18s|%s] %-18s %d  %s\n', sp.actuatorId, ...
        sp.name, pri, sec, sp.kin.dofs{1}, sp.kin.dofs{2}, smap, nEx, pmap);
end

function s = iif(c, a, b)
if c, s = a; else, s = b; end
end
