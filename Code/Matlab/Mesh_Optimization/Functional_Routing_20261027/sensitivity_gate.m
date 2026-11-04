function sensitivity_gate(round)
%SENSITIVITY_GATE Shim: addpath the campaign machinery, then run the gate.
%
% The gate's contract command runs with cd('.../Functional_Routing_20261027'),
% so MATLAB resolves sensitivity_gate from THIS folder; the full
% implementation lives in campaign_20261103/sensitivity_gate.m (everything
% of the machinery itself is under campaign_20261103/, per the campaign
% folder rule - this shim is the minimal entry point the exact runnable
% command requires).
addpath(fullfile(fileparts(mfilename('fullpath')), ...
    'campaign_20261103'));
fr_sensitivity_gate(round);
end
