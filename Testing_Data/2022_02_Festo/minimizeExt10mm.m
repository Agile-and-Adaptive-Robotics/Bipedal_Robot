%% minimizeExt10mm function
function [f_all, bpa_all] = minimizeExt10mm(Xi0, Xi1, Xi2, varargin)
% minimizeExt10mm: Wrapper for evaluating multiple BPAs, optionally with validation split
% Inputs:
%   Xi0, Xi1, Xi2 - optimization parameters
%   varargin - optional list of BPA indices to use as validation set
% Outputs:
%   f_all - Nx3 matrix of [RMSE, FVU, MaxResidual]
%   bpa_all - structure array of BPA structs with predicted fields

%% Load BPA dataset
load ExtPinBPASet.mat ke  % loads the ke(1:4) structs
nBPA = numel(ke);

% Determine which indices to validate vs optimize
if ~isempty(varargin)
    idx_val = varargin{1};
else
    idx_val = 1:nBPA;  % default: treat all as validation
end

% Call core minimization
% disp('--- Re-evaluating each BPA with new parameters ---');
% disp(['Xi0 = ', num2str(Xi0), ', Xi1 = ', num2str(Xi1), ', Xi2 = ', num2str(Xi2)]);
[f_all, bpa_all] = minimizeExt(Xi0, Xi1, Xi2, idx_val);
% assert(any(any(diff(f_all(~isnan(f_all(:,1)), :)))), 'All BPAs evaluated identically — likely loop bug.');
end
