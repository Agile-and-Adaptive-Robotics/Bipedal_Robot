function [tau, T] = fr_loadSignedTargets(groupKey, grid, targetsDir)
%FR_LOADSIGNEDTARGETS Load a group's SIGNED per-pose human torque targets.
%
% Directive item 2 deliverable (task step 5). Loads
% targets_v2/<groupKey>_signed_targets.json (built by
% build_signed_targets.py from the STOCK human model:
% tau_target(cell) = SUM_m F_m(q) * r_m,dof(q), r = OpenSim
% computeMomentArm == -dL/dq, F_m = Fmax * active force-length at the
% rigid-tendon CE length) and returns the per-cell tau_target vector
% AFTER validating it against the caller's grid so a consumer can never
% silently score against a misaligned curve.
%
% INPUTS
%   groupKey   group name of record (e.g. 'vas_int'; the contract's
%              "group" key, which is also the fr_cases27 case name).
%   grid       EITHER a numeric vector of the pass-1 swept-coordinate
%              values (radians, fr_cases cell order) - typically
%              [cs.passes{1}.th1] - OR an fr_cases27 case struct with
%              .passes{1}.th1/.qStructs and .passes{1}.sweptCoord.
%   targetsDir optional folder holding the JSONs (default: the
%              targets_v2 folder next to this file).
%
% OUTPUTS
%   tau  nCells x 1 signed target torques (N.m). Positive tau drives the
%        primary coordinate in its positive direction per the model
%        (knee_angle_r: positive = extension, i.e. flexion NEGATIVE).
%   T    the full parsed group record (cells with per-constituent
%        breakdown, convention and force-length documentation strings,
%        peak fields). cs.target.tau := tau is the drop-in for fr_objective.
%
% ERRORS LOUDLY on: missing JSON, cell-count mismatch, swept-coordinate
% mismatch, or any per-cell coordinate difference > 1e-9 rad. The grid
% compared against is the fr_cases27 REQUESTED value (the qStructs store
% the request; the OpenSim-side realized/clamped values are recorded in
% the JSON per cell - only the knee -120 deg request is clamped, onto the
% model's stored range endpoint).

if nargin < 3
    [here, ~] = fileparts(mfilename('fullpath'));
    targetsDir = fullfile(here, 'targets_v2');
end

pth = fullfile(targetsDir, sprintf('%s_signed_targets.json', groupKey));
if exist(pth, 'file') ~= 2
    error('fr_loadSignedTargets:Missing', ...
        ['no signed-target JSON for group %s at %s - run ' ...
        'build_signed_targets.py first'], groupKey, pth);
end
T = jsondecode(fileread(pth));

% ---- resolve the caller's grid into (coord, values) ----------------------
if isstruct(grid) && isfield(grid, 'passes')
    pass1 = grid.passes{1};
    coord = pass1.sweptCoord;
    vals = reshape(pass1.th1, 1, []);
else
    vals = double(reshape(grid, 1, []));
    coord = T.primaryDof;
end
vals = vals(:);

% ---- validations (all loud) ----------------------------------------------
nJ = numel(T.cells);
if isfield(T, 'primaryDof') && ~strcmp(coord, T.primaryDof)
    error('fr_loadSignedTargets:Dof', ...
        ['group %s: caller grid sweeps coordinate %s but the target ' ...
        'record is about %s'], groupKey, coord, T.primaryDof);
end
if numel(vals) ~= nJ
    error('fr_loadSignedTargets:Count', ...
        ['group %s: caller grid has %d cells but the target record ' ...
        'has %d - the consumers index targets by cell, so this is a ' ...
        'hard mismatch'], groupKey, numel(vals), nJ);
end
req = zeros(nJ, 1);
for k = 1:nJ
    req(k) = T.cells(k).requested_rad;
end
dTol = 1e-9;   % rad
dMax = max(abs(vals - req));
if dMax > dTol
    [~, wi] = max(abs(vals - req));
    error('fr_loadSignedTargets:GridMismatch', ...
        ['group %s: caller grid does not match the fr_cases27 pass-1 ' ...
        'cells (max |dq| = %.3e rad at cell %d; tolerance %.1e). ' ...
        'Expected [rad]: %s; got: %s'], groupKey, dMax, wi, dTol, ...
        mat2str(req', 12), mat2str(vals', 12));
end

% ---- return ---------------------------------------------------------------
tau = T.tau_vec_Nm(:);
if numel(tau) ~= nJ
    error('fr_loadSignedTargets:Corrupt', ...
        'group %s: tau vector length %d != cell count %d', ...
        groupKey, numel(tau), nJ);
end
end
