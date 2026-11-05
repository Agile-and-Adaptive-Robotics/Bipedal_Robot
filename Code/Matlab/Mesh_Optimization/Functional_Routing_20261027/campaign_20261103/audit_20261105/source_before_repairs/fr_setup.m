function C = fr_setup()
%FR_SETUP Campaign setup: paths, physics constants, hardware, Xi context.
%
% Stage-1 machinery for the Functional Routing 2026-11-03 campaign
% (MASTER_DIRECTIVE_20261103.md items 3/5/6/7/8/9 as corrected by
% HARDWARE_AND_ZCODE_HANDOFF.md items 8/9/10).
%
% Everything downstream reads its constants from the struct returned here.
% No path here is invented: every folder was read in this session.
%
% fr_setup() also addpaths the repo helpers the machinery reuses
% (maxBPAforce, festo4, bpaR, BiPam_pulley2_20261022) WITHOUT modifying any
% existing file.

repo = fileparts(mfilename('fullpath'));
while ~strcmpi(string(lastFolder(repo)), "Bipedal_Robot")
    parent = fileparts(repo);
    assert(~strcmp(parent, repo), 'fr_setup:RepoRoot', 'Cannot locate repository root');
    repo = parent;
end
camp = fullfile(repo, 'Code', 'Matlab', 'Mesh_Optimization', ...
    'Functional_Routing_20261027', 'campaign_20261103');

% ---- repo helper path (read-only reuse) --------------------------------
addpath(genpath(fullfile(repo, 'Code', 'Matlab', 'Functions')));
addpath(fullfile(repo, 'Code', 'Matlab', 'Mesh_Optimization'), '-end');
addpath(fullfile(repo, 'Code', 'Matlab', 'Robot_Data'), '-end');
addpath(camp, '-end');
for sub = {'kin', 'routes', 'mech', 'obj', 'pack', 'demo'}
    addpath(fullfile(camp, sub{1}), '-end');
end

% ---- paths of record ----------------------------------------------------
C.repo = repo;
C.camp = camp;
C.geom = fullfile(repo, 'Code', 'Matlab', 'Mesh_Optimization', ...
    'Functional_Routing_20261027', 'geometry_20261103');
C.mechDir = fullfile(repo, 'Code', 'Matlab', 'Mesh_Optimization', ...
    'Functional_Routing_20261027', 'mechanics_20261103');
C.targetsDir = fullfile(repo, 'Code', 'Matlab', 'Mesh_Optimization', ...
    'Functional_Routing_20261027', 'targets_20261103');
C.osimStock = fullfile(repo, 'Solid_Models', 'OpenSim', ...
    'Gait2392_Robotbody', 'gait2392_simbody.osim');
C.osimRobot = fullfile(repo, 'Solid_Models', 'OpenSim', ...
    'Gait2392_Robotbody', 'gait2392_robotbody.osim');
C.hardware = fullfile(repo, 'Neuromechanical_Models', ...
    'AH001_Simscape_Repair', 'packaging', 'hardware', ...
    'hardware_candidates.json');
C.boneCloudDir = fullfile(repo, 'Code', 'Matlab', 'Bone_Mesh_Plots', ...
    'Open_Sim_Bone_Geometry');
C.outputDir = fullfile(camp, 'output');
if ~exist(C.outputDir, 'dir'); mkdir(C.outputDir); end

% ---- physics of record (mechanics_20261103 verification) ---------------
C.pressure = 620;                 % kPa (targetPressure of the knee contexts)
C.diameters = [10 20];            % mm, handoff item 1 (no 40 mm)
% Contraction limits, handoff item 2 (16.5% / 24%; Lmin fractions):
C.maxContractionPct = struct('d10', 16.5, 'd20', 24.0);
C.LminFraction = struct('d10', 0.835, 'd20', 0.76);
% Diameter endpoints, handoff item 2 + mechanics_20261103/mech2 findings:
%   20 mm: 20 mm at rest -> 38 mm fully contracted (calibrated endpoints
%   carried by THIS campaign; the repo braid law gives 23.6 mm Dout at
%   k=0 and does NOT reproduce the handoff pair as one relation - see
%   mech2_diameters_20261103.log).
%   10 mm: own braid law (bpaR 10 mm branch, dtheta0 = -0.1 deg): Dout
%   13.6 mm at rest, 17.837 mm at 16.5% - the 20 mm law is NEVER applied.
C.diameterEndpoints = struct( ...
    'd20_rest_mm', 20.0, 'd20_contracted_mm', 38.0, ...
    'd10_rest_mm', 13.6, 'd10_contracted_mm', 17.837, ...
    'd10_law', 'bpaR own 10 mm branch (mech2 log)');
C.diameterLaw = 'linear in contraction between calibrated endpoints';
C.wRap = 0.025;                   % m, bend-arc radius used by Xi3 loss
                                  % (buildGeoExclusion.m geo.wRap)
C.fittingLength = 0.025;          % m, one end fitting (probe_straight fixture)
C.kneeRangeDeg = [-120 0];        % robot knee range, no hyperextension

% ---- adopted Xi context (target_contract.json xiSource; mech6 log) -----
C.Xi = struct('Xi0_flx', 0.00393737, ...     % m  (flexor pick 77)
    'Xi1', 39978.7, ...                      % N/m (shared pair)
    'Xi2', 14733.8, ...                      % N/m (shared pair)
    'Xi3', 0.158289, ...                     % dimensionless wrapping loss
    'Xi0_ext', -0.00637802, ...              % m  (extensor pick 32)
    'source', ['target_contract.json xiSource; mats of record ' 'minimizeFlxPin10_results.mat row 77 / minimizeExt10mmX3_results.mat ' 'row 32; Xi0-Xi3 only, Xi4 waived (handoff item 9)']);

% ---- hardware families (handoff table; hardware_candidates.json) -------
hw = jsondecode(fileread(C.hardware));
C.hardwareCandidates = hw;
C.families = containers.Map();
for i = 1:numel(hw.families)
    f = hw.families(i);
    C.families(f.id) = struct('part', f.sheave_part, ...
        'od_mm', f.od_mm, 'width_mm', f.width_mm, ...
        'capacity_N', f.catalog_capacity_N, ...
        'capacity_basis', f.capacity_basis, ...
        'plate_mm', f.proposed_plate_thickness_mm);
end
% Bowden packaging placeholders (handoff; ferrule length is the library's
% 12 mm PLACEHOLDER, not a vendor number):
C.bowden = struct('housingOD_mm', 5.0, 'ferrule_mm', 12.0, ...
    'ferrule_note', 'PLACEHOLDER from the Blender library, unverified', ...
    'stopAllowance_mm', 20.0, 'bendRadius_mm', 75.0, ...
    'bendRadius_note', 'preview assumption, NOT a vendor rating');

% ---- mount standoff and patch conventions ------------------------------
% Surface-patch search: each movable attachment moves on its owning bone's
% point cloud within a ball of radius patchRadius_mm around the OpenSim
% seed, then is projected to the nearest cloud point plus a mount standoff
% along the local normal, so a Cartesian step cannot carry the mount far
% off the bone (master directive item 3).
C.patchRadiusDefault_mm = 30.0;
C.standoff_mm = [1.0 12.0];       % [min max] mount standoff, millimeters
C.clearanceBone_m = 0.001;        % rigid clearance (buildGeoExclusion.m:20)
C.clearanceEnvelope_m = 0.003;    % hard inter-layer envelope clearance

% ---- solver budgets (stage-1 demonstrator scale) -----------------------
C.solver = struct('MaxFunctionEvaluations', 600, ...
    'MeshTolerance', 1e-3, 'StepTolerance', 1e-3, ...
    'gateMaxFunctionEvaluations', 150, ...
    'refineMaxFunctionEvaluations', 150);
C.refinementRounds = 1;           % coordinated refinement rounds after pass 1
C.refineTolerance = 1e-6;         % allowed degradation of a retained layer's
                                  % objective in a refinement round
C.parallelMaxWorkers = min(10, feature('numcores'));
if strcmpi(getenv('COMPUTERNAME'), 'DESKTOP-5Q16KE9')
    C.parallelMaxWorkers = min(6, C.parallelMaxWorkers);
end
C.rngSeed = 20261103;

% ---- kinematics validation tolerance -----------------------------------
C.kinTolMm = 0.5;   % vs geometry_20261103 bodyTransformsInGround (6 poses)

% ---- logging ------------------------------------------------------------
C.verbose = true;
end

function name = lastFolder(path)
[~, name] = fileparts(path);
end
