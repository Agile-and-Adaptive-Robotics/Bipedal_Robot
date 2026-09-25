%% Run this before the reverse-pulley optimizer. NO optimization here.
% Opt_sanity_pulley.m -- sanity gate for MonoPam_pulley (Ben, 2026-09-24).
% Builds a 92-orientation knee geometry (Knee_Flexor_Data_20mm.m's T_Pam
% block) with a straight-tendon route whose crossing segment IS the tendon
% line (3 rows: BPA origin and tackle exit on the proximal body, insertion
% on the distal body; CrossPoint = 3, PulleyExitIndex = 2), then asserts:
%   (a) nPulleyBPA = 1, G = 1 reproduces MonoPamDataExplicit_balance to
%       1e-8 RELATIVE on F_p, mA_p, Torque_p (the regression identity);
%   (b) at G = 2, nPulleyBPA = 1 the solved equilibrium satisfies
%       F_t = F_BPA/G and delta_t = 2*s - DeltaL at each orientation
%       (residuals recomputed from the solved state, tol 1e-6; s = the
%       exposed tackle input travel PulleyTravel);
%   (c) at nPulleyBPA = 2, G = 1 the tendon tension equals the SUM of the
%       two BPA forces (total pull doubles, travel gain does not, tol
%       1e-6);
%   (d) the reaction force magnitude equals |nBPA*F_BPA*u_bpa + F_t*u_t|
%       (tol 1e-9);
%   (e) in bowden mode the insertion-side unit direction u_t is constant
%       in the tibia frame across orientations (tol 1e-12) while the
%       closure still holds (tol 1e-6);
%   (f) an infeasible case sets PulleyInfeasible and NaN torque without
%       erroring. By the closure the required contraction is
%       s0ref + (DeltaL + delta_t)/G + cb*F, so infeasibility arrives when
%       the SPAN travel exceeds G*(KMAX*Rest - s0ref) -- it can NOT be
%       reached by making G large (larger G divides the span demand);
%       the case is therefore constructed by sizing rest so that boundary
%       sits inside the sweep.
% No figures, no parpool; matlab -batch runnable in under ~3 minutes.

clear
clear functions
clc
rehash

%% Add paths to the muscle and pam calculators
% Repo root and path setup (same block as Results/Knee_Flexor_Data_20mm.m).
scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root)
        error('Could not locate the Bipedal_Robot repo root from %s', scriptDir)
    end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
% Mesh_Optimization must win any shadowing contest against data subfolders.
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
% Append (do not prepend) so Code\Matlab keeps winning name collisions.
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

%% Joint rotation transformation matrices (92-orientation T_Pam block)
% Same fits and construction as Knee_Flexor_Data_20mm.m, with the position
% count set to 92 (the orientation count the base class's Lok hard-codes;
% MonoPam_pulley derives it from size(Location,3) instead).
positions = 92;
fprintf('OPT_SANITY_PULLEY: building a %d-orientation knee geometry.\n', positions)

R_Pam = zeros(3, 3, positions);
T_Pam = zeros(4, 4, positions);
T_t1_ICR = zeros(4, 4, positions);
T_ICR_t1 = zeros(4, 4, positions);

c = pi/180; %Convert from degrees to radians

%Robot
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam =     ([23.30	22.22	21.55	21.09	19.91	18.70	17.48	13.82	10.44	7.60	5.52	4.35	4.16	5.01	7.04	10.47]')/1000;
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     ([-416.65	-417.03	-417.19	-417.28	-417.41	-417.41	-417.30	-416.28	-414.36	-411.72	-408.62	-405.32	-402.08	-399.16	-396.85	-395.66]')/1000;
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');

%Theta1 to ICR
t1_ICR_x = ([29.66	28.54	27.86	27.40	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]')/1000;
fcn13 = fit(knee_angle,t1_ICR_x,'cubicspline');
t1_ICR_y = ([25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	21.32	20.99	21.2	22.33]')/1000;
fcn14 = fit(knee_angle,t1_ICR_y,'cubicspline');

% Reduced sweep: deep flexion is not needed for the transmission checks,
% and staying inside the strain window keeps every frame a clean interior
% fzero root (no escape paths on either class).
kneeMin = -45*c;
kneeMax = 10*c;
phi = linspace(kneeMin, kneeMax, positions);
%We want one of our positions to be home position, so let's make the
%smallest value of phi equal to 0
[val, pos] = min(abs(phi));
phi(pos) = 0;

for i = 1:positions
    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    R_Pam(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;   %Rotation matrix for robot
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];

    T_Pam(:, :, i) = RpToTrans(R_Pam(:, :, i), hipToKnee_Pam');     %Transformation matrix for robot

    t1toICR = [fcn13(phi(i)), fcn14(phi(i)), 0]; %distance from theta1 to ICR
    T_t1_ICR(:, :, i) = RpToTrans(eye(3), t1toICR');    %transform from the ICR frame to theta1
    T_ICR_t1(:, :, i) = RpToTrans(eye(3), -t1toICR');    %transform from t1 frame to ICR
end

phiD = phi*180/pi;

%% Straight-tendon route: crossing segment = the tendon line
% Three rows: row 1 = BPA origin (proximal body, femur frame), row 2 =
% tackle exit (proximal body, femur frame), row 3 = insertion (tibia,
% theta1 frame), CrossPoint = 3, PulleyExitIndex = 2. Rows 1-2 are held
% CONSTANT across orientations so the rigid route change is exactly the
% tendon-segment (exit-to-insertion) span change -- the condition that
% makes the G = 1 regression identity algebraic. The route is sized so
% the modeled strain stays inside (0, KMAX) over the sweep.
p1 = [-0.050, 0.390, 0.050];        %BPA origin, femur frame
pExit = [-0.040, 0.360, 0.048];     %tackle exit, femur frame
p2 = [0.0574, 0.0355, 0.005];       %insertion, theta1 frame

Location = zeros(3, 3, positions);
routeLength = zeros(positions, 1);
for i = 1:positions
    p2ICR = RowVecTrans(T_ICR_t1(:,:,i), p2);
    Location(:,:,i) = [p1; pExit; p2ICR];
    % Total length in ONE frame (the class's segment convention: the
    % distal row transformed into the proximal frame by T_Pam).
    routeLength(i) = norm(p1 - pExit) + ...
        norm(pExit - RowVecTrans(T_Pam(:,:,i), p2ICR));
end
CrossPoint = 3;
exitIndex = 2;

%% BPA parameters
Name = 'Bicep Femoris (Short Head)';
Dia = 20;
tendon0 = 0.015;
fitting = 0.021;
pres = 620;
wraps = 6;
KMAX = 0.255;

% Identified stiffness triple (same source as buildKneeFlexorContext20mm:
% 20260908 2brkt 2trans noT3 front, pick 77), with a deterministic
% fallback so the gate still runs if the mat moves.
Xi0 = 0.00394; Xi1 = 3.998e4; Xi2 = 1.473e4;
try
    S = load('minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat', ...
        'filtered_results', 'xCols');
    g = S.filtered_results(77, S.xCols);
    Xi0 = g(1); Xi1 = g(2); Xi2 = g(3);
    clear S g
    fprintf('Stiffness from the 20260908 noT3 front, pick 77.\n')
catch
    fprintf('Stiffness mat not found; using fallback Xi (identity is Xi-independent).\n')
end

% Size the rest length so the reference pose (longest span) is exactly the
% BPA's rest configuration (s0ref = 0, the contract's convention) and the
% sweep stays inside (0, KMAX) on both classes.
rest0 = max(routeLength) - Xi0 - tendon0 - 2*fitting;
kmax0 = (1 - KMAX)*rest0;
fprintf(['Route length range = %.6f to %.6f m; rest = %.6f m; ' ...
    'kmax = %.6f m\n'], min(routeLength), max(routeLength), rest0, kmax0)

%% Pulley configs under test
cfgG1 = struct('nPulleyBPA', 1, 'tackleLineParts', 1, ...
    'pulleyExitIndex', exitIndex);                    %pulley disabled
cfgG2 = struct('nPulleyBPA', 1, 'tackleLineParts', 2, ...
    'pulleyExitIndex', exitIndex);                    %2:1 tackle
cfg2BPA = struct('nPulleyBPA', 2, 'tackleLineParts', 1, ...
    'pulleyExitIndex', exitIndex);                    %two parallel BPAs
cfgBowden = struct('nPulleyBPA', 1, 'tackleLineParts', 2, ...
    'routingMode', 'bowden', 'pulleyExitIndex', exitIndex);

%% Build the straight-tendon baseline and the G = 1 pulley twin
baseBPA = MonoPamDataExplicit_balance(Name, Location, CrossPoint, Dia, ...
    T_Pam, rest0, kmax0, tendon0, fitting, pres, Xi0, Xi1, Xi2, wraps);
pulley1 = MonoPam_pulley(Name, Location, CrossPoint, Dia, ...
    T_Pam, rest0, kmax0, tendon0, fitting, pres, Xi0, Xi1, Xi2, wraps, cfgG1);

%% Assertion (a): G = 1 regression identity, rel tol 1e-8
relTol = 1e-8;
fields = {'F_p', 'mA_p', 'Torque_p'};
for k = 1:numel(fields)
    a = pulley1.(fields{k});
    b = baseBPA.(fields{k});
    if any(~isfinite(a(:))) || any(~isfinite(b(:)))
        error('Opt_sanity_pulley:NonFinite', ...
            ['Assertion (a): %s contains nonfinite values (pulley max ' ...
             'abs %g, base max abs %g) -- the sweep left the clean strain ' ...
             'window.'], fields{k}, max(abs(a(:))), max(abs(b(:))))
    end
    relErr = norm(a(:) - b(:))/max(norm(b(:)), eps);
    fprintf('Assertion (a): %8s relative difference = %.3e\n', fields{k}, relErr)
    if ~(relErr < relTol)
        error('Opt_sanity_pulley:RegressionIdentity', ...
            ['Assertion (a) FAILED: G = 1 %s relative difference %.3e ' ...
             'exceeds tol %.1e vs MonoPamDataExplicit_balance.'], ...
            fields{k}, relErr, relTol)
    end
end
fprintf('Assertion (a) PASSED: nBPA = 1, G = 1 reproduces the balance class (rel tol %.1e).\n', relTol)

%% Assertion (b): G = 2 equilibrium closure, tol 1e-6
Gtest = 2.0;
pulley2 = MonoPam_pulley(Name, Location, CrossPoint, Dia, ...
    T_Pam, rest0, kmax0, tendon0, fitting, pres, Xi0, Xi1, Xi2, wraps, cfgG2);

if any(pulley2.PulleyInfeasible)
    error('Opt_sanity_pulley:UnexpectedInfeasible', ...
        'Assertion (b): the G = 2 sweep must be feasible everywhere to test the closure.')
end
if any(pulley2.PulleySlack)
    error('Opt_sanity_pulley:UnexpectedSlack', ...
        'Assertion (b): the G = 2 sweep must keep the tendon taut everywhere.')
end

% Recompute both sides independently from the solved state.
Fbpa2 = festo4(Dia, pulley2.sContraction(:)/rest0/KMAX, pres) .* pulley2.Fmax;
tolClosure = 1e-6;
% F_t = F_BPA/G with both sides recomputed independently
resFt = pulley2.kSpr .* pulley2.gama(:) - Fbpa2/Gtest;
% delta_t = G*s - DeltaL with s = PulleyTravel (the tackle input travel;
% delta_t is the absolute stretch kSpr*gama, referenced to the blocked
% rope state at theta0).
resClosure = (Gtest .* pulley2.PulleyTravel(:) - pulley2.deltaL(:)) ...
    - pulley2.gama(:);
fprintf(['Assertion (b): max |kSpr*gama - F_BPA/G| = %.3e N, ' ...
    'max |(G*PulleyTravel - DeltaL) - gama| = %.3e m\n'], ...
    max(abs(resFt)), max(abs(resClosure)))
if ~(max(abs(resFt)) < tolClosure) || ~(max(abs(resClosure)) < tolClosure)
    error('Opt_sanity_pulley:EquilibriumClosure', ...
        ['Assertion (b) FAILED: equilibrium residuals %.3e N / %.3e m ' ...
         'exceed tol %.1e at G = %.1f.'], ...
        max(abs(resFt)), max(abs(resClosure)), tolClosure, Gtest)
end
fprintf('Assertion (b) PASSED: F_t = F_BPA/G and delta_t = 2*s - DeltaL (tol %.1e).\n', ...
    tolClosure)

%% Assertion (c): nPulleyBPA = 2, G = 1 -- tension sums, travel gain not
pulley2n = MonoPam_pulley(Name, Location, CrossPoint, Dia, ...
    T_Pam, rest0, kmax0, tendon0, fitting, pres, Xi0, Xi1, Xi2, wraps, cfg2BPA);

if any(pulley2n.PulleyInfeasible) || any(pulley2n.PulleySlack)
    error('Opt_sanity_pulley:Unexpected2BPA', ...
        ['Assertion (c): the nPulleyBPA = 2, G = 1 sweep must be feasible ' ...
         'and taut everywhere to test the tension sum.'])
end

Fsingle2n = festo4(Dia, pulley2n.sContraction(:)/rest0/KMAX, pres) ...
    .* pulley2n.Fmax;
tolPull = 1e-6;
% Total pull doubles: F_t = 2*F_single (recomputed independently).
resPull = pulley2n.kSpr .* pulley2n.gama(:) - 2*Fsingle2n;
% Travel gain does not: the same closure holds with G = 1.
resTravel = (pulley2n.PulleyTravel(:) - pulley2n.deltaL(:)) - pulley2n.gama(:);
fprintf(['Assertion (c): max |kSpr*gama - 2*F_single| = %.3e N, ' ...
    'max |(PulleyTravel - DeltaL) - gama| = %.3e m\n'], ...
    max(abs(resPull)), max(abs(resTravel)))
if ~(max(abs(resPull)) < tolPull) || ~(max(abs(resTravel)) < tolPull)
    error('Opt_sanity_pulley:TensionSum', ...
        ['Assertion (c) FAILED: tension-sum residual %.3e N or travel ' ...
         'residual %.3e m exceeds tol %.1e at nPulleyBPA = 2, G = 1.'], ...
        max(abs(resPull)), max(abs(resTravel)), tolPull)
end
fprintf(['Assertion (c) PASSED: F_tendon = F_BPA{1} + F_BPA{2} and the ' ...
    'travel gain is unchanged (tol %.1e).\n'], tolPull)

%% Assertion (d): reaction force magnitude identity, tol 1e-9
% Checked on BOTH the G = 2 (nPulleyBPA = 1) and the nPulleyBPA = 2
% (G = 1) solves: ||ReactionF|| must equal
% ||nBPA*F_BPA*u_bpa + F_t*u_t|| pointwise.
tolReaction = 1e-9;
resReaction = zeros(positions, 2);
for i = 1:positions
    rG2 = norm(pulley2.FsingleBPA(i)*pulley2.u_bpa(i,:) ...
        + pulley2.Ftendon(i)*pulley2.u_t(i,:));
    r2n = norm(2*Fsingle2n(i)*pulley2n.u_bpa(i,:) ...
        + pulley2n.Ftendon(i)*pulley2n.u_t(i,:));
    resReaction(i, 1) = abs(pulley2.ReactionFmag(i) - rG2);
    resReaction(i, 2) = abs(pulley2n.ReactionFmag(i) - r2n);
end
fprintf(['Assertion (d): max | ||ReactionF|| - ||F_BPA*u_bpa + ' ...
    'F_t*u_t|| | = %.3e (G = 2) / %.3e (nBPA = 2) N\n'], ...
    max(resReaction(:,1)), max(resReaction(:,2)))
if ~(max(resReaction(:)) < tolReaction)
    error('Opt_sanity_pulley:ReactionIdentity', ...
        ['Assertion (d) FAILED: reaction identity residual %.3e N ' ...
         'exceeds tol %.1e.'], max(resReaction(:)), tolReaction)
end
fprintf('Assertion (d) PASSED: reaction magnitude identity (tol %.1e).\n', tolReaction)

%% Assertion (e): bowden mode -- u_t constant in the tibia frame
pulleyB = MonoPam_pulley(Name, Location, CrossPoint, Dia, ...
    T_Pam, rest0, kmax0, tendon0, fitting, pres, Xi0, Xi1, Xi2, wraps, cfgBowden);

u_tDeviation = vecnorm(pulleyB.u_t - pulleyB.u_t(1,:), 2, 2);
resClosureB = (Gtest .* pulleyB.PulleyTravel(:) - pulleyB.deltaL(:)) ...
    - pulleyB.gama(:);
fprintf(['Assertion (e): max u_t deviation across orientations = %.3e ' ...
    '(tol 1e-12), max closure residual = %.3e m (tol %.1e)\n'], ...
    max(u_tDeviation), max(abs(resClosureB)), tolClosure)
if ~(max(u_tDeviation) < 1e-12)
    error('Opt_sanity_pulley:BowdenDirection', ...
        ['Assertion (e) FAILED: bowden u_t varies by %.3e across ' ...
         'orientations; it must be constant in the tibia frame (tol 1e-12).'], ...
        max(u_tDeviation))
end
if ~(max(abs(resClosureB)) < tolClosure)
    error('Opt_sanity_pulley:BowdenClosure', ...
        ['Assertion (e) FAILED: the bowden closure residual %.3e m ' ...
         'exceeds tol %.1e -- routing mode must not touch the balance.'], ...
        max(abs(resClosureB)), tolClosure)
end
fprintf('Assertion (e) PASSED: bowden u_t constant in the tibia frame, closure intact.\n')

%% Assertion (f): infeasible case flags and NaNs without erroring
% The infeasible boundary is DeltaL > G*(KMAX*Rest - s0ref) (see header):
% size rest so that boundary sits INSIDE the sweep -- half the span travel
% -- so the deep-flexion frames (largest DeltaL) go infeasible while theta0
% stays feasible, proving the flag discriminates instead of blanket-failing.
mLref = routeLength(pulley2.RefIndex);
restInf = (mLref - Xi0 - tendon0 - 2*fitting - 0.5*max(pulley2.deltaL)/Gtest) ...
    /(1 - KMAX);
kmaxInf = (1 - KMAX)*restInf;
pulleyInf = MonoPam_pulley(Name, Location, CrossPoint, Dia, ...
    T_Pam, restInf, kmaxInf, tendon0, fitting, pres, Xi0, Xi1, Xi2, wraps, cfgG2);

nInfeasible = nnz(pulleyInf.PulleyInfeasible);
s0refInf = restInf - (mLref - Xi0 - tendon0 - 2*fitting);
fprintf(['Assertion (f): rest = %.6f m puts %d of %d frames infeasible ' ...
    '(condition: DeltaL > G*(KMAX*Rest - s0ref) = %.6f m; ' ...
    'max DeltaL = %.6f m).\n'], ...
    restInf, nInfeasible, positions, ...
    Gtest*(KMAX*restInf - s0refInf), max(pulleyInf.deltaL))
if nInfeasible == 0
    error('Opt_sanity_pulley:NoInfeasibleCase', ...
        ['Assertion (f) FAILED: the infeasible geometry produced no ' ...
         'PulleyInfeasible frames (rest = %.6f m, G = %.1f).'], ...
        restInf, Gtest)
end
if nInfeasible == positions
    error('Opt_sanity_pulley:BlanketInfeasible', ...
        ['Assertion (f): the infeasible geometry flagged ALL %d frames -- ' ...
         'size rest so theta0 stays feasible and the flag discriminates.'], ...
        positions)
end
if ~all(isnan(pulleyInf.Torque_ins(pulleyInf.PulleyInfeasible, :)), 'all')
    error('Opt_sanity_pulley:InfeasibleTorqueNotNaN', ...
        'Assertion (f) FAILED: Torque_ins is not NaN on every PulleyInfeasible frame.')
end
if ~all(isnan(pulleyInf.Torque_p(pulleyInf.PulleyInfeasible, :)), 'all')
    error('Opt_sanity_pulley:InfeasibleTorquePNotNaN', ...
        'Assertion (f) FAILED: Torque_p is not NaN on every PulleyInfeasible frame.')
end
feas = ~pulleyInf.PulleyInfeasible;
if ~all(isfinite(pulleyInf.F_ins(feas, :)), 'all') || ...
        ~all(isfinite(pulleyInf.Torque_ins(feas, :)), 'all')
    error('Opt_sanity_pulley:FeasibleFramesNotFinite', ...
        ['Assertion (f) FAILED: frames NOT flagged PulleyInfeasible must ' ...
         'produce finite F_ins and Torque_ins (the flag must discriminate).'])
end
fprintf('Assertion (f) PASSED: infeasible frames flagged, torque NaN, feasible frames finite, no error.\n')

%% Report
fprintf('\n========== OPT_SANITY_PULLEY SUMMARY ==========\n')
fprintf('orientations                = %d\n', positions)
fprintf('route rows / Cross / exit   = %d / %d / %d\n', ...
    size(Location, 1), CrossPoint, exitIndex)
fprintf('rest / kmax                 = %.6f / %.6f m\n', rest0, kmax0)
fprintf('Xi0 / Xi1 / Xi2             = %.6g / %.6g / %.6g\n', Xi0, Xi1, Xi2)
fprintf('reference frame (longest span) = %d (%.3f deg)\n', ...
    pulley1.RefIndex, phiD(pulley1.RefIndex))
fprintf('DeltaL range                = %.6f to %.6f m\n', ...
    min(pulley2.deltaL), max(pulley2.deltaL))
fprintf('G = 2: gama range           = %.3e to %.3e m\n', ...
    min(pulley2.gama), max(pulley2.gama))
fprintf('G = 2: |F_ins| range        = %.3f to %.3f N\n', ...
    min(vecnorm(pulley2.F_ins, 2, 2)), max(vecnorm(pulley2.F_ins, 2, 2)))
fprintf('G = 2: ||ReactionF|| range  = %.3f to %.3f N\n', ...
    min(pulley2.ReactionFmag), max(pulley2.ReactionFmag))
fprintf('nBPA = 2, G = 1: Ftendon range = %.3f to %.3f N\n', ...
    min(pulley2n.Ftendon), max(pulley2n.Ftendon))
fprintf('===============================================\n')

fprintf('OPT_SANITY_PULLEY PASS\n')
