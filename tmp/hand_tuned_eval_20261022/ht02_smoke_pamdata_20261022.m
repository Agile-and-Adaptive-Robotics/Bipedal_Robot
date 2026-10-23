% Hand-tuned BPA evaluation, step 2: SMOKE TEST of the staged
% PamData_20261022 mechanics on eb475ws4 (2026-10-22).
% Questions this run answers:
%  (a) MATLAB semantics: A(:,k) on 3-D arrays; mtimes with N-D operand
%      (decides whether PamData's computeMomentArm line 140 can handle
%      the 3-D Locations Connor's Bi_Hip builds for rect_fem / vas_int).
%  (b) Does PamData_20261022 construct bifemsh (2-D Location), rect_fem
%      and vas_int (3-D Locations) on a small Bi_Hip-style T?
%  (c) Do Connor's rebuilt fit functions (connor_Bi_Hip.m data, verbatim)
%      match the fcn1..fcn8 stored in Ben's MovePoints.mat?

root = 'D:\Github\Bipedal_Robot';
here = fullfile(root, 'tmp', 'capability_table_20261022');
oldc = fullfile(root, 'tmp', 'old_code_copy_20261022');
addpath(here);
addpath(oldc);   % READ-ONLY: ForceStrainTable.mat (festo), MovePoints.mat

fprintf('=== (a) MATLAB semantics ===\n');
A = rand(3, 5, 4);
fprintf('size(A(:,2)) for 3x5x4 array: %s\n', mat2str(size(A(:,2))));
try
    z = rand(4,4) * rand(4,1,3);
    fprintf('mtimes 4x4 * 4x1x3: OK, size %s\n', mat2str(size(z)));
catch e
    fprintf('mtimes 4x4 * 4x1x3 FAIL: %s\n', e.message);
end

%% Joint primitives and fits (connor_Bi_Hip.m, VERBATIM data)
p_h = [-0.0707; -0.0661; 0.0835];                     % connor_Bi_Hip.m:27
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];  % :39
knee_x       = [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];  % :40
KneeXFunc = fit(knee_angle_x, knee_x, 'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];  % :42
knee_y       = [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];  % :43
KneeYFunc = fit(knee_angle_y, knee_y, 'cubicspline');
rect_fem_x  = [0.0156367; 0.0179948; 0.0274274; 0.029683; 0.0306; 0.0366; 0.0422; 0.0451; 0.0484; 0.0533; 0.0617; 0.0634; 0.067; 0.0733];        % :108
rect_fem_xD = pi/180*[-120.118; -114.871; -90.068; -83.532; -80; -60; -40; -30; -20; -10; 0; 1.6; 5; 10];                                        % :109
fcn3 = fit(rect_fem_xD, rect_fem_x, 'smoothingspline');                                                                                        % :110
rect_fem_y  = [0.0234; 0.0238; 0.0251; 0.0253; 0.025284; 0.0249; 0.0243; 0.0239; 0.0234; 0.0228; 0.0210; 0.0206; 0.0192; 0.0160];               % :111
rect_fem_yD = pi/180*[-120; -114.6; -90; -83.5; -80.01; -60; -40; -30; -20; -10; 0; 1.6; 5; 10 ];                                              % :112
fcn4 = fit(rect_fem_yD, rect_fem_y, 'smoothingspline');                                                                                        % :113
fcn7 = fit(rect_fem_xD, rect_fem_x-(0.0156367-0.008285), 'smoothingspline');                                                                   % :152
fcn8 = fit(rect_fem_yD, rect_fem_y+(0.0256239-0.02346), 'smoothingspline');                                                                    % :153

% Rodrigues rotation (JointData.computeRotMat formula, verbatim, unit axes)
Rax = @(ax, th) [cos(th)+ax(1)^2*(1-cos(th)), ax(1)*ax(2)*(1-cos(th))-ax(3)*sin(th), ax(1)*ax(3)*(1-cos(th))+ax(2)*sin(th);
                 ax(1)*ax(2)*(1-cos(th))+ax(3)*sin(th), cos(th)+ax(2)^2*(1-cos(th)), ax(2)*ax(3)*(1-cos(th))-ax(1)*sin(th);
                 ax(1)*ax(3)*(1-cos(th))-ax(2)*sin(th), ax(2)*ax(3)*(1-cos(th))+ax(1)*sin(th), cos(th)+ax(3)^2*(1-cos(th))];

%% Small Bi_Hip-style T: 2 hip x 3 knee poses (connor_Bi_Hip.m:51-59 layout)
div1 = 2; div2 = 3;
thH = linspace(-15, 110, div1) * pi/180;
thK = linspace(-120, 10, div2) * pi/180;
T = zeros(4, 4, 2, div1*div2);
ii = 1;
for i = 1:div1
    for j = 1:div2
        T(:, :, 1, ii) = [Rax([0;0;1], thH(i)), p_h; 0 0 0 1];
        kneeHome = [KneeXFunc(thK(j)); KneeYFunc(thK(j)); 0];
        T(:, :, 2, ii) = [Rax([0;0;1], thK(j)), kneeHome; 0 0 0 1];
        ii = ii + 1;
    end
end
fprintf('\n=== (b) PamData_20261022 constructions ===\n');

% bifemsh (connor_Bi_Hip.m:141-147, VERBATIM)
LocationBif = [-0.007, -0.03, -0.023;
              -0.096, -0.036, -0.056;
               0.032,  0.029,  0.029];
try
    m = PamData_20261022('Bicep Femoris, Short Head', LocationBif, [0, 2], 804, T, [0 0 0; 0, 0, 30]);
    fprintf('bifemsh (2-D): OK  Dia=%g Rest=%.6g Long=%.6g  L=%s  MA(2,:,3)=%s\n', ...
        m.Diameter, m.RestingL, m.LongestL, mat2str(m.MuscleLength, 5), ...
        mat2str(squeeze(m.MomentArm(2, :, 3))', 5));
catch e
    fprintf('bifemsh (2-D) FAIL: %s\n', e.message);
end

% rect_fem (connor_Bi_Hip.m:115-128, VERBATIM, 3-D Location)
LocationRec = zeros(3, 5, size(T, 4));
ii = 1;
for i = 1:div1
    for j = 1:div2
        LocationRec(:, :, ii) = [-0.072, 0.072, 0.065, 0.058, fcn3(thK(j));
                                  0.125, -0.117, -0.156, -0.329, fcn4(thK(j));
                                  0.079, 0.015, 0.005, -0.003, 0.002];
        ii = ii + 1;
    end
end
try
    m5 = PamData_20261022('Rectus Femoris', LocationRec, [2, 5], 1169, T, [10 20 30; 0, 0, 30]);
    fprintf('rect_fem (3-D): OK  Dia=%g Rest=%.6g Long=%.6g  L=%s\n', ...
        m5.Diameter, m5.RestingL, m5.LongestL, mat2str(m5.MuscleLength, 5));
catch e
    fprintf('rect_fem (3-D) FAIL: %s\n', e.message);
end

% vas_int (connor_Bi_Hip.m:154-167, VERBATIM, 3-D Location)
LocationVas = zeros(3, 4, size(T, 4));
ii = 1;
for i = 1:div1
    for j = 1:div2
        LocationVas(:, :, ii) = [0.0, 0.057, 0.054, fcn7(thK(j));
                                -0.043, -0.216, -0.332, fcn8(thK(j));
                                 0.042, 0.029, 0.013, 0.002];
        ii = ii + 1;
    end
end
try
    m8 = PamData_20261022('Vastus Intermedius', LocationVas, [0, 4], 4530, T, [0 0 0; 0, 0, 30]);
    fprintf('vas_int (3-D): OK  Dia=%g Rest=%.6g Long=%.6g  L=%s\n', ...
        m8.Diameter, m8.RestingL, m8.LongestL, mat2str(m8.MuscleLength, 5));
catch e
    fprintf('vas_int (3-D) FAIL: %s\n', e.message);
end

%% (c) Connor's rebuilt fits vs Ben's MovePoints.mat fcn1..fcn8
fprintf('\n=== (c) rebuilt fits vs MovePoints.mat fcns ===\n');
MP = load(fullfile(oldc, 'MovePoints.mat'));
qs = [-2.0 -1.2 -0.5 0 0.5 1.5];
pairs = {'KneeXFunc','fcn1'; 'KneeYFunc','fcn2'; 'fcn3','fcn3'; 'fcn4','fcn4'; 'fcn7','fcn7'; 'fcn8','fcn8'};
for p = 1:size(pairs, 1)
    rebuilt = eval(pairs{p, 1});         %#ok<EVLCOM>  (this smoke script only)
    stored  = MP.(pairs{p, 2});
    dr = max(abs(arrayfun(@(q) rebuilt(q), qs) - arrayfun(@(q) stored(q), qs)));
    fprintf('%-10s vs mat %-6s : max|diff| over [%s] = %.3g\n', pairs{p,1}, pairs{p,2}, ...
        strjoin(arrayfun(@(x) sprintf('%.2f', x), qs, 'UniformOutput', false), ' '), dr);
end
fprintf('SMOKE DONE\n');
