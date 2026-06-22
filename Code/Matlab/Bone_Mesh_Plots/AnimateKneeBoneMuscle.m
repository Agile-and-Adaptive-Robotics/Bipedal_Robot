function AnimateKneeBoneMuscle(T, T_ICR_t1, phi, pos, p1, p2, Bifemsh, varargin)
% AnimateKneeBoneMuscle
%
% Human/OpenSim skeleton animation.
%
% T(:,:,i):
%   tibia frame -> femur frame
%
% p1:
%   femur frame
%
% p2:
%   optimized robot/theta1 point
%
% v2:
%   p2 converted into the human tibia/ICR frame at phi(pos)=0
%
% w2:
%   v2 transformed into femur frame at each knee angle
%
% RowVecTrans only accepts a single 4x4 matrix, so every call below uses
% Ti = T(:,:,i), not the whole T array.

    %% Optional inputs
    % inputParser lets you call the function with name/value pairs, e.g.
    %
    %   AnimateKneeBoneMuscle(..., 'PauseTime', 0.25, 'Loop', true)
    %
    % p.Results is the struct holding those parsed option values.
    % I store p.Results as opt so the code is shorter.
    %
    % opt.DisplayRotation is only for how the plot is viewed on screen.
    % It does NOT change the model geometry, path length, or transforms.
    % Rplot is just a short local name for opt.DisplayRotation.

    p = inputParser;

    addParameter(p, 'FemurFile', 'Femur_Mesh_Points.xlsx');
    addParameter(p, 'TibiaFile', 'Tibia_Mesh_Points.xlsx');
    addParameter(p, 'UnitScale', 1);

    % Same display offset used in MuscleBonePlotting.
    addParameter(p, 'Hip', [-0.0707, -0.0661, 0.0835]);

    % Same display rotation used in MuscleBonePlotting.
    addParameter(p, 'DisplayRotation', [1, 0, 0;
                                        0, 0, 1;
                                        0,-1, 0]);

    addParameter(p, 'PauseTime', 0.25);  % slower
    addParameter(p, 'FrameStep', 1);
    addParameter(p, 'Loop', true);

    parse(p, varargin{:});
    opt = p.Results;
    Rplot = opt.DisplayRotation;

    %% Load raw bone point clouds
    % femur0 is in femur frame.
    % tibia0 is in tibia frame.
    %
    % Do NOT add Knee to tibia0 here.
    % T(:,:,i) already contains the tibia->femur pose, including the knee
    % translation and rotation for each frame.

    femur0 = readmatrix(opt.FemurFile) * opt.UnitScale;
    tibia0 = readmatrix(opt.TibiaFile) * opt.UnitScale;

    %% Convert optimized p2 into the human tibia/ICR frame at home pose
    % Since human and robot share the ICR at phi(pos)=0:
    %
    %   v2 = RowVecTrans(T_ICR_t1(:,:,pos), p2)
    %
    % For your pure translation T_ICR_t1, this is equivalent to:
    %
    %   v2 = p2 - t1toICR(1,:,pos)

    v2 = RowVecTrans(T_ICR_t1(:,:,pos), p2);

    %% Animation order
    % Start near home/extension, flex, then return.
    idxFlex = pos:-opt.FrameStep:1;
    idxReturn = 2:opt.FrameStep:pos;
    idxList = [idxFlex, idxReturn];

    %% Initial frame
    i = pos;
    Ti = T(:,:,i);

    [tibia_i_f, w2, humanPath_f] = transformFrame(Ti, tibia0, p1, v2, Bifemsh, i);

    femur_plot = (femur0 + opt.Hip) * Rplot;
    tibia_plot = (tibia_i_f + opt.Hip) * Rplot;
    bpa_plot   = ([p1; w2] + opt.Hip) * Rplot;
    human_plot = (humanPath_f + opt.Hip) * Rplot;

    %% Plot setup
    fig = figure;
    hold on
    grid on
    axis equal
    view(3)

    hFemur = plot3( ...
        femur_plot(:,1), femur_plot(:,2), femur_plot(:,3), ...
        '.', 'Color', [0.2 0.2 1.0]);

    hTibia = plot3( ...
        tibia_plot(:,1), tibia_plot(:,2), tibia_plot(:,3), ...
        '.', 'Color', [0.1 0.1 0.1]);

    hBPA = plot3( ...
        bpa_plot(:,1), bpa_plot(:,2), bpa_plot(:,3), ...
        '.-', 'Color', [0 0.6 0], 'LineWidth', 3, 'MarkerSize', 24);

    hHuman = plot3( ...
        human_plot(:,1), human_plot(:,2), human_plot(:,3), ...
        '.-', 'Color', [0.8 0 0], 'LineWidth', 3, 'MarkerSize', 20);

    hP1 = plot3( ...
        bpa_plot(1,1), bpa_plot(1,2), bpa_plot(1,3), ...
        'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

    hP2 = plot3( ...
        bpa_plot(2,1), bpa_plot(2,2), bpa_plot(2,3), ...
        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    legend('Femur', 'Tibia', 'Optimized BPA path', 'Human muscle path', ...
        'p1', 'p2', 'Location', 'best')

    allPts = [femur_plot; tibia_plot; bpa_plot; human_plot];
    pad = 0.08;
    xlim([min(allPts(:,1))-pad, max(allPts(:,1))+pad])
    ylim([min(allPts(:,2))-pad, max(allPts(:,2))+pad])
    zlim([min(allPts(:,3))-pad, max(allPts(:,3))+pad])

    %% Loop animation
    keepGoing = true;

    while keepGoing && ishandle(fig)

        for i = idxList

            if ~ishandle(fig)
                return
            end

            Ti = T(:,:,i);  % 4x4 only

            [tibia_i_f, w2, humanPath_f] = transformFrame(Ti, tibia0, p1, v2, Bifemsh, i);

            tibia_plot = (tibia_i_f + opt.Hip) * Rplot;
            bpa_plot   = ([p1; w2] + opt.Hip) * Rplot;
            human_plot = (humanPath_f + opt.Hip) * Rplot;

            set(hTibia, ...
                'XData', tibia_plot(:,1), ...
                'YData', tibia_plot(:,2), ...
                'ZData', tibia_plot(:,3));

            set(hBPA, ...
                'XData', bpa_plot(:,1), ...
                'YData', bpa_plot(:,2), ...
                'ZData', bpa_plot(:,3));

            set(hHuman, ...
                'XData', human_plot(:,1), ...
                'YData', human_plot(:,2), ...
                'ZData', human_plot(:,3));

            set(hP1, ...
                'XData', bpa_plot(1,1), ...
                'YData', bpa_plot(1,2), ...
                'ZData', bpa_plot(1,3));

            set(hP2, ...
                'XData', bpa_plot(2,1), ...
                'YData', bpa_plot(2,2), ...
                'ZData', bpa_plot(2,3));

            Lbpa = norm(p1 - w2);

            title(sprintf('Knee angle = %.1f deg, BPA path length = %.3f m', ...
                phi(i)*180/pi, Lbpa))

            drawnow
            pause(opt.PauseTime)
        end

        keepGoing = opt.Loop;
    end

end


function [tibia_i_f, w2, humanPath_f] = transformFrame(Ti, tibia0, p1, v2, Bifemsh, i)
% transformFrame
%
% Ti is T(:,:,i), a single 4x4 matrix.
%
% tibia0:
%   raw tibia mesh in tibia frame
%
% v2:
%   optimized insertion point in tibia frame
%
% w2:
%   v2 transformed into femur frame
%
% humanPath_f:
%   human muscle-tendon path transformed into femur frame

    %% Tibia mesh: tibia frame -> femur frame
    tibia_i_f = zeros(size(tibia0));

    for j = 1:size(tibia0,1)
        tibia_i_f(j,:) = RowVecTrans(Ti, tibia0(j,:));
    end

    %% Optimized insertion point: tibia frame -> femur frame
    w2 = RowVecTrans(Ti, v2);

    %% Human muscle path
    HLoc = Bifemsh.Location(:,:,i);
    HCross = Bifemsh.Cross;

    humanPath_f = zeros(size(HLoc));

    for j = 1:size(HLoc,1)

        if j < HCross
            % Femur-side human muscle path point.
            humanPath_f(j,:) = HLoc(j,:);
        else
            % Tibia-side human muscle path point.
            humanPath_f(j,:) = RowVecTrans(Ti, HLoc(j,:));
        end

    end

end