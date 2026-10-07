function data = readserialnumbers2()
%READSERIALNUMBERS Live Arduino acquisition and valve control.
%
% Click the live plot window and press:
%
%   V  Turn both valve outputs HIGH and begin recording
%   S  Save the recorded dynamic data and stop recording
%   O  Turn both valve outputs LOW
%   Q  Turn both outputs LOW and quit
%
% Live serial columns (encoder removed 2026-10-04):
%   1. Time (ms)
%   2. Force (N)
%   3. Pressure (kPa)
%   4. Fill valve output state
%   5. Exhaust valve output state
%
% Live torque guard (2026-10-04): enter the KNEE ANGLE and LOAD-CELL
% ANGLE in the boxes at the bottom of the window. The live force is
% converted to knee torque with the same Adjoint transform as the
% ExtTest20mm_1 section of Knee_Extensor_20mm.m (robot t1->ICR
% kinematics, reaction point dLC/angLC), and compared with the OpenSim
% human vasti torque target at that knee angle. The third plot shows the
% measured torque against the target line, and the status text reports
% whether the human torque magnitude is met (>=) and by what margin.
%
% Saved/returned columns add the guard values per row:
%   6. Knee angle as entered (deg)
%   7. Load-cell angle as entered (deg)
%   8. Measured knee torque (N*m, Adjoint)
%   9. Human torque target at that knee angle (N*m)
%
% run this to clear ports
  % ports = serialportfind("Port", "COM10");
  %
  % if ~isempty(ports)
  %   delete(ports);
  % end

    %% Serial settings

    port = "COM10";
    baudRate = 115200;

    expectedNumCols = 5;
    savedNumCols = 9;

    %% Torque-guard geometry (matches Knee_Extensor_20mm.m dExt1/angExt1)
    dLC  = 215.05/1000;    % theta1 origin -> load-cell arm, m
    angLC = -92.41;       % arm angle in the tibia frame, deg

    %% Save settings

    functionFolder = fileparts(mfilename("fullpath"));
    saveFolder = fullfile(functionFolder, "Ext_20mm");
    baseName = "ExtTest_2_";

    if ~isfolder(saveFolder)
        mkdir(saveFolder);
    end

    %% Load the torque-guard model (robot knee kinematics + human target)

    torqueOK = false;
    ctx = [];
    phiV = [];
    txV = [];
    tyV = [];
    try
        root = functionFolder;
        for k = 1:8
            [parent, name] = fileparts(root);
            if strcmpi(name, "Bipedal_Robot")
                break;
            end
            if strcmp(parent, root)
                error("Could not locate the Bipedal_Robot repo root.");
            end
            root = parent;
        end
        addpath(genpath(fullfile(root, "Code", "Matlab")));
        addpath(fullfile(root, "Code", "Matlab", "Mesh_Optimization"));
        addpath(fullfile(root, "Testing_Data", "2022_02_Festo"), "-end");

        % buildKneeExtContext20mm prints its route-seed geometry (tendon
        % lengths etc.) -- irrelevant here, so capture and discard it.
        ctx = loadTorqueGuardContext();
        pRF = [dLC*cosd(angLC), dLC*sind(angLC), 0];
        phiV = ctx.phi(:);
        txV = squeeze(ctx.T_t1_ICR(1, 4, :)) - pRF(1);
        tyV = squeeze(ctx.T_t1_ICR(2, 4, :)) - pRF(2);
        txV = txV(:);
        tyV = tyV(:);
        torqueOK = true;
        fprintf("Torque guard model loaded (knee kinematics + human vasti target).\n");
    catch torqueModelME
        fprintf("Torque guard unavailable (continuing without it):\n  %s\n", ...
            torqueModelME.message);
    end

    %% Live-display settings

    % Only this much recent data is retained for the rolling plots.
    % This is separate from recordedData.
    plotWindowSeconds = 60;

    % Limit figure redraw rate without limiting serial acquisition.
    plotUpdatePeriod = 0.10;

    % Prevent continuous serial traffic from blocking figure callbacks
    maxLinesPerPass = 20;

    %% Open the serial port

    % Free the port if a stale MATLAB handle from an earlier crashed run
    % still holds it (the manual snippet from the header, done for you).
    stale = serialportfind("Port", port);
    if ~isempty(stale)
        fprintf("Deleting stale MATLAB handle(s) on %s...\n", char(port));
        delete(stale);
        pause(1);
    end

    try
        s = serialport(port, baudRate);
    catch portME
        avail = string(serialportlist("available"));
        if ~any(strcmpi(avail, port))
            error("readserialnumbers2:PortAbsent", ...
                "Port %s is not visible to MATLAB. Available: %s.\nPlug the board in (or unplug/replug it), then retry.", ...
                char(port), strjoin(avail, ", "));
        else
            error("readserialnumbers2:PortBusy", ...
                "Port %s is visible but BUSY -- close the Arduino IDE Serial Monitor\n(it holds the port after a sketch upload), wait a second, and retry.\n(%s)", ...
                char(port), portME.message);
        end
    end
    configureTerminator(s, "LF");
    s.Timeout = 1;

    % Opening a serial port normally resets this Arduino-compatible board.
    pause(3);

    % Remove startup text or incomplete lines.
    flush(s, "input");

    %% Data storage

    % recordedData contains only data collected between V and S
    % (serial columns + entered angles + torque guard columns).
    recordedData = zeros(0, savedNumCols);

    % liveData contains only the recent rolling display window.
    liveData = zeros(0, savedNumCols);

    % Function output. This is populated when S is pressed.
    data = zeros(0, savedNumCols);

    recording = false;
    recordingPending = false;
    stopRequested = false;
    testSaved = false;
    warnedSixCol = false;

    fillState = 0;
    exhaustState = 0;

    %% Create the live display

    % Explicit on-screen position: the MATLAB default figure position can
    % land off-screen (display scaling / stale monitor layouts), which
    % shows as a taskbar entry that never becomes a visible window.
    ss = get(groot, "ScreenSize");
    figW = 780;
    figH = 660;
    figPos = [max(20, floor((ss(3) - figW)/2)), ...
              max(40, floor((ss(4) - figH)/2)), figW, figH];

    fig = figure( ...
        "Name", "Arduino Live Data", ...
        "NumberTitle", "off", ...
        "Position", figPos, ...
        "WindowKeyReleaseFcn", @keyReleased, ...
        "CloseRequestFcn", @closeRequested);
    movegui(fig, "center");

    layout = tiledlayout(fig, 3, 1, ...
        "TileSpacing", "compact", ...
        "Padding", "compact");

    axForce = nexttile(layout);
    forceLine = plot(axForce, NaN, NaN);
    ylabel(axForce, "Force (N)");
    grid(axForce, "on");

    axPressure = nexttile(layout);
    pressureLine = plot(axPressure, NaN, NaN);
    ylabel(axPressure, "Pressure (kPa)");
    grid(axPressure, "on");

    axTorque = nexttile(layout);
    torqueLine = plot(axTorque, NaN, NaN, "LineWidth", 1.5);
    hold(axTorque, "on");
    targetLine = plot(axTorque, NaN, NaN, "r--", "LineWidth", 1.5);
    hold(axTorque, "off");
    ylabel(axTorque, "Knee torque (N\cdotm)");
    xlabel(axTorque, "Time relative to latest sample (s)");
    grid(axTorque, "on");

    statusTitle = title(layout, ...
        "V = valves on/start | S = save | O = valves off | Q = quit");

    % Knee/load-cell angle entry + live torque verdict (bottom strip)
    uKneeLabel = uicontrol(fig, "Style", "text", ...
        "String", "Knee angle [deg]", "Units", "normalized", ...
        "Position", [0.01 0.055 0.13 0.030], "HorizontalAlignment", "left"); %#ok<NASGU>
    uKneeEdit = uicontrol(fig, "Style", "edit", "String", "-30", ...
        "Units", "normalized", "Position", [0.01 0.010 0.07 0.040]);
    uLCLabel = uicontrol(fig, "Style", "text", ...
        "String", "LC angle [deg, from tibia axis]", "Units", "normalized", ...
        "Position", [0.15 0.055 0.22 0.030], "HorizontalAlignment", "left"); %#ok<NASGU>
    uLCEdit = uicontrol(fig, "Style", "edit", "String", "30", ...
        "Units", "normalized", "Position", [0.15 0.010 0.07 0.040]);
    uVerdict = uicontrol(fig, "Style", "text", "String", ...
        "Enter knee + load-cell angles for the torque guard", ...
        "Units", "normalized", "Position", [0.35 0.010 0.63 0.045], ...
        "HorizontalAlignment", "left", "FontWeight", "bold");

    % Ensure the valves are switched off when the function exits,
    % including exits caused by an error.
    cleanupObject = onCleanup(@()shutdownSerial(s, fig)); %#ok<NASGU>

    fprintf("\nArduino acquisition started.\n");
    fprintf("Click the live plot and release each key once:\n");
    fprintf("  V = valves HIGH and begin recording\n");
    fprintf("  S = save dynamic data and stop recording\n");
    fprintf("  O = valves LOW\n");
    fprintf("  Q = valves LOW and quit\n\n");

    % Request the current valve state after the startup buffer was cleared.
    write(s, uint8('?'), "uint8");

    %% Continuous acquisition loop

    lastPlotUpdate = tic;

    while ~stopRequested && isgraphics(fig)

        % Process figure keyboard callbacks.
        drawnow;

        %% Read a limited number of complete serial lines per pass
        linesRead = 0;

        while s.NumBytesAvailable > 0 && linesRead < maxLinesPerPass
            linesRead = linesRead + 1;

            try
                lineText = strtrim(readline(s));
            catch
                % A partial line may occasionally time out.
                continue;
            end

            if strlength(lineText) == 0
                continue;
            end

            %% Process board-startup messages

            if startsWith(lineText, "BOOT,")
                parts = split(lineText, ",");

                if numel(parts) == 3
                    stateValues = str2double(parts(2:3));

                    if all(isfinite(stateValues))
                        fillState = stateValues(1);
                        exhaustState = stateValues(2);

                        fprintf( ...
                            "Arduino boot state: fill = %d, exhaust = %d\n", ...
                            fillState, exhaustState);
                    end
                end

                continue;
            end

            %% Process immediate valve-state messages

            if startsWith(lineText, "STATE,")
                parts = split(lineText, ",");

                if numel(parts) == 3
                    stateValues = str2double(parts(2:3));

                    if all(isfinite(stateValues))
                        fillState = stateValues(1);
                        exhaustState = stateValues(2);

                        fprintf( ...
                            "Arduino valve state: fill = %d, exhaust = %d\n", ...
                            fillState, exhaustState);
                    end
                end

                continue;
            end

            %% Process command acknowledgements

            if startsWith(lineText, "ACK,")
                parts = split(lineText, ",");

                if numel(parts) == 4
                    acknowledgedCommand = upper(parts(2));
                    stateValues = str2double(parts(3:4));

                    if all(isfinite(stateValues))
                        fillState = stateValues(1);
                        exhaustState = stateValues(2);

                        fprintf( ...
                            "Arduino acknowledged %s: fill = %d, exhaust = %d\n", ...
                            char(acknowledgedCommand), fillState, exhaustState);

                        if acknowledgedCommand == "V"
                            if fillState == 1 && exhaustState == 1
                                recordingPending = false;
                                recording = true;
                                fprintf("Valve HIGH state confirmed. Recording started.\n");
                            else
                                recordingPending = false;
                                recording = false;
                                fprintf( ...
                                    "Valve command failed: Arduino did not report both outputs HIGH.\n");
                            end

                        elseif acknowledgedCommand == "O"
                            recordingPending = false;
                            recording = false;
                        end
                    end
                end

                continue;
            end

            %% Process numeric measurement lines
            % Accept BOTH firmware formats: the new 5-column sketch
            % (time, force, pressure, fill, exhaust) and the old
            % 6-column encoder sketch (angle in position 2, ignored).

            parts = split(lineText, ",");

            switch numel(parts)
                case expectedNumCols
                    numericValues = str2double(parts).';
                case 6
                    if ~warnedSixCol
                        warnedSixCol = true;
                        fprintf(['Board is streaming the OLD 6-column (encoder) ', ...
                            'format -- angle column ignored.\nRe-upload ', ...
                            'ValveDataAcquisition.ino to switch to 5 columns.\n']);
                    end
                    numericValues = str2double(parts([1 3 4 5 6])).';
                otherwise
                    continue;
            end

            if any(~isfinite(numericValues))
                continue;
            end

            fillState = numericValues(4);
            exhaustState = numericValues(5);

            %% Torque guard for this row (entered angles, live force)

            Kdeg = str2double(strtrim(uKneeEdit.String));
            Ldeg = str2double(strtrim(uLCEdit.String));
            tz = computeTorqueZ(Kdeg, Ldeg, numericValues(2));
            tgt = humanTargetAt(Kdeg);

            guardRow = [Kdeg, Ldeg, tz, tgt];

            %% Start recording with the first confirmed valves-HIGH row
            % This is a fallback in case an ACK line was missed.

            if recordingPending && ...
                    fillState == 1 && exhaustState == 1

                recordingPending = false;
                recording = true;

                fprintf("Valve HIGH state confirmed. Recording started.\n");
            end

            %% Store dynamic data only while recording

            if recording
                recordedData(end + 1, :) = [numericValues, guardRow]; %#ok<AGROW>
            end

            %% Store only a limited rolling window for the plots

            liveData(end + 1, :) = [numericValues, guardRow]; %#ok<AGROW>

            latestTimeMs = numericValues(1);
            cutoffTimeMs = latestTimeMs - 1000 * plotWindowSeconds;

            liveData(liveData(:, 1) < cutoffTimeMs, :) = [];
        end

        %% Update the plots without slowing serial collection

        if ~isempty(liveData) && ...
                toc(lastPlotUpdate) >= plotUpdatePeriod

            % Show time relative to the most recent measurement.
            relativeTime = ...
                (liveData(:, 1) - liveData(end, 1)) / 1000;

            set(forceLine, ...
                "XData", relativeTime, ...
                "YData", liveData(:, 2));

            set(pressureLine, ...
                "XData", relativeTime, ...
                "YData", liveData(:, 3));

            set(torqueLine, ...
                "XData", relativeTime, ...
                "YData", liveData(:, 8));

            tgtNow = liveData(end, 9);
            set(targetLine, ...
                "XData", relativeTime([1 end]), ...
                "YData", [tgtNow tgtNow]);

            if recording
                recordingText = sprintf( ...
                    "RECORDING: %d rows", size(recordedData, 1));
            elseif recordingPending
                recordingText = "WAITING FOR VALVE CONFIRMATION";
            else
                recordingText = "NOT RECORDING";
            end

            statusTitle.String = sprintf( ...
                '%s | Fill = %d | Exhaust = %d\nForce = %.3f N | Pressure = %.3f kPa', ...
                char(recordingText), ...
                fillState, ...
                exhaustState, ...
                liveData(end, 2), ...
                liveData(end, 3));

            updateVerdict(liveData(end, 8), liveData(end, 9), ...
                liveData(end, 7), liveData(end, 2));

            drawnow limitrate;
            lastPlotUpdate = tic;

        else
            pause(0.005);
        end
    end

    %% Return unsaved recorded data if the function was quit before S

    if isempty(data) && ~isempty(recordedData)
        data = recordedData;

        % Make saved/returned time begin at zero.
        data(:, 1) = data(:, 1) - data(1, 1);
    end

    %% Nested helper: Adjoint knee torque from live force + entered angles

    function tz = computeTorqueZ(Kdeg, Ldeg, forceN)
        tz = NaN;
        if ~torqueOK
            return;
        end
        if ~isfinite(Kdeg) || ~isfinite(Ldeg) || ~isfinite(forceN)
            return;
        end
        Kr = deg2rad(Kdeg);
        % Ben's angle convention (2026-10-04): the LC angle is measured
        % FROM THE TIBIA AXIS -- torque about t1 = F*sin(LC+0.83 deg)*d.
        % The Adjoint machinery expects the force angle from the tibia
        % x-axis, i.e. 90 - LC. Convert here.
        Lr = deg2rad(Ldeg);
        tx = interp1(phiV, txV, Kr, "pchip");
        ty = interp1(phiV, tyV, Kr, "pchip");
        if ~isfinite(tx) || ~isfinite(ty)
            return;   % knee angle outside the robot kinematics table
        end
        Trk = RpToTrans(eye(3), [tx; ty; 0]);
        Fr = -[0; 0; 0; forceN*cos(pi - Lr); forceN*sin(pi - Lr); 0];
        Fk = Adjoint(Trk)' * Fr;
        tz = Fk(3);
    end

    function tgt = humanTargetAt(Kdeg)
        tgt = NaN;
        if torqueOK && isfinite(Kdeg)
            tgt = interp1(ctx.humanAngleD(:), ctx.humanTorque(:), ...
                Kdeg, "pchip");
        end
    end

    function updateVerdict(tz, tgt, Ldeg, forceN)
        if ~torqueOK
            uVerdict.String = "Torque guard unavailable (model load failed)";
            uVerdict.ForegroundColor = [0.5 0.5 0.5];
            return;
        end
        if ~isfinite(tz) || ~isfinite(tgt)
            uVerdict.String = sprintf( ...
                "Enter valid knee + load-cell angles (LC %.1f deg, force %.2f N)", ...
                Ldeg, forceN);
            uVerdict.ForegroundColor = [0 0 0];
            return;
        end
        margin = tz - tgt;
        if margin >= 0
            uVerdict.String = sprintf( ...
                "MET: %.2f N*m vs human %.2f N*m (+%.2f margin)", ...
                tz, tgt, margin);
            uVerdict.ForegroundColor = [0 0.55 0];
        else
            uVerdict.String = sprintf( ...
                "BELOW: %.2f N*m vs human %.2f N*m (%.2f short)", ...
                tz, tgt, margin);
            uVerdict.ForegroundColor = [0.8 0 0];
        end
    end

    %% Nested callback functions

    function keyReleased(~, event)

        key = lower(string(event.Key));

        switch key

            case "v"

                % Ignore duplicate V commands while already starting/recording
                if recordingPending || recording
                    fprintf("A test is already recording.\n");
                    return;
                end

                % Begin a new dynamic recording
                recordedData = zeros(0, savedNumCols);
                data = zeros(0, savedNumCols);

                recording = false;
                recordingPending = true;
                testSaved = false;

                % Remove any queued startup/state lines, then send exactly one byte.
                flush(s, "input");
                write(s, uint8('V'), "uint8");

                fprintf( ...
                    "Sent V: requesting both valve outputs HIGH.\n");

            case "s"

                % Prevent saving the same recording more than once
                if testSaved
                    fprintf([ ...
                        "This recording has already been saved.\n" ...
                        "Press V to begin a new test.\n"]);
                    return;
                end

                if isempty(recordedData)
                    fprintf( ...
                        "Nothing was saved because no recorded rows exist.\n");
                    return;
                end

                % Stop adding rows after this point
                recording = false;
                recordingPending = false;

                data = recordedData;

                % Make the first recorded sample time equal to zero
                data(:, 1) = data(:, 1) - data(1, 1);

                [matFileName, csvFileName] = getNextFileNames();

                save(matFileName, "data");
                writematrix(data, csvFileName);

                % Lock this recording against additional S commands
                testSaved = true;

                fprintf("\nDynamic data saved as:\n");
                fprintf("%s\n", char(matFileName));
                fprintf("%s\n\n", char(csvFileName));

                fprintf([ ...
                    "Live readings will continue, but new rows " ...
                    "are not being stored.\n"]);

            case "o"

                write(s, uint8('O'), "uint8");

                % Stop recording when the valves are switched off
                recording = false;
                recordingPending = false;

                fprintf( ...
                    "Sent O: requesting both valve outputs LOW.\n");

            case "q"

                try
                    write(s, uint8('O'), "uint8");
                catch
                end

                stopRequested = true;

            otherwise
                % Ignore all other keys
        end
    end


    function closeRequested(source, ~)

        % Hide immediately, then allow the acquisition loop to exit.
        source.Visible = "off";
        stopRequested = true;
    end


    function [matFileName, csvFileName] = getNextFileNames()

        testNumber = 1;

        while true
            testLabel = sprintf("%02d", testNumber);
            fileStem = baseName + string(testLabel);

            matFileName = fullfile( ...
                saveFolder, fileStem + ".mat");

            csvFileName = fullfile( ...
                saveFolder, fileStem + ".csv");

            if ~isfile(matFileName) && ~isfile(csvFileName)
                break;
            end

            testNumber = testNumber + 1;
        end
    end
end


function shutdownSerial(s, fig)
%SHUTDOWNSERIAL Force valve outputs LOW and release the serial port.

    try
        write(s, uint8('O'), "uint8");
        pause(0.05);
    catch
    end

    try
        delete(s);
    catch
    end

    if isgraphics(fig)
        delete(fig);
    end
end


function ctx = loadTorqueGuardContext()
%LOADTORQUEGUARDCONTEXT Build the extensor context with its chatter
% suppressed (evalc is not allowed in the main function because it
% contains nested functions, so the call lives here).

    evalc("ctx = buildKneeExtContext20mm();");
end
