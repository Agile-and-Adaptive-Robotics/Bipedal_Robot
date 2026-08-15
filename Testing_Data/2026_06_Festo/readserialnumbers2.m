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
% Returned/saved columns:
%   1. Time (ms)
%   2. Position (degrees)
%   3. Force (N)
%   4. Pressure (kPa)
%   5. Fill valve output state
%   6. Exhaust valve output state
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

    expectedNumCols = 6;

    %% Save settings

    functionFolder = fileparts(mfilename("fullpath"));
    saveFolder = fullfile(functionFolder, "Flx_10mm_pinned");
    baseName = "FlxTest06_";

    if ~isfolder(saveFolder)
        mkdir(saveFolder);
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

    s = serialport(port, baudRate);
    configureTerminator(s, "LF");
    s.Timeout = 1;

    % Opening a serial port normally resets this Arduino-compatible board.
    pause(3);

    % Remove startup text or incomplete lines.
    flush(s, "input");

    %% Data storage

    % recordedData contains only data collected between V and S.
    recordedData = zeros(0, expectedNumCols);

    % liveData contains only the recent rolling display window.
    liveData = zeros(0, expectedNumCols);

    % Function output. This is populated when S is pressed.
    data = zeros(0, expectedNumCols);

    recording = false;
    recordingPending = false;
    stopRequested = false;
    testSaved = false;

    fillState = 0;
    exhaustState = 0;

    %% Create the live display

    fig = figure( ...
        "Name", "Arduino Live Data", ...
        "NumberTitle", "off", ...
        "WindowKeyReleaseFcn", @keyReleased, ...
        "CloseRequestFcn", @closeRequested);

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

    axAngle = nexttile(layout);
    angleLine = plot(axAngle, NaN, NaN);
    ylabel(axAngle, "Position (deg)");
    xlabel(axAngle, "Time relative to latest sample (s)");
    grid(axAngle, "on");

    statusTitle = title(layout, ...
        "V = valves on/start | S = save | O = valves off | Q = quit");

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

            parts = split(lineText, ",");

            if numel(parts) ~= expectedNumCols
                continue;
            end

            numericValues = str2double(parts).';

            if any(~isfinite(numericValues))
                continue;
            end

            fillState = numericValues(5);
            exhaustState = numericValues(6);

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
                recordedData(end + 1, :) = numericValues; %#ok<AGROW>
            end

            %% Store only a limited rolling window for the plots

            liveData(end + 1, :) = numericValues; %#ok<AGROW>

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
                "YData", liveData(:, 3));

            set(pressureLine, ...
                "XData", relativeTime, ...
                "YData", liveData(:, 4));

            set(angleLine, ...
                "XData", relativeTime, ...
                "YData", liveData(:, 2));

            if recording
                recordingText = sprintf( ...
                    "RECORDING: %d rows", size(recordedData, 1));
            elseif recordingPending
                recordingText = "WAITING FOR VALVE CONFIRMATION";
            else
                recordingText = "NOT RECORDING";
            end

            statusTitle.String = sprintf( ...
                '%s | Fill = %d | Exhaust = %d\nAngle = %.3f deg | Force = %.3f N | Pressure = %.3f kPa', ...
                char(recordingText), ...
                fillState, ...
                exhaustState, ...
                liveData(end, 2), ...
                liveData(end, 3), ...
                liveData(end, 4));

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
                recordedData = zeros(0, expectedNumCols);
                data = zeros(0, expectedNumCols);

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
