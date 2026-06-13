classdef HX711 < matlab.apps.AppBase

    % HX711 data-collection app for bipedal robot flexion tests.
    % MATLAB R2025-compatible programmatic App Designer class.
    %
    % Save convention:
    %   <repo>\Testing_Data\2026_06_Festo\Flx_20mm\FlxTest$_#.mat
    % where $ = test series selected by user, # = run number 00..99.
    %
    % MAT-file variables:
    %   Data  - 750 x A numeric array. Columns are listed in ColumnNames.
    %   Stats - table with rows Mean, Median, Mode, Min, Max, StdDev and
    %           columns Force, Pressure.
    %   Metadata - struct with angles, pins, calibration coefficients, units.
    %   ColumnNames - names for Data columns.

    % Properties that correspond to app components
    properties (Access = public)
        MatlabArduinoHX711UIFigure     matlab.ui.Figure
        PressureGauge                  matlab.ui.control.SemicircularGauge
        PressureGaugeLabel             matlab.ui.control.Label
        kPaLabel                       matlab.ui.control.Label
        Pressure                       matlab.ui.control.NumericEditField
        PressureLabel                  matlab.ui.control.Label
        CleanPanel                     matlab.ui.container.Panel
        Clean                          matlab.ui.control.Button
        Rate                           matlab.ui.control.NumericEditField
        SamplingRLabel                 matlab.ui.control.Label
        measure3                       matlab.ui.control.EditField
        measure2                       matlab.ui.control.EditField
        measure                        matlab.ui.control.EditField
        DataEditField                  matlab.ui.control.NumericEditField
        DataEditFieldLabel             matlab.ui.control.Label
        TimeEdit                       matlab.ui.control.NumericEditField
        TimeEditFieldLabel             matlab.ui.control.Label
        ForceEdit                      matlab.ui.control.NumericEditField
        ForceEditFieldLabel            matlab.ui.control.Label
        ContinuosDataAcquisitionPanel  matlab.ui.container.Panel
        Axes1                          matlab.ui.control.UIAxes
        GlobalSettingsPanel            matlab.ui.container.Panel
        TabGroup                       matlab.ui.container.TabGroup
        ConnectionTab                  matlab.ui.container.Tab
        PressureEdit                   matlab.ui.control.EditField
        PressurePinLabel               matlab.ui.control.Label
        DataEdit                       matlab.ui.control.EditField
        DataPinEditFieldLabel          matlab.ui.control.Label
        BoardEdit                      matlab.ui.control.DropDown
        ArduinoDropDownLabel           matlab.ui.control.Label
        ClockEdit                      matlab.ui.control.EditField
        ClockPinEditFieldLabel         matlab.ui.control.Label
        SerialEdit                     matlab.ui.control.EditField
        SerialportEditFieldLabel       matlab.ui.control.Label
        ValveIncEdit                   matlab.ui.control.EditField
        ValveMaintainEdit              matlab.ui.control.EditField
        ValveIncPinLabel               matlab.ui.control.Label
        ValveMaintainPinLabel          matlab.ui.control.Label
        DataAcquisitionTab             matlab.ui.container.Tab
        Add_time                       matlab.ui.control.Spinner
        SamplingRateSpinnerLabel       matlab.ui.control.Label
        SessionV                       matlab.ui.control.Spinner
        SessionTimeSpinnerLabel        matlab.ui.control.Label
        SetSession                     matlab.ui.control.CheckBox
        Unit                           matlab.ui.control.DropDown
        ForceDropDownLabel             matlab.ui.control.Label
        Nsamples                       matlab.ui.control.Spinner
        SampleCountSpinnerLabel        matlab.ui.control.Label
        SaveDataTab                    matlab.ui.container.Tab
        Reset                          matlab.ui.control.CheckBox
        Name                           matlab.ui.control.EditField
        NameEditFieldLabel             matlab.ui.control.Label
        SaveFolder                     matlab.ui.control.EditField
        SaveFolderLabel                matlab.ui.control.Label
        SeriesSpinner                  matlab.ui.control.Spinner
        SeriesSpinnerLabel             matlab.ui.control.Label
        RunSpinner                     matlab.ui.control.Spinner
        RunSpinnerLabel                matlab.ui.control.Label
        MetadataTab                    matlab.ui.container.Tab
        KneeAngle                      matlab.ui.control.NumericEditField
        KneeAngleLabel                 matlab.ui.control.Label
        LoadCellAngle                  matlab.ui.control.NumericEditField
        LoadCellAngleLabel             matlab.ui.control.Label
        CalibrationResultPanel         matlab.ui.container.Panel
        TabGroup2                      matlab.ui.container.TabGroup
        ValueTab                       matlab.ui.container.Tab
        measure4                       matlab.ui.control.EditField
        Raw                            matlab.ui.control.NumericEditField
        RawreadingLabel                matlab.ui.control.Label
        StdDisp                        matlab.ui.control.NumericEditField
        StdDeviationgLabel             matlab.ui.control.Label
        AverageDisp                    matlab.ui.control.NumericEditField
        AveragegLabel                  matlab.ui.control.Label
        ScaleDisp                      matlab.ui.control.NumericEditField
        ScalefactorLabel               matlab.ui.control.Label
        TareDisp                       matlab.ui.control.NumericEditField
        TareLabel                      matlab.ui.control.Label
        Gauge                          matlab.ui.control.LinearGauge
        SettingsTab                    matlab.ui.container.Tab
        UnitLabel_2                    matlab.ui.control.Label
        UnitLabel                      matlab.ui.control.Label
        UnitCal2                       matlab.ui.control.DropDown
        MaxLoad                        matlab.ui.control.NumericEditField
        MaxLoadLabel                   matlab.ui.control.Label
        UnitCal                        matlab.ui.control.DropDown
        Known                          matlab.ui.control.NumericEditField
        KnownWeightLabel               matlab.ui.control.Label
        n                              matlab.ui.control.NumericEditField
        NumreadingsLabel               matlab.ui.control.Label
        KnownCalTab                    matlab.ui.container.Tab
        ApplyKnownLoadCellButton       matlab.ui.control.Button
        KnownTare                      matlab.ui.control.NumericEditField
        KnownTareLabel                 matlab.ui.control.Label
        KnownScale                     matlab.ui.control.NumericEditField
        KnownScaleLabel                matlab.ui.control.Label
        PressureCalTab                 matlab.ui.container.Tab
        PressureA                      matlab.ui.control.NumericEditField
        PressureALabel                 matlab.ui.control.Label
        PressureB                      matlab.ui.control.NumericEditField
        PressureBLabel                 matlab.ui.control.Label
        PressureCalN                   matlab.ui.control.NumericEditField
        PressureCalNLabel              matlab.ui.control.Label
        ApplyKnownPressureButton       matlab.ui.control.Button
        PressureCalButton              matlab.ui.control.Button
        LicenseTab                     matlab.ui.container.Tab
        TextArea                       matlab.ui.control.TextArea
        Axes2                          matlab.ui.control.UIAxes
        CalibrationPanel               matlab.ui.container.Panel
        RawRead                        matlab.ui.control.Button
        CalibrationButton              matlab.ui.control.Button
        ScaleFactorButton              matlab.ui.control.Button
        TareButton                     matlab.ui.control.Button
        ArduinoHX711Panel              matlab.ui.container.Panel
        SaveButton                     matlab.ui.control.Button
        PauseButton                    matlab.ui.control.Button
        GetData                        matlab.ui.control.Button
        Connect                        matlab.ui.control.Button
        StatusPanel                    matlab.ui.container.Panel
        Cyan                           matlab.ui.control.Lamp
        DataAcquisitionLabel           matlab.ui.control.Label
        Yellow                         matlab.ui.control.Lamp
        InpauseLabel                   matlab.ui.control.Label
        Red                            matlab.ui.control.Lamp
        NotConnectedLabel              matlab.ui.control.Label
        Green                          matlab.ui.control.Lamp
        ConnectedLabel                 matlab.ui.control.Label
        Message                        matlab.ui.control.EditField
        MessageEditFieldLabel          matlab.ui.control.Label
        ValvePanel                     matlab.ui.container.Panel
        IncreasePressureButton         matlab.ui.control.Button
        MaintainPressureButton         matlab.ui.control.Button
        DecreasePressureButton         matlab.ui.control.Button
    end

    properties (Access = private)
        a              % Arduino object
        serial         % Serial port
        board          % Arduino board
        data           % Data pin
        clock          % Clock pin
        pressurePin    % Pressure pin
        valveIncPin    % Valve pin D11 default
        valveMaintainPin % Valve pin D6 default
        HX711_obj      % HX711 class object
        tare = NaN     % zero offset/raw tare; NaN means not calibrated yet
        scale = NaN    % raw-counts per gram equivalent; NaN means not calibrated yet
        pressureA = 155.61 % kPa/V slope, default from old app
        pressureB = -126.99 % kPa intercept, default from old app
        g = 9.80665
        check_connection = false
        v = 1
        t = 0
        session_time = Inf
        known_weight = 0
        Max_Load
        last_time1 = 0
        last_time2 = 0
        last_data = 0
        file_name1 = ""
        file_name2 = ""
        rst_time = 0
        get_true = false
        start_time = []
    end

    properties (Access = public)
        i = 1           % Index of force/time/pressure arrays
        time = []       % Array of time [s]
        force = []      % Array of force in currently selected display units
        pressure = []   % Array of pressure [kPa]
        rawForce = []   % Raw HX711 readings
        forceN = []     % Force converted to Newtons for saving/stats
        pressureV = []  % Pressure sensor voltage [V]
        r = 0
        k1 = 1
        k2 = 1
    end

    methods (Access = private)

        function Xaxis(app)
            if app.time(app.i) > 100*app.v || app.t == 1
                app.v = app.v + 1;
                app.Axes1.XLim = [0 100*app.v];
                app.Axes1.XTick = [0 10*app.v 20*app.v 30*app.v 40*app.v 50*app.v ...
                    60*app.v 70*app.v 80*app.v 90*app.v 100*app.v];
            end
        end

        function Xaxis2(app)
            app.Axes1.XLim = [0 100];
            app.Axes1.XTick = [0 10 20 30 40 50 60 70 80 90 100];
        end

        function grams = rawToGrams(app, rawValue)
            grams = (double(rawValue) - app.tare)./app.scale;
        end

        function value = gramsToSelectedUnit(app, grams)
            switch app.Unit.Value
                case '[ g ]'
                    value = grams;
                case '[ kg ]'
                    value = grams/1000;
                case '[ N ]'
                    value = grams*app.g/1000;
                case '[ kN ]'
                    value = grams*app.g/1000000;
                otherwise
                    value = grams*app.g/1000;
            end
        end

        function unitText = selectedUnitText(app)
            switch app.Unit.Value
                case '[ kg ]'
                    unitText = 'kg';
                case '[ N ]'
                    unitText = 'N';
                case '[ g ]'
                    unitText = 'g';
                case '[ kN ]'
                    unitText = 'kN';
                otherwise
                    unitText = 'N';
            end
        end

        function forceN = gramsToNewtons(app, grams)
            forceN = grams*app.g/1000;
        end

        function kPa = pressureVoltageToKPa(app, voltage)
            kPa = app.pressureA*double(voltage) + app.pressureB;
        end

        function tf = isLoadCellCalibrated(app)
            tf = isfinite(app.tare) && isfinite(app.scale) && app.scale ~= 0;
        end

        function updateStatusConnected(app, connected)
            app.check_connection = connected;
            if connected
                app.Red.Color = 'white';
                app.Green.Color = 'green';
                app.Yellow.Color = 'yellow';
            else
                app.Red.Color = 'red';
                app.Green.Color = 'white';
                app.Yellow.Color = 'white';
                app.Cyan.Color = 'white';
            end
        end

        function setValves(app, d11State, d6State, msg)
            if ~app.check_connection
                app.Message.Value = 'Error: You are not connected yet.';
                return;
            end
            writeDigitalPin(app.a, app.valveIncPin, logical(d11State));
            writeDigitalPin(app.a, app.valveMaintainPin, logical(d6State));
            app.Message.Value = msg;
        end

        function forceTitle(app)
            ylabel(app.Axes1, ['Force ', app.Unit.Value]);
        end

        function stats = buildStatsTable(~, forceN, pressureKPa)
            forceN = forceN(:);
            pressureKPa = pressureKPa(:);
            rowNames = {'Mean';'Median';'Mode';'Min';'Max';'StdDev'};
            Force = [mean(forceN,'omitnan'); median(forceN,'omitnan'); mode(forceN); ...
                min(forceN,[],'omitnan'); max(forceN,[],'omitnan'); std(forceN,0,'omitnan')];
            Pressure = [mean(pressureKPa,'omitnan'); median(pressureKPa,'omitnan'); mode(pressureKPa); ...
                min(pressureKPa,[],'omitnan'); max(pressureKPa,[],'omitnan'); std(pressureKPa,0,'omitnan')];
            stats = table(Force, Pressure, 'RowNames', rowNames);
        end

        function filePath = buildSavePath(app)
            saveDir = string(app.SaveFolder.Value);
            if strlength(saveDir) == 0
                saveDir = fullfile('GitHub','Bipedal_Robot','Testing_Data','2026_06_Festo','Flx_20mm');
            end
            if ~isfolder(saveDir)
                mkdir(saveDir);
            end
            fileName = sprintf('FlxTest%d_%02d.mat', round(app.SeriesSpinner.Value), round(app.RunSpinner.Value));
            filePath = fullfile(saveDir, fileName);
        end

        function readPressureOnce(app)
            app.pressureV(end+1) = readVoltage(app.a, app.pressurePin);
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: Connect
        function ConnectButtonPushed(app, event)
            app.Message.Value = 'Please wait...';
            drawnow;
            try
                app.serial = app.SerialEdit.Value;
                app.board = app.BoardEdit.Value;
                app.data = app.DataEdit.Value;
                app.clock = app.ClockEdit.Value;
                app.pressurePin = app.PressureEdit.Value;
                app.valveIncPin = app.ValveIncEdit.Value;
                app.valveMaintainPin = app.ValveMaintainEdit.Value;

                app.a = arduino(app.serial, app.board, 'libraries', {'basicHX711/basic_HX711'});
                app.HX711_obj = addon(app.a, 'basicHX711/basic_HX711', {app.data, app.clock});
                configurePin(app.a, app.pressurePin, 'AnalogInput');
                configurePin(app.a, app.valveIncPin, 'DigitalOutput');
                configurePin(app.a, app.valveMaintainPin, 'DigitalOutput');
                writeDigitalPin(app.a, app.valveIncPin, 0);
                writeDigitalPin(app.a, app.valveMaintainPin, 0);
                app.pressureV = readVoltage(app.a, app.pressurePin);
                updateStatusConnected(app, true);
                app.Message.Value = 'Connected.';
            catch ME
                updateStatusConnected(app, false);
                app.Message.Value = ['Connection error: ', ME.message];
            end
        end

        % Button pushed function: GetData
        function GetDataButtonPushed(app, event)
            if ~app.check_connection
                app.Message.Value = 'Error: You are not connected yet.';
                return;
            end
            if ~isLoadCellCalibrated(app)
                app.Message.Value = 'Error: enter/apply load-cell tare and scale, or calibrate first.';
                return;
            end

            app.r = 0;
            app.t = 0;
            app.measure.Value = selectedUnitText(app);
            forceTitle(app);

            sampleTarget = round(app.Nsamples.Value);
            if isempty(app.time) || app.i == 1
                app.time(app.i) = 0;
                app.start_time = tic;
            else
                app.time(app.i) = app.time(app.i-1);
            end

            if app.SetSession.Value == 1
                app.session_time = app.SessionV.Value*60 + app.time(app.i);
            else
                app.session_time = Inf;
            end

            app.Message.Value = 'Data acquisition...';
            app.Cyan.Color = 'cyan';
            app.Yellow.Color = 'white';
            app.get_true = true;
            drawnow;

            while (app.r == 0) && (app.time(app.i) < app.session_time) && ((app.i - app.last_data) < sampleTarget)
                loopTimer = tic;
                raw = read_HX711(app.HX711_obj);
                grams = rawToGrams(app, raw);
                forceSelected = gramsToSelectedUnit(app, grams);
                forceNValue = gramsToNewtons(app, grams);

                voltage = readVoltage(app.a, app.pressurePin);
                p = pressureVoltageToKPa(app, voltage);

                if app.MaxLoad.Value ~= 0
                    switch app.UnitCal2.Value
                        case '[ g ]'
                            app.Max_Load = app.MaxLoad.Value;
                        case '[ kg ]'
                            app.Max_Load = app.MaxLoad.Value*1000;
                        case '[ N ]'
                            app.Max_Load = app.MaxLoad.Value/app.g*1000;
                        case '[ kN ]'
                            app.Max_Load = app.MaxLoad.Value/app.g*1000000;
                    end
                    app.Gauge.Value = min(100, abs((grams/app.Max_Load))*100);
                end

                app.rawForce(app.i) = raw;
                app.force(app.i) = forceSelected;
                app.forceN(app.i) = forceNValue;
                app.pressure(app.i) = p;
                app.pressureV(app.i) = voltage;

                plot(app.Axes1, app.time(1:app.i), app.force(1:app.i));
                app.Pressure.Value = app.pressure(app.i);
                app.PressureGauge.Value = max(app.PressureGauge.Limits(1), min(app.PressureGauge.Limits(2), app.pressure(app.i)));
                app.ForceEdit.Value = app.force(app.i);
                app.TimeEdit.Value = app.time(app.i) - app.rst_time;
                app.DataEditField.Value = app.i - app.last_data;
                drawnow limitrate;

                pause(max(0, app.Add_time.Value - toc(loopTimer)));
                app.Rate.Value = toc(loopTimer);
                Xaxis(app);
                app.i = app.i + 1;
                app.time(app.i) = toc(app.start_time);
            end
            app.Cyan.Color = 'white';
            app.Yellow.Color = 'yellow';
            app.Message.Value = sprintf('Acquisition stopped. Samples in current run: %d.', app.i - 1 - app.last_data);
        end

        % Button pushed function: PauseButton
        function PauseButtonPushed(app, event)
            app.r = 1;
            app.Message.Value = 'In pause';
        end

        % Button pushed function: SaveButton
        function SaveButtonPushed(app, event)
            if ~app.get_true
                app.Message.Value = 'Error: You have not got any data yet.';
                return;
            end

            idx = (app.last_data + 1):(app.i - 1);
            if isempty(idx)
                app.Message.Value = 'Error: no unsaved samples in current run.';
                return;
            end

            filePath = buildSavePath(app);
            if isfile(filePath)
                answer = uiconfirm(app.MatlabArduinoHX711UIFigure, ...
                    sprintf('%s already exists. Overwrite?', filePath), ...
                    'Overwrite existing MAT-file?', 'Options', {'Overwrite','Cancel'}, ...
                    'DefaultOption', 'Cancel', 'CancelOption', 'Cancel');
                if ~strcmp(answer, 'Overwrite')
                    app.Message.Value = 'Save cancelled.';
                    return;
                end
            end

            % A = 6 columns below: time, raw HX711, force in N, pressure kPa,
            % pressure voltage, force in user-selected display unit.
            ColumnNames = {'Time_s','RawHX711_counts','Force_N','Pressure_kPa','PressureVoltage_V','Force_SelectedUnit'};
            A = numel(ColumnNames);
            Data = NaN(750, A);
            nRows = min(750, numel(idx));
            idx = idx(1:nRows);
            Data(1:nRows,:) = [app.time(idx).', app.rawForce(idx).', app.forceN(idx).', ...
                app.pressure(idx).', app.pressureV(idx).', app.force(idx).'];

            Stats = buildStatsTable(app, Data(1:nRows,3), Data(1:nRows,4));

            Metadata = struct();
            Metadata.Created = datetime('now');
            Metadata.Series = round(app.SeriesSpinner.Value);
            Metadata.Run = round(app.RunSpinner.Value);
            Metadata.FileName = char(filePath);
            Metadata.NsamplesSaved = nRows;
            Metadata.DataRows = 750;
            Metadata.DataColumns = A;
            Metadata.SelectedForceUnit = selectedUnitText(app);
            Metadata.KneeAngle_deg = app.KneeAngle.Value;
            Metadata.LoadCellAngle_deg = app.LoadCellAngle.Value;
            Metadata.LoadCellTare = app.tare;
            Metadata.LoadCellScale = app.scale;
            Metadata.PressureCalibration = struct('a_kPaPerV', app.pressureA, 'b_kPa', app.pressureB, 'Equation', 'Pressure_kPa = a*Voltage_V + b');
            Metadata.Pins = struct('Clock', app.clock, 'Data', app.data, 'Pressure', app.pressurePin, 'ValveD11Function', app.valveIncPin, 'ValveD6Function', app.valveMaintainPin);
            Metadata.ValveLogic = struct('Increasing','D11 High, D6 High', 'Maintain','D11 Low, D6 High', 'DecreasingOrZero','D11 Low, D6 Low');

            save(filePath, 'Data', 'Stats', 'Metadata', 'ColumnNames');
            app.Name.Value = char(filePath);
            app.last_data = app.i - 1;
            app.rst_time = app.time(app.i);
            app.Message.Value = ['Saved: ', char(filePath)];
        end

        % Button pushed function: TareButton
        function TareButtonPushed(app, event)
            if ~app.check_connection
                app.Message.Value = 'Error: You are not connected yet.';
                return;
            end
            if app.n.Value == 0
                app.Message.Value = 'Error: set number of readings.';
                return;
            end
            app.Message.Value = 'Please wait...';
            app.Cyan.Color = 'cyan';
            app.Yellow.Color = 'white';
            drawnow;
            z = round(app.n.Value);
            x = zeros(1,z);
            for j = 1:z
                x(j) = read_HX711(app.HX711_obj);
                pause(1/1000);
            end
            app.tare = mean(x);
            app.TareDisp.Value = app.tare;
            app.KnownTare.Value = app.tare;
            app.Message.Value = 'Tare phase is completed.';
            app.Cyan.Color = 'white';
            app.Yellow.Color = 'yellow';
        end

        % Button pushed function: ScaleFactorButton
        function ScaleFactorButtonPushed(app, event)
            if ~app.check_connection
                app.Message.Value = 'Error: You are not connected yet.';
                return;
            end
            if ~isfinite(app.tare)
                app.Message.Value = 'Error: perform tare phase first or enter known tare.';
                return;
            end
            switch app.UnitCal.Value
                case '[ g ]'
                    app.known_weight = app.Known.Value;
                case '[ kg ]'
                    app.known_weight = app.Known.Value*1000;
                case '[ N ]'
                    app.known_weight = app.Known.Value/app.g*1000;
                case '[ kN ]'
                    app.known_weight = app.Known.Value/app.g*1000000;
            end
            if app.known_weight == 0
                app.Message.Value = 'Error: check the known weight for calibration.';
                return;
            end
            app.Message.Value = 'Please wait...';
            app.Cyan.Color = 'cyan';
            app.Yellow.Color = 'white';
            drawnow;
            z = round(app.n.Value);
            x = zeros(1,z);
            for j = 1:z
                x(j) = read_HX711(app.HX711_obj);
                pause(1/1000);
            end
            app.scale = (mean(x) - app.tare)/app.known_weight;
            app.ScaleDisp.Value = app.scale;
            app.KnownScale.Value = app.scale;
            app.Message.Value = 'Scale factor is determined.';
            app.Cyan.Color = 'white';
            app.Yellow.Color = 'yellow';
        end

        % Button pushed function: ApplyKnownLoadCellButton
        function ApplyKnownLoadCellButtonPushed(app, event)
            if app.KnownScale.Value == 0
                app.Message.Value = 'Error: load-cell scale cannot be zero.';
                return;
            end
            app.tare = app.KnownTare.Value;
            app.scale = app.KnownScale.Value;
            app.TareDisp.Value = app.tare;
            app.ScaleDisp.Value = app.scale;
            app.Message.Value = 'Known load-cell tare and scale applied.';
        end

        % Button pushed function: PressureCalButton
        function PressureCalButtonPushed(app, event)
            if ~app.check_connection
                app.Message.Value = 'Error: You are not connected yet.';
                return;
            end
            setpoints = [0 200 300 400 500 600 620];
            actual = zeros(size(setpoints));
            voltage = zeros(size(setpoints));
            z = max(1, round(app.PressureCalN.Value));
            for j = 1:numel(setpoints)
                prompt = sprintf('Set regulator near %d kPa. Enter actual gauge kPa after pressure stabilizes:', setpoints(j));
                answer = inputdlg(prompt, 'Pressure calibration', [1 70], {num2str(setpoints(j))});
                if isempty(answer)
                    app.Message.Value = 'Pressure calibration cancelled.';
                    return;
                end
                actual(j) = str2double(answer{1});
                if ~isfinite(actual(j))
                    app.Message.Value = 'Pressure calibration error: actual kPa must be numeric.';
                    return;
                end
                readings = zeros(1,z);
                app.Message.Value = sprintf('Reading A0 at actual %.2f kPa...', actual(j));
                drawnow;
                for k = 1:z
                    readings(k) = readVoltage(app.a, app.pressurePin);
                    pause(0.02);
                end
                voltage(j) = mean(readings);
            end
            coeff = polyfit(voltage, actual, 1);
            app.pressureA = coeff(1);
            app.pressureB = coeff(2);
            app.PressureA.Value = app.pressureA;
            app.PressureB.Value = app.pressureB;
            plot(app.Axes2, voltage, actual, 'o', voltage, polyval(coeff, voltage), '-');
            xlabel(app.Axes2, 'Voltage [V]');
            ylabel(app.Axes2, 'Pressure [kPa]');
            app.Message.Value = sprintf('Pressure calibration complete: y = %.6g*x %+.6g', app.pressureA, app.pressureB);
        end

        % Button pushed function: ApplyKnownPressureButton
        function ApplyKnownPressureButtonPushed(app, event)
            app.pressureA = app.PressureA.Value;
            app.pressureB = app.PressureB.Value;
            app.Message.Value = sprintf('Known pressure calibration applied: y = %.6g*x %+.6g', app.pressureA, app.pressureB);
        end

        % Button pushed function: CalibrationButton
        function CalibrationButtonPushed(app, event)
            if ~app.check_connection
                app.Message.Value = 'Error: You are not connected yet.';
                return;
            end
            if ~isLoadCellCalibrated(app)
                app.Message.Value = 'Error: enter/apply load-cell tare and scale, or calibrate first.';
                return;
            end
            app.Message.Value = 'Please wait...';
            app.Cyan.Color = 'cyan';
            app.Yellow.Color = 'white';
            drawnow;
            switch app.UnitCal.Value
                case '[ g ]'
                    app.known_weight = app.Known.Value;
                case '[ kg ]'
                    app.known_weight = app.Known.Value*1000;
                case '[ N ]'
                    app.known_weight = app.Known.Value/app.g*1000;
                case '[ kN ]'
                    app.known_weight = app.Known.Value/app.g*1000000;
            end
            z = round(app.n.Value);
            x = zeros(1,z);
            for j = 1:z
                x(j) = read_HX711(app.HX711_obj);
                pause(1/1000);
            end
            x = rawToGrams(app, x);
            M = mean(x);
            S = std(x);
            app.AverageDisp.Value = M;
            app.StdDisp.Value = S;
            w = linspace(M-4*S, M+4*S, 5000);
            if S > 0
                y = (S*sqrt(2*pi))^-1*exp(-((w-M).^2)./(2*S^2));
                app.Axes2.LineWidth = 2;
                plot(app.Axes2, w, y, [app.known_weight app.known_weight], [0 max(y)*1.1], ...
                    [M-2*S M-2*S], [0 max(y)]/2, [M+2*S M+2*S], [0 max(y)]/2);
            else
                plot(app.Axes2, x, zeros(size(x)), 'o');
            end
            app.Message.Value = 'Calibration phase is completed.';
            app.Cyan.Color = 'white';
            app.Yellow.Color = 'yellow';
        end

        % Button pushed function: RawRead
        function RawReadButtonPushed(app, event)
            if ~app.check_connection
                app.Message.Value = 'Error: You are not connected yet.';
                return;
            end
            if ~isLoadCellCalibrated(app)
                app.Message.Value = 'Error: enter/apply load-cell tare and scale, or calibrate first.';
                return;
            end
            app.measure4.Value = selectedUnitText(app);
            app.Message.Value = 'Data acquisition';
            app.Cyan.Color = 'cyan';
            app.Yellow.Color = 'white';
            raw = read_HX711(app.HX711_obj);
            grams = rawToGrams(app, raw);
            x = gramsToSelectedUnit(app, grams);
            app.Raw.Value = x;
            app.Cyan.Color = 'white';
            app.Yellow.Color = 'yellow';
            app.Message.Value = 'Done';
        end

        % Button pushed function: Clean
        function CleanButtonPushed(app, event)
            cla(app.Axes1);
            cla(app.Axes2);
            app.time = [];
            app.force = [];
            app.forceN = [];
            app.rawForce = [];
            app.pressure = [];
            app.pressureV = [];
            app.i = 1;
            app.v = 1;
            app.t = 1;
            Xaxis2(app);
            app.last_time1 = 0;
            app.last_time2 = 0;
            app.rst_time = 0;
            app.k1 = app.i;
            app.k2 = app.i;
            app.last_data = 0;
            app.DataEditField.Value = 0;
            app.PressureGauge.Value = 0;
            app.ForceEdit.Value = 0;
            app.TimeEdit.Value = 0;
            app.AverageDisp.Value = 0;
            app.StdDisp.Value = 0;
            app.Message.Value = 'All cleaned up.';
        end

        % Valve-control callbacks
        function IncreasePressureButtonPushed(app, event)
            setValves(app, 1, 1, 'Valves: pressure increasing (D11 High, D6 High).');
        end

        function MaintainPressureButtonPushed(app, event)
            setValves(app, 0, 1, 'Valves: pressure maintained (D11 Low, D6 High).');
        end

        function DecreasePressureButtonPushed(app, event)
            setValves(app, 0, 0, 'Valves: pressure decreasing / 0 kPa maintain (D11 Low, D6 Low).');
        end
    end

    % Component initialization
    methods (Access = private)

        function createComponents(app)
            app.MatlabArduinoHX711UIFigure = uifigure('Visible', 'off');
            app.MatlabArduinoHX711UIFigure.Color = [0.9412 0.9412 0.9412];
            app.MatlabArduinoHX711UIFigure.Position = [100 100 1120 900];
            app.MatlabArduinoHX711UIFigure.Name = 'Matlab - Arduino - HX711';
            app.MatlabArduinoHX711UIFigure.Resize = 'off';

            app.MessageEditFieldLabel = uilabel(app.MatlabArduinoHX711UIFigure);
            app.MessageEditFieldLabel.HorizontalAlignment = 'right';
            app.MessageEditFieldLabel.FontName = 'Verdana';
            app.MessageEditFieldLabel.FontAngle = 'italic';
            app.MessageEditFieldLabel.Position = [15 286 55 15];
            app.MessageEditFieldLabel.Text = 'Message';

            app.Message = uieditfield(app.MatlabArduinoHX711UIFigure, 'text');
            app.Message.Editable = 'off';
            app.Message.FontName = 'Verdana';
            app.Message.FontWeight = 'bold';
            app.Message.FontAngle = 'italic';
            app.Message.FontColor = [1 0 0];
            app.Message.Position = [13 241 350 38];

            app.StatusPanel = uipanel(app.MatlabArduinoHX711UIFigure);
            app.StatusPanel.TitlePosition = 'centertop';
            app.StatusPanel.Title = 'Status';
            app.StatusPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.StatusPanel.FontName = 'Verdana';
            app.StatusPanel.FontAngle = 'italic';
            app.StatusPanel.FontWeight = 'bold';
            app.StatusPanel.Position = [13 685 185 189];

            app.NotConnectedLabel = uilabel(app.StatusPanel);
            app.NotConnectedLabel.HorizontalAlignment = 'center';
            app.NotConnectedLabel.FontName = 'Verdana';
            app.NotConnectedLabel.FontSize = 14;
            app.NotConnectedLabel.FontAngle = 'italic';
            app.NotConnectedLabel.Position = [16 140 168 19];
            app.NotConnectedLabel.Text = 'Not Connected';
            app.Red = uilamp(app.StatusPanel);
            app.Red.Position = [138 134 30 30];
            app.Red.Color = [1 0 0];

            app.ConnectedLabel = uilabel(app.StatusPanel);
            app.ConnectedLabel.HorizontalAlignment = 'center';
            app.ConnectedLabel.FontName = 'Verdana';
            app.ConnectedLabel.FontSize = 14;
            app.ConnectedLabel.FontAngle = 'italic';
            app.ConnectedLabel.Position = [41 100 114 19];
            app.ConnectedLabel.Text = 'Connected';
            app.Green = uilamp(app.StatusPanel);
            app.Green.Position = [138 94 30 30];
            app.Green.Color = [1 1 1];

            app.InpauseLabel = uilabel(app.StatusPanel);
            app.InpauseLabel.HorizontalAlignment = 'center';
            app.InpauseLabel.FontName = 'Verdana';
            app.InpauseLabel.FontSize = 14;
            app.InpauseLabel.FontAngle = 'italic';
            app.InpauseLabel.Position = [60 59 76 19];
            app.InpauseLabel.Text = 'In pause';
            app.Yellow = uilamp(app.StatusPanel);
            app.Yellow.Position = [138 53 30 30];
            app.Yellow.Color = [1 1 1];

            app.DataAcquisitionLabel = uilabel(app.StatusPanel);
            app.DataAcquisitionLabel.HorizontalAlignment = 'center';
            app.DataAcquisitionLabel.FontName = 'Verdana';
            app.DataAcquisitionLabel.FontSize = 14;
            app.DataAcquisitionLabel.FontAngle = 'italic';
            app.DataAcquisitionLabel.Position = [7 17 185 19];
            app.DataAcquisitionLabel.Text = 'Data Acquisition';
            app.Cyan = uilamp(app.StatusPanel);
            app.Cyan.Position = [139 11 30 30];
            app.Cyan.Color = [1 1 1];

            app.ArduinoHX711Panel = uipanel(app.MatlabArduinoHX711UIFigure);
            app.ArduinoHX711Panel.TitlePosition = 'centertop';
            app.ArduinoHX711Panel.Title = 'Arduino - HX711';
            app.ArduinoHX711Panel.FontName = 'Verdana';
            app.ArduinoHX711Panel.FontAngle = 'italic';
            app.ArduinoHX711Panel.FontWeight = 'bold';
            app.ArduinoHX711Panel.Position = [216 685 147 189];

            app.Connect = uibutton(app.ArduinoHX711Panel, 'push');
            app.Connect.ButtonPushedFcn = createCallbackFcn(app, @ConnectButtonPushed, true);
            app.Connect.FontName = 'Verdana';
            app.Connect.FontWeight = 'bold';
            app.Connect.FontAngle = 'italic';
            app.Connect.Position = [21 139 105 21];
            app.Connect.Text = 'Connect';

            app.GetData = uibutton(app.ArduinoHX711Panel, 'push');
            app.GetData.ButtonPushedFcn = createCallbackFcn(app, @GetDataButtonPushed, true);
            app.GetData.FontName = 'Verdana';
            app.GetData.FontWeight = 'bold';
            app.GetData.FontAngle = 'italic';
            app.GetData.Position = [21 98 105 22];
            app.GetData.Text = 'Get Data';

            app.PauseButton = uibutton(app.ArduinoHX711Panel, 'push');
            app.PauseButton.ButtonPushedFcn = createCallbackFcn(app, @PauseButtonPushed, true);
            app.PauseButton.FontName = 'Verdana';
            app.PauseButton.FontWeight = 'bold';
            app.PauseButton.FontAngle = 'italic';
            app.PauseButton.Position = [21 57 105 22];
            app.PauseButton.Text = 'Pause';

            app.SaveButton = uibutton(app.ArduinoHX711Panel, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.FontName = 'Verdana';
            app.SaveButton.FontWeight = 'bold';
            app.SaveButton.FontAngle = 'italic';
            app.SaveButton.Position = [21 15 105 22];
            app.SaveButton.Text = 'Save';

            app.CalibrationPanel = uipanel(app.MatlabArduinoHX711UIFigure);
            app.CalibrationPanel.TitlePosition = 'centertop';
            app.CalibrationPanel.Title = 'Calibration';
            app.CalibrationPanel.FontName = 'Verdana';
            app.CalibrationPanel.FontAngle = 'italic';
            app.CalibrationPanel.FontWeight = 'bold';
            app.CalibrationPanel.Position = [216 482 147 188];

            app.TareButton = uibutton(app.CalibrationPanel, 'push');
            app.TareButton.ButtonPushedFcn = createCallbackFcn(app, @TareButtonPushed, true);
            app.TareButton.FontName = 'Verdana';
            app.TareButton.FontWeight = 'bold';
            app.TareButton.FontAngle = 'italic';
            app.TareButton.Position = [21 137 105 22];
            app.TareButton.Text = 'Tare';

            app.ScaleFactorButton = uibutton(app.CalibrationPanel, 'push');
            app.ScaleFactorButton.ButtonPushedFcn = createCallbackFcn(app, @ScaleFactorButtonPushed, true);
            app.ScaleFactorButton.FontName = 'Verdana';
            app.ScaleFactorButton.FontWeight = 'bold';
            app.ScaleFactorButton.FontAngle = 'italic';
            app.ScaleFactorButton.Position = [21 97 105 22];
            app.ScaleFactorButton.Text = 'Scale Factor';

            app.CalibrationButton = uibutton(app.CalibrationPanel, 'push');
            app.CalibrationButton.ButtonPushedFcn = createCallbackFcn(app, @CalibrationButtonPushed, true);
            app.CalibrationButton.FontName = 'Verdana';
            app.CalibrationButton.FontWeight = 'bold';
            app.CalibrationButton.FontAngle = 'italic';
            app.CalibrationButton.Position = [21 56 105 22];
            app.CalibrationButton.Text = 'Calibration';

            app.RawRead = uibutton(app.CalibrationPanel, 'push');
            app.RawRead.ButtonPushedFcn = createCallbackFcn(app, @RawReadButtonPushed, true);
            app.RawRead.FontName = 'Verdana';
            app.RawRead.FontWeight = 'bold';
            app.RawRead.FontAngle = 'italic';
            app.RawRead.Position = [21 15 105 22];
            app.RawRead.Text = 'Get Data';

            app.GlobalSettingsPanel = uipanel(app.MatlabArduinoHX711UIFigure);
            app.GlobalSettingsPanel.TitlePosition = 'centertop';
            app.GlobalSettingsPanel.Title = 'Global Settings';
            app.GlobalSettingsPanel.FontName = 'Verdana';
            app.GlobalSettingsPanel.FontAngle = 'italic';
            app.GlobalSettingsPanel.FontWeight = 'bold';
            app.GlobalSettingsPanel.Position = [13 482 185 188];

            app.TabGroup = uitabgroup(app.GlobalSettingsPanel);
            app.TabGroup.Position = [1 -25 185 193];

            app.ConnectionTab = uitab(app.TabGroup);
            app.ConnectionTab.Title = 'Connection';
            app.ArduinoDropDownLabel = uilabel(app.ConnectionTab, 'Text', 'Arduino', 'Position', [5 139 49 15]);
            app.BoardEdit = uidropdown(app.ConnectionTab, 'Items', {'Uno','Mega2560'}, 'Position', [77 135 104 22], 'Value', 'Uno');
            app.SerialportEditFieldLabel = uilabel(app.ConnectionTab, 'Text', 'Serial port', 'Position', [0 106 66 15]);
            app.SerialEdit = uieditfield(app.ConnectionTab, 'text', 'Position', [112 102 69 22], 'Value', 'Com4');
            app.DataPinEditFieldLabel = uilabel(app.ConnectionTab, 'Text', 'Data Pin', 'Position', [0 74 55 15]);
            app.DataEdit = uieditfield(app.ConnectionTab, 'text', 'Position', [112 70 69 22], 'Value', 'D3');
            app.ClockPinEditFieldLabel = uilabel(app.ConnectionTab, 'Text', 'Clock Pin', 'Position', [1 42 58 15]);
            app.ClockEdit = uieditfield(app.ConnectionTab, 'text', 'Position', [112 38 69 22], 'Value', 'D2');
            app.PressurePinLabel = uilabel(app.ConnectionTab, 'Text', 'Pressure Pin', 'Position', [4 10 74 22]);
            app.PressureEdit = uieditfield(app.ConnectionTab, 'text', 'Position', [111 12 69 22], 'Value', 'A0');

            app.DataAcquisitionTab = uitab(app.TabGroup);
            app.DataAcquisitionTab.Title = 'Data Acquisition';
            app.ForceDropDownLabel = uilabel(app.DataAcquisitionTab, 'Text', 'Force', 'Position', [7 139 37 15]);
            app.Unit = uidropdown(app.DataAcquisitionTab, 'Items', {'[ N ]','[ g ]','[ kg ]','[ kN ]'}, 'Position', [104 135 77 22], 'Value', '[ N ]');
            app.SetSession = uicheckbox(app.DataAcquisitionTab, 'Text', 'Set Session Time [min]', 'Position', [9 109 161 15]);
            app.SessionTimeSpinnerLabel = uilabel(app.DataAcquisitionTab, 'Text', 'Session Time', 'Position', [2 73 83 15]);
            app.SessionV = uispinner(app.DataAcquisitionTab, 'Limits', [0 Inf], 'Position', [103 69 77 22]);
            app.SamplingRateSpinnerLabel = uilabel(app.DataAcquisitionTab, 'Text', 'Sample Period [s]', 'Position', [2 44 98 15]);
            app.Add_time = uispinner(app.DataAcquisitionTab, 'Step', 0.1, 'Limits', [0.001 Inf], 'ValueDisplayFormat', '%.3f', 'Position', [104 40 76 22], 'Value', 0.1);
            app.SampleCountSpinnerLabel = uilabel(app.DataAcquisitionTab, 'Text', 'Samples', 'Position', [2 14 60 15]);
            app.Nsamples = uispinner(app.DataAcquisitionTab, 'Limits', [1 750], 'ValueDisplayFormat', '%.0f', 'Position', [104 10 76 22], 'Value', 750);

            app.SaveDataTab = uitab(app.TabGroup);
            app.SaveDataTab.Title = 'Save Data';
            app.SaveFolderLabel = uilabel(app.SaveDataTab, 'Text', 'Folder', 'Position', [7 139 40 15]);
            app.SaveFolder = uieditfield(app.SaveDataTab, 'text', 'Position', [52 135 129 22], 'Value', fullfile('GitHub','Bipedal_Robot','Testing_Data','2026_06_Festo','Flx_20mm'));
            app.SeriesSpinnerLabel = uilabel(app.SaveDataTab, 'Text', 'Series $', 'Position', [7 104 52 15]);
            app.SeriesSpinner = uispinner(app.SaveDataTab, 'Limits', [0 999], 'ValueDisplayFormat', '%.0f', 'Position', [72 100 45 22], 'Value', 1);
            app.RunSpinnerLabel = uilabel(app.SaveDataTab, 'Text', 'Run #', 'Position', [7 73 42 15]);
            app.RunSpinner = uispinner(app.SaveDataTab, 'Limits', [0 99], 'ValueDisplayFormat', '%02.0f', 'Position', [72 69 45 22], 'Value', 0);
            app.NameEditFieldLabel = uilabel(app.SaveDataTab, 'Text', 'Last file', 'Position', [7 42 50 15]);
            app.Name = uieditfield(app.SaveDataTab, 'text', 'Editable', 'off', 'Position', [60 38 120 22]);
            app.Reset = uicheckbox(app.SaveDataTab, 'Text', 'Reset Time after Save', 'Position', [12 13 160 15], 'Value', true);

            app.MetadataTab = uitab(app.TabGroup);
            app.MetadataTab.Title = 'Metadata';
            app.KneeAngleLabel = uilabel(app.MetadataTab, 'Text', 'Knee angle [deg]', 'Position', [7 119 115 15]);
            app.KneeAngle = uieditfield(app.MetadataTab, 'numeric', 'Position', [118 115 62 22], 'ValueDisplayFormat', '%.2f');
            app.LoadCellAngleLabel = uilabel(app.MetadataTab, 'Text', 'Load cell angle [deg]', 'Position', [7 82 130 15]);
            app.LoadCellAngle = uieditfield(app.MetadataTab, 'numeric', 'Position', [118 78 62 22], 'ValueDisplayFormat', '%.2f');
            app.ValveIncPinLabel = uilabel(app.MetadataTab, 'Text', 'Valve D11 pin', 'Position', [7 45 84 15]);
            app.ValveIncEdit = uieditfield(app.MetadataTab, 'text', 'Position', [118 41 62 22], 'Value', 'D11');
            app.ValveMaintainPinLabel = uilabel(app.MetadataTab, 'Text', 'Valve D6 pin', 'Position', [7 14 84 15]);
            app.ValveMaintainEdit = uieditfield(app.MetadataTab, 'text', 'Position', [118 10 62 22], 'Value', 'D6');

            app.ContinuosDataAcquisitionPanel = uipanel(app.MatlabArduinoHX711UIFigure);
            app.ContinuosDataAcquisitionPanel.Title = 'Continuous Data Acquisition';
            app.ContinuosDataAcquisitionPanel.FontName = 'Verdana';
            app.ContinuosDataAcquisitionPanel.FontAngle = 'italic';
            app.ContinuosDataAcquisitionPanel.Position = [384 591 720 282];
            app.Axes1 = uiaxes(app.ContinuosDataAcquisitionPanel);
            xlabel(app.Axes1, 't [s]');
            ylabel(app.Axes1, 'Force [N]');
            app.Axes1.XLim = [0 100];
            app.Axes1.XGrid = 'on';
            app.Axes1.YGrid = 'on';
            app.Axes1.Box = 'on';
            app.Axes1.Position = [12 9 695 249];

            app.CalibrationResultPanel = uipanel(app.MatlabArduinoHX711UIFigure);
            app.CalibrationResultPanel.Title = 'Calibration Result';
            app.CalibrationResultPanel.FontName = 'Verdana';
            app.CalibrationResultPanel.FontAngle = 'italic';
            app.CalibrationResultPanel.Position = [384 282 720 300];
            app.Axes2 = uiaxes(app.CalibrationResultPanel);
            app.Axes2.XGrid = 'on'; app.Axes2.YGrid = 'on'; app.Axes2.Box = 'on';
            app.Axes2.Position = [8 8 350 265];
            app.TabGroup2 = uitabgroup(app.CalibrationResultPanel);
            app.TabGroup2.Position = [367 8 344 257];

            app.ValueTab = uitab(app.TabGroup2); app.ValueTab.Title = 'Value';
            app.Gauge = uigauge(app.ValueTab, 'linear', 'Orientation', 'vertical', 'Position', [285 10 51 218]);
            app.TareLabel = uilabel(app.ValueTab, 'Text', 'Tare', 'Position', [11 207 35 15]);
            app.TareDisp = uieditfield(app.ValueTab, 'numeric', 'ValueDisplayFormat', '%.0f', 'Editable', 'off', 'Position', [102 203 110 22]);
            app.ScalefactorLabel = uilabel(app.ValueTab, 'Text', 'Scale factor', 'Position', [11 171 84 15]);
            app.ScaleDisp = uieditfield(app.ValueTab, 'numeric', 'ValueDisplayFormat', '%.6g', 'Editable', 'off', 'Position', [102 167 110 22]);
            app.AveragegLabel = uilabel(app.ValueTab, 'Text', 'Average [g]', 'Position', [11 135 86 15]);
            app.AverageDisp = uieditfield(app.ValueTab, 'numeric', 'ValueDisplayFormat', '%.2f', 'Editable', 'off', 'Position', [102 131 110 22]);
            app.StdDeviationgLabel = uilabel(app.ValueTab, 'Text', 'Std Deviation [g]', 'Position', [11 101 121 15]);
            app.StdDisp = uieditfield(app.ValueTab, 'numeric', 'ValueDisplayFormat', '%.2f', 'Editable', 'off', 'Position', [131 97 81 22]);
            app.RawreadingLabel = uilabel(app.ValueTab, 'Text', 'Raw reading', 'Position', [11 70 89 15]);
            app.Raw = uieditfield(app.ValueTab, 'numeric', 'ValueDisplayFormat', '%.2f', 'Editable', 'off', 'Position', [122 66 90 22]);
            app.measure4 = uieditfield(app.ValueTab, 'text', 'Editable', 'off', 'Position', [176 33 36 22]);

            app.SettingsTab = uitab(app.TabGroup2); app.SettingsTab.Title = 'Settings';
            app.NumreadingsLabel = uilabel(app.SettingsTab, 'Text', 'Num. readings', 'Position', [11 207 103 15]);
            app.n = uieditfield(app.SettingsTab, 'numeric', 'ValueDisplayFormat', '%.0f', 'Position', [167 203 75 22], 'Value', 20);
            app.KnownWeightLabel = uilabel(app.SettingsTab, 'Text', 'Known Weight', 'Position', [10 171 104 15]);
            app.Known = uieditfield(app.SettingsTab, 'numeric', 'Limits', [0 Inf], 'ValueDisplayFormat', '%.3f', 'Position', [134 167 107 22]);
            app.UnitLabel = uilabel(app.SettingsTab, 'Text', 'Unit', 'Position', [11 139 33 15]);
            app.UnitCal = uidropdown(app.SettingsTab, 'Items', {'[ g ]','[ kg ]','[ N ]','[ kN ]'}, 'Position', [164 135 77 22], 'Value', '[ g ]');
            app.MaxLoadLabel = uilabel(app.SettingsTab, 'Text', 'Max Load', 'Position', [10 104 68 15]);
            app.MaxLoad = uieditfield(app.SettingsTab, 'numeric', 'Limits', [0 Inf], 'ValueDisplayFormat', '%.2f', 'Position', [134 100 107 22]);
            app.UnitLabel_2 = uilabel(app.SettingsTab, 'Text', 'Unit', 'Position', [11 70 33 15]);
            app.UnitCal2 = uidropdown(app.SettingsTab, 'Items', {'[ g ]','[ kg ]','[ N ]','[ kN ]'}, 'Position', [164 66 77 22], 'Value', '[ g ]');

            app.KnownCalTab = uitab(app.TabGroup2); app.KnownCalTab.Title = 'Known LC Cal';
            app.KnownTareLabel = uilabel(app.KnownCalTab, 'Text', 'Zero offset / tare', 'Position', [10 190 120 15]);
            app.KnownTare = uieditfield(app.KnownCalTab, 'numeric', 'ValueDisplayFormat', '%.6g', 'Position', [150 186 170 22]);
            app.KnownScaleLabel = uilabel(app.KnownCalTab, 'Text', 'Calibration slope', 'Position', [10 150 120 15]);
            app.KnownScale = uieditfield(app.KnownCalTab, 'numeric', 'ValueDisplayFormat', '%.12g', 'Position', [150 146 170 22], 'Value', 1);
            app.ApplyKnownLoadCellButton = uibutton(app.KnownCalTab, 'push', 'Text', 'Apply Known Load-Cell Cal', 'Position', [70 92 200 28]);
            app.ApplyKnownLoadCellButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyKnownLoadCellButtonPushed, true);

            app.PressureCalTab = uitab(app.TabGroup2); app.PressureCalTab.Title = 'Pressure Cal';
            app.PressureALabel = uilabel(app.PressureCalTab, 'Text', 'a [kPa/V]', 'Position', [10 205 80 15]);
            app.PressureA = uieditfield(app.PressureCalTab, 'numeric', 'ValueDisplayFormat', '%.12g', 'Position', [120 201 130 22], 'Value', 155.61);
            app.PressureBLabel = uilabel(app.PressureCalTab, 'Text', 'b [kPa]', 'Position', [10 172 80 15]);
            app.PressureB = uieditfield(app.PressureCalTab, 'numeric', 'ValueDisplayFormat', '%.12g', 'Position', [120 168 130 22], 'Value', -126.99);
            app.PressureCalNLabel = uilabel(app.PressureCalTab, 'Text', 'Readings/point', 'Position', [10 139 100 15]);
            app.PressureCalN = uieditfield(app.PressureCalTab, 'numeric', 'Limits', [1 Inf], 'ValueDisplayFormat', '%.0f', 'Position', [120 135 130 22], 'Value', 25);
            app.ApplyKnownPressureButton = uibutton(app.PressureCalTab, 'push', 'Text', 'Apply Known Pressure Cal', 'Position', [65 86 210 28]);
            app.ApplyKnownPressureButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyKnownPressureButtonPushed, true);
            app.PressureCalButton = uibutton(app.PressureCalTab, 'push', 'Text', 'Run 7-Point Pressure Cal', 'Position', [65 45 210 28]);
            app.PressureCalButton.ButtonPushedFcn = createCallbackFcn(app, @PressureCalButtonPushed, true);

            app.LicenseTab = uitab(app.TabGroup2); app.LicenseTab.Title = '*License*';
            app.TextArea = uitextarea(app.LicenseTab, 'Editable', 'off', 'Position', [10 9 320 215]);
            app.TextArea.Value = {'Original HX711 app copyright 2018, Nicholas Giacoboni'; 'Modified for Bipedal_Robot data collection.'};

            app.CleanPanel = uipanel(app.MatlabArduinoHX711UIFigure);
            app.CleanPanel.TitlePosition = 'centertop';
            app.CleanPanel.Title = 'Clean';
            app.CleanPanel.FontName = 'Verdana';
            app.CleanPanel.FontAngle = 'italic';
            app.CleanPanel.FontWeight = 'bold';
            app.CleanPanel.Position = [216 359 147 111];
            app.Clean = uibutton(app.CleanPanel, 'push');
            app.Clean.ButtonPushedFcn = createCallbackFcn(app, @CleanButtonPushed, true);
            app.Clean.FontName = 'Verdana';
            app.Clean.FontWeight = 'bold';
            app.Clean.FontAngle = 'italic';
            app.Clean.Position = [21 41 105 22];
            app.Clean.Text = 'Plot & Data';

            app.ValvePanel = uipanel(app.MatlabArduinoHX711UIFigure);
            app.ValvePanel.TitlePosition = 'centertop';
            app.ValvePanel.Title = 'Valves';
            app.ValvePanel.FontName = 'Verdana';
            app.ValvePanel.FontAngle = 'italic';
            app.ValvePanel.FontWeight = 'bold';
            app.ValvePanel.Position = [216 241 147 105];
            app.IncreasePressureButton = uibutton(app.ValvePanel, 'push', 'Text', 'Increase', 'Position', [21 58 105 22]);
            app.IncreasePressureButton.ButtonPushedFcn = createCallbackFcn(app, @IncreasePressureButtonPushed, true);
            app.MaintainPressureButton = uibutton(app.ValvePanel, 'push', 'Text', 'Maintain', 'Position', [21 33 105 22]);
            app.MaintainPressureButton.ButtonPushedFcn = createCallbackFcn(app, @MaintainPressureButtonPushed, true);
            app.DecreasePressureButton = uibutton(app.ValvePanel, 'push', 'Text', 'Decrease / 0', 'Position', [21 8 105 22]);
            app.DecreasePressureButton.ButtonPushedFcn = createCallbackFcn(app, @DecreasePressureButtonPushed, true);

            app.ForceEditFieldLabel = uilabel(app.MatlabArduinoHX711UIFigure, 'Text', 'Force', 'Position', [14 426 42 15]);
            app.ForceEdit = uieditfield(app.MatlabArduinoHX711UIFigure, 'numeric', 'ValueDisplayFormat', '%.2f', 'Editable', 'off', 'Position', [67 422 82 22]);
            app.measure = uieditfield(app.MatlabArduinoHX711UIFigure, 'text', 'Editable', 'off', 'Position', [159 421 36 22], 'Value', 'N');
            app.PressureLabel = uilabel(app.MatlabArduinoHX711UIFigure, 'Text', 'Pressure', 'Position', [15 391 57 22]);
            app.Pressure = uieditfield(app.MatlabArduinoHX711UIFigure, 'numeric', 'ValueDisplayFormat', '%.2f', 'Editable', 'off', 'Position', [68 394 82 22]);
            app.kPaLabel = uilabel(app.MatlabArduinoHX711UIFigure, 'Text', 'kPa', 'Position', [161 387 30 22]);
            app.TimeEditFieldLabel = uilabel(app.MatlabArduinoHX711UIFigure, 'Text', 'Time', 'Position', [13 373 38 15]);
            app.TimeEdit = uieditfield(app.MatlabArduinoHX711UIFigure, 'numeric', 'ValueDisplayFormat', '%.1f', 'Editable', 'off', 'Position', [67 369 82 22]);
            app.measure2 = uieditfield(app.MatlabArduinoHX711UIFigure, 'text', 'Editable', 'off', 'Position', [158 366 36 22], 'Value', 's');
            app.SamplingRLabel = uilabel(app.MatlabArduinoHX711UIFigure, 'Text', 'Sampling R.', 'Position', [13 343 85 15]);
            app.Rate = uieditfield(app.MatlabArduinoHX711UIFigure, 'numeric', 'ValueDisplayFormat', '%.3f', 'Editable', 'off', 'Position', [105 339 44 22]);
            app.measure3 = uieditfield(app.MatlabArduinoHX711UIFigure, 'text', 'Editable', 'off', 'Position', [158 339 36 22], 'Value', 's');
            app.DataEditFieldLabel = uilabel(app.MatlabArduinoHX711UIFigure, 'Text', '# Data', 'Position', [18 317 50 15]);
            app.DataEditField = uieditfield(app.MatlabArduinoHX711UIFigure, 'numeric', 'ValueDisplayFormat', '%.0f', 'Editable', 'off', 'Position', [77 313 121 22]);

            app.PressureGaugeLabel = uilabel(app.MatlabArduinoHX711UIFigure, 'HorizontalAlignment', 'center', 'Position', [139 21 93 22]);
            app.PressureGaugeLabel.Text = 'Pressure Gauge';
            app.PressureGauge = uigauge(app.MatlabArduinoHX711UIFigure, 'semicircular');
            app.PressureGauge.Limits = [0 700];
            app.PressureGauge.Position = [61 58 249 135];

            app.MatlabArduinoHX711UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)
        function app = HX711
            createComponents(app)
            registerApp(app, app.MatlabArduinoHX711UIFigure)
            if nargout == 0
                clear app
            end
        end

        function delete(app)
            delete(app.MatlabArduinoHX711UIFigure)
        end
    end
end
