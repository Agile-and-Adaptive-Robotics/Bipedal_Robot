classdef HX711_BPA < matlab.apps.AppBase

    % HX711_BPA - BPA force + pressure test app (AARL, Bipedal_Robot).
    %
    % Rebuild of Ben's customized HX711 app (encoder-angle metadata +
    % pressure sensor). Two ways to calibrate, both always available:
    %   1. Original workflow: Tare -> Scale Factor (known weight) ->
    %      Calibration check (Gaussian overlay).
    %   2. Known-factor entry: type the zero offset/tare and calibration
    %      slope into the "Known LC Cal" tab, or the pressure-sensor a/b
    %      into the "Pressure Cal" tab, and click Apply - no hardware
    %      step needed. A guided 7-point pressure calibration is also
    %      provided on the same tab.
    %
    % Dynamic pressure calibration ("Pressure Ctrl" tab): a PID pressure
    % controller drives the valves by time-proportioning - PID duty
    % u in [-100,100]% becomes a FILL pulse (u>0, D11+D6 High), a VENT
    % pulse (u<0, both Low), or HOLD (D11 Low, D6 High) each control
    % period. "Run Step Test" steps a free BPA to the setpoint, plots the
    % response with the deadband, reports overshoot/undershoot (kPa and
    % % of step), 10-90% rise time, settling time, and steady-state
    % error, and saves DPC_S##_R##.mat traces for tuning comparisons.
    % The same PID (checkbox on) holds pressure during Get Data.
    %
    % The factors entered (or measured) are remembered across sessions in
    % hx711_bpa_last_cal.mat next to this file (machine-local, gitignored),
    % including the PID gains.
    %
    % Launch with Start_HX711_BPA (adds this folder to the path so the
    % Arduino add-on package +arduinoioaddons/+basicHX711 resolves on any
    % machine and any clone location - no install step, no absolute paths).
    %
    % Save convention (Save button):
    %   <save folder>\<Prefix><Series>_<Run>.mat  plus a matching .txt
    %   sidecar in the original app's 2-column format (Time, Force), so
    %   existing text-file tooling keeps working.
    %   Default save folder: <repo>\Testing_Data\2026_06_Festo\Flx_20mm,
    %   found RELATIVE to this file (walks up to the folder containing
    %   Testing_Data); falls back to <this folder>\saved_data.
    %
    % MAT-file variables:
    %   Data       - 750 x 6 numeric array (columns named in ColumnNames:
    %                Time_s, RawHX711_counts, Force_N, Pressure_kPa,
    %                PressureVoltage_V, Force_SelectedUnit).
    %   Stats      - table (Mean/Median/Mode/Min/Max/StdDev x Force/Pressure).
    %   Metadata   - struct: angles, pins, calibration factors, valve logic,
    %                pressure-servo settings, series/run, MATLAB version.
    %   ColumnNames- cell array naming the Data columns.

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
        EnablePressureControl          matlab.ui.control.CheckBox
        DesiredPressure                matlab.ui.control.Spinner
        SetpointkPaLabel               matlab.ui.control.Label
        PressureDeadband               matlab.ui.control.Spinner
        PressureCtrlTab                matlab.ui.container.Tab
        Kp                             matlab.ui.control.NumericEditField
        KpLabel                        matlab.ui.control.Label
        Ki                             matlab.ui.control.NumericEditField
        KiLabel                        matlab.ui.control.Label
        Kd                             matlab.ui.control.NumericEditField
        KdLabel                        matlab.ui.control.Label
        CtrlPeriod                     matlab.ui.control.Spinner
        CtrlPeriodLabel                matlab.ui.control.Label
        RunStepTestButton              matlab.ui.control.Button
        SaveDataTab                    matlab.ui.container.Tab
        Name                           matlab.ui.control.EditField
        NameEditFieldLabel             matlab.ui.control.Label
        SaveFolder                     matlab.ui.control.EditField
        SaveFolderLabel                matlab.ui.control.Label
        Prefix                         matlab.ui.control.EditField
        PrefixLabel                    matlab.ui.control.Label
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
        KnownTare                      matlab.ui.control.EditField
        KnownTareLabel                 matlab.ui.control.Label
        KnownScale                     matlab.ui.control.EditField
        KnownScaleLabel                matlab.ui.control.Label
        KnownCalHint                   matlab.ui.control.Label
        PressureCalTab                 matlab.ui.container.Tab
        PressureA                      matlab.ui.control.NumericEditField
        PressureALabel                 matlab.ui.control.Label
        PressureB                      matlab.ui.control.NumericEditField
        PressureBLabel                 matlab.ui.control.Label
        PressureCalN                   matlab.ui.control.NumericEditField
        PressureCalNLabel              matlab.ui.control.Label
        ApplyKnownPressureButton       matlab.ui.control.Button
        PressureCalButton              matlab.ui.control.Button
        PressureCalHint                matlab.ui.control.Label
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
        a = []            % Arduino object
        serial            % Serial port
        board             % Arduino board
        data              % Data pin
        clock             % Clock pin
        pressurePin       % Pressure pin
        valveIncPin       % Increase-pressure valve pin (default D11)
        valveMaintainPin  % Maintain-pressure valve pin (default D6)
        HX711_obj = []    % basic_HX711 add-on object
        tare = NaN        % zero offset in raw HX711 counts; NaN = not set
        scale = NaN       % raw counts per gram-equivalent; NaN = not set
        pressureA = 155.61 % kPa/V slope (default = old app's hard-coded line)
        pressureB = -126.99 % kPa intercept
        g = 9.80665
        check_connection = false
        v = 1
        t = 0
        session_time = Inf
        known_weight = 0
        Max_Load = 0
        last_data = 0
        rst_time = 0
        get_true = false
        start_time = []
        appRoot = ''      % folder containing this file (anchoring point)
        pidI = 0          % PID integrator state [valve duty %]
        pidPrevP = NaN    % previous pressure [kPa] for derivative-on-measurement
        dpcRun = 0        % dynamic pressure calibration run counter
    end

    properties (Access = public)
        i = 1           % Index of force/time/pressure arrays
        time = []       % Time [s]
        force = []      % Force in currently selected display units
        pressure = []   % Pressure [kPa]
        rawForce = []   % Raw HX711 counts
        forceN = []     % Force [N] (saving/stats)
        pressureV = []  % Pressure sensor voltage [V]
        r = 0
    end

    methods (Access = private)

        function Xaxis(app)
            if app.time(app.i) > 100*app.v || app.t == 1
                app.v = app.v + 1;
                app.Axes1.XLim = [0 100*app.v];
                app.Axes1.XTick = [0 10*app.v 20*app.v 30*app.v 40*app.v ...
                    50*app.v 60*app.v 70*app.v 80*app.v 90*app.v 100*app.v];
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
            % Either a real calibration or typed-in factors unlock acquisition.
            tf = isfinite(app.tare) && isfinite(app.scale) && app.scale > 0 ...
                && ~(app.tare == 0 && app.scale == 1);
        end

        function setWeightToGramsFromKnown(app)
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
        end

        function updateMaxLoadGaugeWithGrams(app, grams)
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
                app.Gauge.Value = min(100, abs(grams/app.Max_Load)*100);
            end
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

        function setValves(app, incState, maintainState, msg)
            if ~app.check_connection
                app.Message.Value = 'Error: You are not connected yet.';
                return;
            end
            writeDigitalPin(app.a, app.valveIncPin, logical(incState));
            writeDigitalPin(app.a, app.valveMaintainPin, logical(maintainState));
            app.Message.Value = msg;
        end

        function pidReset(app)
            app.pidI = 0;
            app.pidPrevP = NaN;
        end

        function u = pidCompute(app, p, dt, sp)
            % Discrete PID -> valve duty u in [-100, +100] %.
            % Positive u fills (FILL), negative u vents, ~0 holds.
            % Derivative acts on the MEASUREMENT (no setpoint kick);
            % integrator uses conditional integration (anti-windup).
            Kp = app.Kp.Value;
            Ki = app.Ki.Value;
            Kd = app.Kd.Value;
            e = sp - p;
            P = Kp*e;
            if isfinite(app.pidPrevP)
                D = -Kd*(p - app.pidPrevP)/max(dt, 1e-3);
            else
                D = 0;
            end
            Icand = app.pidI + Ki*e*dt;
            uUnsat = P + Icand + D;
            u = min(100, max(-100, uUnsat));
            saturated = (u ~= uUnsat);
            if ~saturated || (uUnsat > 100 && e < 0) || (uUnsat < -100 && e > 0)
                app.pidI = Icand;  % only integrate when it helps, not into the rail
            end
            app.pidPrevP = p;
        end

        function [st, tFill, tVent, tHold] = dutyPlan(app, u, Tc)
            % Map duty [-100,100] % onto one control period Tc:
            % +: FILL for the fraction, then HOLD; -: VENT, then HOLD.
            % |u| below MinDutyPct keeps the valves on HOLD (no chatter).
            minDuty = 2;  % % of the period; shortest useful valve pulse
            st = 0;
            tFill = 0;
            tVent = 0;
            tHold = Tc;
            if u >= minDuty
                st = 1;
                tFill = (u/100)*Tc;
                tHold = Tc - tFill;
            elseif u <= -minDuty
                st = -1;
                tVent = (-u/100)*Tc;
                tHold = Tc - tVent;
            end
        end

        function st = applyValveDuty(app, u, Tc)
            % Pulse the valves for one control period. Returns the state
            % applied: +1 FILL (both valves open), 0 HOLD, -1 VENT.
            [st, tFill, tVent, tHold] = dutyPlan(app, u, Tc);
            switch st
                case 1
                    writeDigitalPin(app.a, app.valveIncPin, 1);
                    writeDigitalPin(app.a, app.valveMaintainPin, 1);
                    pause(tFill);
                    writeDigitalPin(app.a, app.valveIncPin, 0);
                    pause(max(0, tHold));
                case -1
                    writeDigitalPin(app.a, app.valveIncPin, 0);
                    writeDigitalPin(app.a, app.valveMaintainPin, 0);
                    pause(tVent);
                    writeDigitalPin(app.a, app.valveMaintainPin, 1);
                    pause(max(0, tHold));
                otherwise
                    writeDigitalPin(app.a, app.valveIncPin, 0);
                    writeDigitalPin(app.a, app.valveMaintainPin, 1);
                    pause(max(0, tHold));
            end
        end

        function stats = buildStatsTable(~, forceN, pressureKPa)
            forceN = forceN(:);
            pressureKPa = pressureKPa(:);
            rowNames = {'Mean';'Median';'Mode';'Min';'Max';'StdDev'};
            forceCol = [mean(forceN,'omitnan'); median(forceN,'omitnan'); mode(forceN); ...
                min(forceN,[],'omitnan'); max(forceN,[],'omitnan'); std(forceN,0,'omitnan')];
            pressureCol = [mean(pressureKPa,'omitnan'); median(pressureKPa,'omitnan'); mode(pressureKPa); ...
                min(pressureKPa,[],'omitnan'); max(pressureKPa,[],'omitnan'); std(pressureKPa,0,'omitnan')];
            stats = table(forceCol, pressureCol, 'RowNames', rowNames, ...
                'VariableNames', {'Force','Pressure'});
        end

        function defaultDir = resolveDefaultSaveDir(app)
            % Find the repo root RELATIVE to this file (never an absolute
            % hard-coded path), so any clone on any machine works.
            defaultDir = fullfile(app.appRoot, 'saved_data');
            cand = app.appRoot;
            for k = 1:8
                if isfolder(fullfile(cand, 'Testing_Data'))
                    defaultDir = fullfile(cand, 'Testing_Data', '2026_06_Festo', 'Flx_20mm');
                    return;
                end
                [parent,~,~] = fileparts(cand);
                if strcmp(parent, cand) || isempty(parent)
                    return;
                end
                cand = parent;
            end
        end

        function filePath = buildSavePath(app)
            saveDir = strtrim(char(string(app.SaveFolder.Value)));
            if isempty(saveDir)
                saveDir = app.resolveDefaultSaveDir();
            end
            if ~isfolder(saveDir)
                mkdir(saveDir);
            end
            prefix = strtrim(app.Prefix.Value);
            if isempty(prefix)
                prefix = 'FlxTest';
            end
            fileName = sprintf('%s%d_%02d.mat', prefix, ...
                round(app.SeriesSpinner.Value), round(app.RunSpinner.Value));
            filePath = fullfile(saveDir, fileName);
        end

        function writeTxtSidecar(app, matFilePath, Data, nRows)
            % Same 2-column format the original app wrote (Time, Force in
            % the selected unit), so existing text-file tooling keeps working.
            [d, base] = fileparts(matFilePath);
            txtPath = fullfile(d, [base, '.txt']);
            fid = fopen(txtPath, 'wt');
            if fid == -1
                error('Could not open %s for writing.', txtPath);
            end
            fprintf(fid, '%11s %15s\n', 'Time [ s ]', ['Force [ ', selectedUnitText(app), ' ]']);
            for j = 1:nRows
                fprintf(fid, '%8.2f %15.2f\n', Data(j,1), Data(j,6));
            end
            fclose(fid);
        end

        function calFile = calCachePath(app)
            calFile = fullfile(app.appRoot, 'hx711_bpa_last_cal.mat');
        end

        function loadCalCache(app)
            % Prefill the last used calibration factors (machine-local file,
            % gitignored). Missing/corrupt file is silently ignored.
            calFile = calCachePath(app);
            if ~isfile(calFile)
                return;
            end
            try
                S = load(calFile, 'cal');
            catch
                return;
            end
            cal = S.cal;
            if isfield(cal, 'pressureA') && isfinite(cal.pressureA)
                app.pressureA = cal.pressureA;
                app.PressureA.Value = cal.pressureA;
            end
            if isfield(cal, 'pressureB') && isfinite(cal.pressureB)
                app.pressureB = cal.pressureB;
                app.PressureB.Value = cal.pressureB;
            end
            if isfield(cal, 'tare') && isfinite(cal.tare)
                app.tare = cal.tare;
                app.TareDisp.Value = cal.tare;
                app.KnownTare.Value = num2str(cal.tare, '%.6g');
            end
            if isfield(cal, 'scale') && isfinite(cal.scale) && cal.scale > 0
                app.scale = cal.scale;
                app.ScaleDisp.Value = cal.scale;
                app.KnownScale.Value = num2str(cal.scale, '%.12g');
            end
            if isfield(cal, 'pidKp') && isfinite(cal.pidKp) && cal.pidKp >= 0
                app.Kp.Value = cal.pidKp;
            end
            if isfield(cal, 'pidKi') && isfinite(cal.pidKi) && cal.pidKi >= 0
                app.Ki.Value = cal.pidKi;
            end
            if isfield(cal, 'pidKd') && isfinite(cal.pidKd) && cal.pidKd >= 0
                app.Kd.Value = cal.pidKd;
            end
            if isfield(cal, 'pidTc') && isfinite(cal.pidTc) && cal.pidTc >= 0.02
                app.CtrlPeriod.Value = cal.pidTc;
            end
            app.Message.Value = 'Loaded saved calibration factors from last session (re-zero for the current setup).';
        end

        function saveCalCache(app)
            try
                cal = struct('tare', app.tare, 'scale', app.scale, ...
                    'pressureA', app.pressureA, 'pressureB', app.pressureB, ...
                    'pidKp', app.Kp.Value, 'pidKi', app.Ki.Value, ...
                    'pidKd', app.Kd.Value, 'pidTc', app.CtrlPeriod.Value, ...
                    'saved', datetime('now'));
                save(calCachePath(app), 'cal');
            catch ME
                app.Message.Value = ['Note: could not save calibration cache (', ME.message, ').'];
            end
        end

        function restoreCalCache(app, hadCache)
            % Put back the real calibration cache after the offline self test.
            calFile = calCachePath(app);
            if hadCache
                if isfile(calFile)
                    delete(calFile);
                end
                if isfile([calFile, '.selftest_bak'])
                    movefile([calFile, '.selftest_bak'], calFile, 'f');
                end
            elseif isfile(calFile)
                delete(calFile);  % test created it; do not leave dummy factors
            end
            if isfile([calFile, '.selftest_bak'])  % belt and suspenders
                delete([calFile, '.selftest_bak']);
            end
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
                pidReset(app);
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
                app.Message.Value = 'Error: calibrate, or enter tare + scale on the Known LC Cal tab.';
                return;
            end

            app.r = 0;
            app.t = 0;
            app.measure.Value = selectedUnitText(app);
            ylabel(app.Axes1, ['Force ', app.Unit.Value]);

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

            servoWasOn = app.EnablePressureControl.Value;
            if servoWasOn
                pidReset(app);  % fresh step for each acquisition run
            end
            while (app.r == 0) && (app.time(app.i) < app.session_time) ...
                    && ((app.i - app.last_data) < sampleTarget)
                loopTimer = tic;
                raw = read_HX711(app.HX711_obj);
                grams = rawToGrams(app, raw);
                forceSelected = gramsToSelectedUnit(app, grams);
                forceNValue = gramsToNewtons(app, grams);

                voltage = readVoltage(app.a, app.pressurePin);
                p = pressureVoltageToKPa(app, voltage);

                dutyU = 0;
                if servoWasOn
                    % PID tick: derivative dt = actual time since last sample.
                    dutyU = pidCompute(app, p, max(toc(loopTimer), 0.01), ...
                        app.DesiredPressure.Value);
                end

                updateMaxLoadGaugeWithGrams(app, grams);

                app.rawForce(app.i) = raw;
                app.force(app.i) = forceSelected;
                app.forceN(app.i) = forceNValue;
                app.pressure(app.i) = p;
                app.pressureV(app.i) = voltage;

                plot(app.Axes1, app.time(1:app.i), app.force(1:app.i));
                app.Pressure.Value = app.pressure(app.i);
                app.PressureGauge.Value = max(app.PressureGauge.Limits(1), ...
                    min(app.PressureGauge.Limits(2), app.pressure(app.i)));
                app.ForceEdit.Value = app.force(app.i);
                app.TimeEdit.Value = app.time(app.i) - app.rst_time;
                app.DataEditField.Value = app.i - app.last_data;
                drawnow limitrate;

                if servoWasOn
                    % Spend the rest of this sample's tick pulsing the
                    % valves (fill/hold/vent fraction of the PID duty).
                    applyValveDuty(app, dutyU, max(0.02, app.Add_time.Value - toc(loopTimer)));
                end
                pause(max(0, app.Add_time.Value - toc(loopTimer)));
                app.Rate.Value = toc(loopTimer);
                Xaxis(app);
                app.i = app.i + 1;
                app.time(app.i) = toc(app.start_time);
            end
            if servoWasOn
                % Park on "maintain" so pressure does not run away after stop.
                writeDigitalPin(app.a, app.valveIncPin, 0);
                writeDigitalPin(app.a, app.valveMaintainPin, 1);
            end
            app.Cyan.Color = 'white';
            app.Yellow.Color = 'yellow';
            app.Message.Value = sprintf('Acquisition stopped. Samples in current run: %d.', ...
                app.i - 1 - app.last_data);
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

            try
                filePath = buildSavePath(app);
            catch ME
                app.Message.Value = ['Save error (folder?): ', ME.message];
                return;
            end
            if isfile(filePath)
                answer = uiconfirm(app.MatlabArduinoHX711UIFigure, ...
                    sprintf('%s already exists. Overwrite?', filePath), ...
                    'Overwrite existing file?', 'Options', {'Overwrite','Cancel'}, ...
                    'DefaultOption', 'Cancel', 'CancelOption', 'Cancel');
                if ~strcmp(answer, 'Overwrite')
                    app.Message.Value = 'Save cancelled.';
                    return;
                end
            end

            nRows = min(750, numel(idx));
            idx = idx(1:nRows);
            ColumnNames = {'Time_s','RawHX711_counts','Force_N','Pressure_kPa', ...
                'PressureVoltage_V','Force_SelectedUnit'};
            A = numel(ColumnNames);
            Data = NaN(750, A);
            Data(1:nRows,:) = [app.time(idx).', app.rawForce(idx).', app.forceN(idx).', ...
                app.pressure(idx).', app.pressureV(idx).', app.force(idx).'];

            Stats = buildStatsTable(app, Data(1:nRows,3), Data(1:nRows,4));

            Metadata = struct();
            Metadata.Created = datetime('now');
            Metadata.Series = round(app.SeriesSpinner.Value);
            Metadata.Run = round(app.RunSpinner.Value);
            Metadata.FileName = char(filePath);
            Metadata.FilePrefix = app.Prefix.Value;
            Metadata.NsamplesSaved = nRows;
            Metadata.DataRows = 750;
            Metadata.DataColumns = A;
            Metadata.SelectedForceUnit = selectedUnitText(app);
            Metadata.KneeAngle_deg = app.KneeAngle.Value;
            Metadata.LoadCellAngle_deg = app.LoadCellAngle.Value;
            Metadata.LoadCellTare = app.tare;
            Metadata.LoadCellScale = app.scale;
            Metadata.PressureCalibration = struct('a_kPaPerV', app.pressureA, ...
                'b_kPa', app.pressureB, 'Equation', 'Pressure_kPa = a*Voltage_V + b');
            Metadata.PressureServo = struct('Enabled', logical(app.EnablePressureControl.Value), ...
                'Setpoint_kPa', app.DesiredPressure.Value, 'Deadband_kPa', app.PressureDeadband.Value);
            Metadata.Pins = struct('Clock', app.clock, 'Data', app.data, ...
                'Pressure', app.pressurePin, 'ValveIncrease', app.valveIncPin, ...
                'ValveMaintain', app.valveMaintainPin);
            Metadata.ValveLogic = struct('Increasing', 'Increase High, Maintain High', ...
                'Maintain', 'Increase Low, Maintain High', ...
                'DecreasingOrZero', 'Increase Low, Maintain Low');
            Metadata.MatlabVersion = version;

            save(filePath, 'Data', 'Stats', 'Metadata', 'ColumnNames');
            try
                writeTxtSidecar(app, filePath, Data, nRows);
                sideNote = ' + .txt';
            catch ME2
                sideNote = [' (txt sidecar failed: ', ME2.message, ')'];
            end

            app.Name.Value = char(filePath);
            app.last_data = app.i - 1;
            app.rst_time = app.time(app.i);
            if app.RunSpinner.Value < app.RunSpinner.Limits(2)
                app.RunSpinner.Value = app.RunSpinner.Value + 1;
            else
                app.RunSpinner.Value = 0;
            end
            app.Message.Value = ['Saved: ', char(filePath), sideNote, ...
                sprintf(' (next run #%02d)', app.RunSpinner.Value)];
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
            x = zeros(1, z);
            for j = 1:z
                x(j) = read_HX711(app.HX711_obj);
                pause(1/1000);
            end
            % Zero offset ONLY. The scale factor (counts per gram) is a
            % property of the load cell and is deliberately NOT touched:
            % the normal procedure re-zeros after mounting (horizontal,
            % tied to the tibia) and must keep the measured scale.
            app.tare = mean(x);
            app.TareDisp.Value = app.tare;
            app.KnownTare.Value = num2str(app.tare, '%.6g');
            saveCalCache(app);
            if isfinite(app.scale) && app.scale > 0
                app.Message.Value = sprintf(['Tare (zero offset) completed. ', ...
                    'Scale factor preserved: %.6g counts/g.'], app.scale);
            else
                app.Message.Value = 'Tare (zero offset) completed. No scale factor yet - run Scale Factor or enter it on the Known LC Cal tab.';
            end
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
            setWeightToGramsFromKnown(app);
            if app.known_weight == 0
                app.Message.Value = 'Error: check the known weight for calibration.';
                return;
            end
            app.Message.Value = 'Please wait...';
            app.Cyan.Color = 'cyan';
            app.Yellow.Color = 'white';
            drawnow;
            z = round(app.n.Value);
            x = zeros(1, z);
            for j = 1:z
                x(j) = read_HX711(app.HX711_obj);
                pause(1/1000);
            end
            app.scale = (mean(x) - app.tare)/app.known_weight;
            app.ScaleDisp.Value = app.scale;
            app.KnownScale.Value = num2str(app.scale, '%.12g');
            saveCalCache(app);
            app.Message.Value = sprintf('Scale factor determined (%.6g counts/g). Remove the weight, mount the cell, then Tare again - the scale factor is kept.', app.scale);
            app.Cyan.Color = 'white';
            app.Yellow.Color = 'yellow';
        end

        % Button pushed function: ApplyKnownLoadCellButton
        function ApplyKnownLoadCellButtonPushed(app, event)
            % Apply each typed field independently: a blank field keeps the
            % current value, so entering just a new zero offset can never
            % erase the calibration slope (and vice versa).
            tareTxt = strtrim(app.KnownTare.Value);
            scaleTxt = strtrim(app.KnownScale.Value);
            wantTare = ~isempty(tareTxt) && ~strcmpi(tareTxt, 'nan');
            wantScale = ~isempty(scaleTxt) && ~strcmpi(scaleTxt, 'nan');
            if ~wantTare && ~wantScale
                app.Message.Value = 'Error: enter a zero offset (tare) and/or a scale factor.';
                return;
            end
            if wantTare
                tVal = str2double(tareTxt);
                if ~isfinite(tVal)
                    app.Message.Value = 'Error: zero offset must be a number (leave blank to keep the current one).';
                    return;
                end
            end
            if wantScale
                sVal = str2double(scaleTxt);
                if ~isfinite(sVal)
                    app.Message.Value = 'Error: scale factor must be a number (leave blank to keep the current one).';
                    return;
                end
                if sVal == 0
                    app.Message.Value = 'Error: load-cell scale cannot be zero.';
                    return;
                end
            end
            if wantTare
                app.tare = tVal;
                app.TareDisp.Value = app.tare;
            end
            if wantScale
                app.scale = sVal;
                app.ScaleDisp.Value = app.scale;
            end
            saveCalCache(app);
            switch true
                case wantTare && wantScale
                    app.Message.Value = sprintf('Known zero offset and scale applied (tare %.6g, scale %.6g).', app.tare, app.scale);
                case wantTare
                    if isfinite(app.scale) && app.scale > 0
                        app.Message.Value = sprintf('Zero offset applied. Scale factor preserved: %.6g counts/g.', app.scale);
                    else
                        app.Message.Value = 'Zero offset applied. No scale factor yet - enter it or run Scale Factor.';
                    end
                otherwise
                    app.Message.Value = sprintf('Scale factor applied (%.6g counts/g); zero offset unchanged (%.6g).', app.scale, app.tare);
            end
        end

        % Button pushed function: PressureCalButton (guided 7-point)
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
                prompt = sprintf(['Set the regulator near %d kPa. Enter the actual gauge kPa ', ...
                    'after the pressure stabilizes:'], setpoints(j));
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
                readings = zeros(1, z);
                app.Message.Value = sprintf('Reading pressure pin at actual %.2f kPa...', actual(j));
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
            saveCalCache(app);
            plot(app.Axes2, voltage, actual, 'o', voltage, polyval(coeff, voltage), '-');
            xlabel(app.Axes2, 'Voltage [V]');
            ylabel(app.Axes2, 'Pressure [kPa]');
            app.Message.Value = sprintf('Pressure calibration complete: y = %.6g*x %+.6g', ...
                app.pressureA, app.pressureB);
        end

        % Button pushed function: ApplyKnownPressureButton
        function ApplyKnownPressureButtonPushed(app, event)
            app.pressureA = app.PressureA.Value;
            app.pressureB = app.PressureB.Value;
            saveCalCache(app);
            app.Message.Value = sprintf('Known pressure calibration applied: y = %.6g*x %+.6g', ...
                app.pressureA, app.pressureB);
        end

        % Button pushed function: CalibrationButton (Gaussian check)
        function CalibrationButtonPushed(app, event)
            if ~app.check_connection
                app.Message.Value = 'Error: You are not connected yet.';
                return;
            end
            if ~isLoadCellCalibrated(app)
                app.Message.Value = 'Error: calibrate, or enter tare + scale on the Known LC Cal tab.';
                return;
            end
            app.Message.Value = 'Please wait...';
            app.Cyan.Color = 'cyan';
            app.Yellow.Color = 'white';
            drawnow;
            setWeightToGramsFromKnown(app);
            z = round(app.n.Value);
            x = zeros(1, z);
            for j = 1:z
                x(j) = read_HX711(app.HX711_obj);
                pause(1/1000);
            end
            x = rawToGrams(app, x);
            M = mean(x);
            S = std(x);
            app.AverageDisp.Value = M;
            app.StdDisp.Value = S;
            w = linspace(M - 4*S, M + 4*S, 5000);
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
                app.Message.Value = 'Error: calibrate, or enter tare + scale on the Known LC Cal tab.';
                return;
            end
            app.Message.Value = 'Data acquisition';
            app.Cyan.Color = 'cyan';
            app.Yellow.Color = 'white';
            drawnow;
            raw = read_HX711(app.HX711_obj);
            grams = rawToGrams(app, raw);
            updateMaxLoadGaugeWithGrams(app, grams);
            x = gramsToSelectedUnit(app, grams);
            app.Raw.Value = x;
            app.measure4.Value = selectedUnitText(app);
            app.Cyan.Color = 'white';
            app.Yellow.Color = 'yellow';
            app.Message.Value = 'Done';
        end

        % Button pushed function: Clean
        function CleanButtonPushed(app, event)
            cla(app.Axes1);
            cla(app.Axes2);
            pidReset(app);
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
            app.rst_time = 0;
            app.last_data = 0;
            app.start_time = [];
            app.get_true = false;
            app.DataEditField.Value = 0;
            app.PressureGauge.Value = 0;
            app.Pressure.Value = 0;
            app.ForceEdit.Value = 0;
            app.TimeEdit.Value = 0;
            app.Rate.Value = 0;
            app.AverageDisp.Value = 0;
            app.StdDisp.Value = 0;
            app.Message.Value = 'All cleaned up.';
        end

        % Valve-control callbacks (manual)
        function IncreasePressureButtonPushed(app, event)
            setValves(app, 1, 1, 'Valves: pressure increasing (Increase High, Maintain High).');
        end

        function MaintainPressureButtonPushed(app, event)
            setValves(app, 0, 1, 'Valves: pressure maintained (Increase Low, Maintain High).');
        end

        function DecreasePressureButtonPushed(app, event)
            setValves(app, 0, 0, 'Valves: pressure decreasing / 0 kPa (Increase Low, Maintain Low).');
        end

        % Button pushed function: RunStepTestButton (dynamic pressure cal)
        function RunStepTestButtonPushed(app, event)
            if ~app.check_connection
                app.Message.Value = 'Error: You are not connected yet.';
                return;
            end
            sp = app.DesiredPressure.Value;
            db = max(0.5, app.PressureDeadband.Value);
            Tc = max(0.02, app.CtrlPeriod.Value);
            settleHold = 2.0;   % s inside the deadband before "settled"
            maxDur = 60;        % s hard cap; click Pause to abort early

            app.pidReset();
            v0 = readVoltage(app.a, app.pressurePin);
            p0 = pressureVoltageToKPa(app, v0);
            app.Message.Value = sprintf('Dynamic pressure cal: stepping %.1f -> %.0f kPa (PID %.4g/%.4g/%.4g)...', ...
                p0, sp, app.Kp.Value, app.Ki.Value, app.Kd.Value);
            app.Cyan.Color = 'cyan';
            app.Yellow.Color = 'white';
            drawnow;

            % FILL (+1), HOLD (0), VENT (-1)
            tLog = [];
            pLog = [];
            vLog = [];
            uLog = [];
            sLog = [];
            settleSince = [];
            settled = false;
            t0 = tic;
            while toc(t0) < maxDur && app.r == 0 && ~settled
                v = readVoltage(app.a, app.pressurePin);
                p = pressureVoltageToKPa(app, v);
                u = pidCompute(app, p, Tc, sp);
                if abs(p - sp) <= db
                    if isempty(settleSince)
                        settleSince = toc(t0);
                    end
                else
                    settleSince = [];
                end
                if ~isempty(settleSince) && (toc(t0) - settleSince) >= settleHold
                    settled = true;
                end
                if settled
                    st = applyValveDuty(app, 0, Tc);  % park on HOLD
                else
                    st = applyValveDuty(app, u, Tc);
                end
                tLog(end+1) = toc(t0);   %#ok<AGROW>
                pLog(end+1) = p;         %#ok<AGROW>
                vLog(end+1) = v;         %#ok<AGROW>
                uLog(end+1) = u;         %#ok<AGROW>
                sLog(end+1) = st;        %#ok<AGROW>
                if mod(numel(tLog), 5) == 0 || settled
                    cla(app.Axes2);
                    hold(app.Axes2, 'on');
                    plot(app.Axes2, tLog, pLog, 'b-', 'LineWidth', 1.5);
                    plot(app.Axes2, [tLog(1) tLog(end)], [sp sp], 'k--');
                    plot(app.Axes2, [tLog(1) tLog(end)], [sp+db sp+db], 'r:');
                    plot(app.Axes2, [tLog(1) tLog(end)], [sp-db sp-db], 'r:');
                    hold(app.Axes2, 'off');
                    xlabel(app.Axes2, 't [s]');
                    ylabel(app.Axes2, 'Pressure [kPa]');
                    title(app.Axes2, sprintf('Step %.1f -> %.0f kPa', p0, sp));
                    drawnow limitrate;
                end
            end
            app.r = 0;  % a Pause press aborts this test, not the next Get Data
            writeDigitalPin(app.a, app.valveIncPin, 0);
            writeDigitalPin(app.a, app.valveMaintainPin, 1);  % park HOLD

            % Metrics relative to the step actually commanded.
            stepSize = abs(sp - p0);
            pmax = max(pLog);
            pmin = min(pLog);
            if sp >= p0  % upward step
                ovKPa = max(0, pmax - sp);
                ic = find(pLog >= sp, 1);
                if isempty(ic)
                    unKPa = NaN;
                else
                    unKPa = max(0, sp - min(pLog(ic:end)));
                end
            else         % downward step
                ovKPa = max(0, sp - pmin);
                ic = find(pLog <= sp, 1);
                if isempty(ic)
                    unKPa = NaN;
                else
                    unKPa = max(0, max(pLog(ic:end)) - sp);
                end
            end
            outIdx = find(abs(pLog - sp) > db, 1, 'last');
            if isempty(outIdx)
                settleTime = 0;
            else
                settleTime = tLog(min(outIdx + 1, numel(tLog)));
            end
            ess = pLog(end) - sp;
            if isempty(settleSince)
                settled = false;
            end

            fprintf(['DPC step %.1f -> %.0f kPa (Kp %.4g, Ki %.4g, Kd %.4g, Tc %.3f s)\n', ...
                '  reached target: %d   overshoot: %.2f kPa (%.1f%% of step)\n', ...
                '  undershoot: %.2f kPa (%.1f%% of step)   settling time: %.2f s\n', ...
                '  steady-state error: %+.2f kPa   duration: %.2f s   settled: %d\n'], ...
                p0, sp, app.Kp.Value, app.Ki.Value, app.Kd.Value, Tc, ...
                ~isempty(ic), ovKPa, 100*ovKPa/max(stepSize, eps), ...
                unKPa, 100*unKPa/max(stepSize, eps), settleTime, ess, tLog(end), settled);

            % Save the trace + gains + metrics for tuning comparisons.
            DPC_Data = [tLog.', pLog.', vLog.', uLog.', sLog.'];
            DPC_ColumnNames = {'Time_s','Pressure_kPa','Voltage_V','Duty_pct','ValveState'};
            DPC_Metrics = struct('StartPressure_kPa', p0, 'Setpoint_kPa', sp, ...
                'StepSize_kPa', sp - p0, 'Overshoot_kPa', ovKPa, ...
                'Overshoot_pctStep', 100*ovKPa/max(stepSize, eps), ...
                'Undershoot_kPa', unKPa, 'Undershoot_pctStep', 100*unKPa/max(stepSize, eps), ...
                'RiseTime10to90_s', rise10to90(tLog, pLog, p0, sp), ...
                'SettlingTime_s', settleTime, 'Ess_kPa', ess, ...
                'Settled', logical(settled), 'Duration_s', tLog(end), ...
                'Deadband_kPa', db, 'SettleHold_s', settleHold);
            DPC_PID = struct('Kp', app.Kp.Value, 'Ki', app.Ki.Value, ...
                'Kd', app.Kd.Value, 'CtrlPeriod_s', Tc, 'MinDuty_pct', 2, ...
                'ValveStates', 'FILL=+1 (D11 High,D6 High); HOLD=0 (D11 Low,D6 High); VENT=-1 (D11 Low,D6 Low)');
            try
                saveDir = strtrim(char(string(app.SaveFolder.Value)));
                if isempty(saveDir)
                    saveDir = app.resolveDefaultSaveDir();
                end
                if ~isfolder(saveDir)
                    mkdir(saveDir);
                end
                app.dpcRun = app.dpcRun + 1;
                dpcPath = fullfile(saveDir, sprintf('DPC_S%02d_R%02d.mat', ...
                    round(app.SeriesSpinner.Value), app.dpcRun));
                save(dpcPath, 'DPC_Data', 'DPC_Metrics', 'DPC_PID', 'DPC_ColumnNames');
                saveNote = [' Saved: ', dpcPath];
            catch ME2
                saveNote = [' (DPC save failed: ', ME2.message, ')'];
            end
            app.Message.Value = sprintf(['Step test done: overshoot %.2f kPa (%.1f%%), undershoot %.2f kPa, ', ...
                'settle %.2f s, ess %+.2f kPa.%s'], ovKPa, 100*ovKPa/max(stepSize, eps), ...
                unKPa, settleTime, ess, saveNote);
            app.Cyan.Color = 'white';
            app.Yellow.Color = 'yellow';
        end
    end

    % Public: hardware-free self test (used by test_HX711_BPA_offline.m)
    methods (Access = public)
        function res = offlineSelfTest(app, tmpDir)
            if nargin < 2
                tmpDir = fullfile(tempdir, 'hx711_bpa_selftest');
            end
            if ~isfolder(tmpDir)
                mkdir(tmpDir);
            end
            res = struct();
            res.allPass = true;

            % Preserve any real calibration cache: this test writes dummy
            % factors, so stash the cache and restore it at the end.
            calFile = calCachePath(app);
            hadCache = isfile(calFile);
            if hadCache
                movefile(calFile, [calFile, '.selftest_bak'], 'f');
            end
            cleanupObj = onCleanup(@() restoreCalCache(app, hadCache));

            % 1) Guarded callbacks before Connect must show the error message.
            guarded = {@GetDataButtonPushed, @TareButtonPushed, @ScaleFactorButtonPushed, ...
                @CalibrationButtonPushed, @RawReadButtonPushed, @IncreasePressureButtonPushed, ...
                @MaintainPressureButtonPushed, @DecreasePressureButtonPushed, @RunStepTestButtonPushed};
            res.guardsPass = true;
            for k = 1:numel(guarded)
                guarded{k}(app);
                if ~startsWith(app.Message.Value, 'Error')
                    res.guardsPass = false;
                end
            end
            res.allPass = res.allPass && res.guardsPass;
            fprintf('guards pre-connect: %d (%d callbacks checked)\n', ...
                res.guardsPass, numel(guarded));

            % 2) Known-factor entry unlocks calibration without hardware.
            app.KnownTare.Value = '1234';
            app.KnownScale.Value = '98.7';
            ApplyKnownLoadCellButtonPushed(app);
            res.knownLC = abs(app.tare - 1234) < 1e-9 && abs(app.scale - 98.7) < 1e-9 ...
                && isLoadCellCalibrated(app);
            res.allPass = res.allPass && res.knownLC;
            fprintf('apply known LC factors: %d\n', res.knownLC);

            % 2b) Re-zeroing must preserve the scale factor (normal
            % procedure: hang -> tare -> weight -> scale factor -> mount
            % horizontally -> tare again).
            app.KnownTare.Value = '2222';
            app.KnownScale.Value = '';  % blank slope = keep current scale
            ApplyKnownLoadCellButtonPushed(app);
            res.retareKeepsScale = abs(app.tare - 2222) < 1e-9 ...
                && abs(app.scale - 98.7) < 1e-9 && isLoadCellCalibrated(app);
            S = load(calCachePath(app), 'cal');
            res.retareKeepsScale = res.retareKeepsScale ...
                && abs(S.cal.tare - 2222) < 1e-9 && abs(S.cal.scale - 98.7) < 1e-9;
            % Garbage in a field must be rejected without touching state.
            app.KnownScale.Value = 'abc';
            ApplyKnownLoadCellButtonPushed(app);
            res.retareKeepsScale = res.retareKeepsScale ...
                && startsWith(app.Message.Value, 'Error') ...
                && abs(app.scale - 98.7) < 1e-9;
            res.allPass = res.allPass && res.retareKeepsScale;
            fprintf('re-tare keeps scale factor (app state + cache, garbage rejected): %d\n', ...
                res.retareKeepsScale);

            % 3) Conversion math (relative to the current tare/scale state,
            % which the re-zero check above intentionally changed).
            g1 = rawToGrams(app, app.tare + app.scale*1000);
            res.mathPass = abs(g1 - 1000) < 1e-6 && abs(gramsToNewtons(app, 1000) - 9.80665) < 1e-3;
            app.PressureA.Value = 156.04;
            app.PressureB.Value = -128.2;
            ApplyKnownPressureButtonPushed(app);
            res.mathPass = res.mathPass && abs(pressureVoltageToKPa(app, 1.0) - 27.84) < 1e-6;
            res.allPass = res.allPass && res.mathPass;
            fprintf('conversion math: %d\n', res.mathPass);

            % 3b) PID math: proportional response, saturation, anti-windup,
            % derivative on measurement, and the duty->valve timing map.
            app.Kp.Value = 2; app.Ki.Value = 0; app.Kd.Value = 0;
            pidReset(app);
            uSat = pidCompute(app, 0, 0.1, 100);   % huge error -> +100 rail
            uP = pidCompute(app, 60, 0.1, 100);    % (100-60)*2 = 80
            pidReset(app);
            app.Kp.Value = 0; app.Ki.Value = 1; app.Kd.Value = 0;
            uI1 = pidCompute(app, 90, 0.1, 100);   % e=10 -> u = 1*10*0.1 = 1
            uI2 = pidCompute(app, 90, 0.1, 100);   % integral accumulates -> 2
            app.pidI = 0;
            app.Ki.Value = 10;
            for kk = 1:50
                uW = pidCompute(app, 0, 0.1, 100); % saturated, e>0: no windup
            end
            app.Kp.Value = 0; app.Ki.Value = 0; app.Kd.Value = 1;
            pidReset(app);
            uD1 = pidCompute(app, 90, 0.1, 100);   % first tick: no derivative
            uD2 = pidCompute(app, 95, 0.1, 100);   % dp/dt = +50 -> D = -50
            [stF, tfF, ~, thF] = dutyPlan(app, 50, 0.1);
            [stV, ~, tvV, thV] = dutyPlan(app, -30, 0.1);
            [stH, ~, ~, thH] = dutyPlan(app, 1, 0.1);
            res.pidPass = (uSat == 100) && abs(uP - 80) < 1e-9 ...
                && abs(uI1 - 1) < 1e-9 && abs(uI2 - 2) < 1e-9 ...
                && app.pidI == 0 && abs(uW) <= 100 ...
                && abs(uD1) < 1e-12 && abs(uD2 + 50) < 1e-9 ...
                && stF == 1 && abs(tfF - 0.05) < 1e-12 && abs(thF - 0.05) < 1e-12 ...
                && stV == -1 && abs(tvV - 0.03) < 1e-12 && abs(thV - 0.07) < 1e-12 ...
                && stH == 0 && abs(thH - 0.1) < 1e-12;
            res.allPass = res.allPass && res.pidPass;
            fprintf('PID math + duty map: %d\n', res.pidPass);

            % 4) Save path, MAT content, txt sidecar, run auto-increment.
            app.SaveFolder.Value = tmpDir;
            app.Prefix.Value = 'SelfTest';
            app.SeriesSpinner.Value = 2;
            app.RunSpinner.Value = 3;
            % Clear leftovers from a previous (possibly crashed) run so the
            % overwrite prompt - which cannot run headless - never triggers.
            delete(fullfile(tmpDir, 'SelfTest2_03.*'));
            ns = 5;
            app.rawForce = (1:ns)*100 + app.tare;
            app.force = (1:ns);
            app.forceN = (1:ns)*app.g/1000;
            app.pressure = 100 + (1:ns);
            app.pressureV = 1 + (1:ns)/app.pressureA;
            app.time = (0:ns)*0.1;  % ns+1 entries: post-acquisition state has time(i)
            app.i = ns + 1;
            app.get_true = true;
            app.last_data = 0;
            SaveButtonPushed(app);
            matPath = fullfile(tmpDir, 'SelfTest2_03.mat');
            txtPath = fullfile(tmpDir, 'SelfTest2_03.txt');
            res.savePass = isfile(matPath) && isfile(txtPath);
            if res.savePass
                S = load(matPath);
                lines = readlines(txtPath);
                lines = lines(strlength(lines) > 0);  % ignore trailing empty line
                chk = struct();
                chk.sizeOk = isequal(size(S.Data), [750 6]);
                chk.colsOk = isequal(S.ColumnNames, {'Time_s','RawHX711_counts','Force_N', ...
                    'Pressure_kPa','PressureVoltage_V','Force_SelectedUnit'});
                chk.rawOk = abs(S.Data(1,2) - (100 + app.tare)) < 1e-6;
                chk.calOk = abs(S.Metadata.LoadCellScale - 98.7) < 1e-9;
                chk.seriesRunOk = S.Metadata.Series == 2 && S.Metadata.Run == 3;
                chk.txtHeadOk = contains(lines(1), 'Force [ N ]');
                chk.txtRowsOk = height(lines) == ns + 1;
                res.savePass = all(structfun(@(v) logical(v), chk));
                if ~res.savePass
                    disp(chk);  % which sub-check failed
                end
            end
            res.runIncrement = app.RunSpinner.Value == 4;
            res.allPass = res.allPass && res.savePass && res.runIncrement;
            fprintf('save mat+txt + run increment: %d / %d\n', res.savePass, res.runIncrement);

            % 5) Default save folder resolves relative to the app file.
            defDir = app.resolveDefaultSaveDir();
            res.saveDirPass = (contains(defDir, 'Testing_Data') && isfolder(fileparts(defDir))) ...
                || endsWith(defDir, 'saved_data');
            res.allPass = res.allPass && res.saveDirPass;
            fprintf('default save folder (%s): %d\n', defDir, res.saveDirPass);

            % 6) Clean resets acquisition state.
            CleanButtonPushed(app);
            res.cleanPass = app.get_true == false && app.i == 1 && isempty(app.time);
            res.allPass = res.allPass && res.cleanPass;
            fprintf('clean resets state: %d\n', res.cleanPass);

            % Remove self-test artifacts (the cal cache is restored/removed
            % by the onCleanup above when this method returns).
            delete(fullfile(tmpDir, 'SelfTest2_03.*'));

            if res.allPass
                app.Message.Value = 'Offline self test PASSED.';
            else
                app.Message.Value = 'Offline self test FAILED - see command window.';
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)
            app.MatlabArduinoHX711UIFigure = uifigure('Visible', 'off');
            app.MatlabArduinoHX711UIFigure.Color = [0.9412 0.9412 0.9412];
            app.MatlabArduinoHX711UIFigure.Position = [100 100 1120 900];
            app.MatlabArduinoHX711UIFigure.Name = 'HX711_BPA - BPA Force & Pressure Test App (AARL)';
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
            app.RawRead.Text = 'Raw Read';

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
            app.ForceDropDownLabel = uilabel(app.DataAcquisitionTab, 'Text', 'Force', 'Position', [4 152 32 15]);
            app.Unit = uidropdown(app.DataAcquisitionTab, 'Items', {'[ N ]','[ g ]','[ kg ]','[ kN ]'}, 'Position', [44 148 76 22], 'Value', '[ N ]');
            app.SetSession = uicheckbox(app.DataAcquisitionTab, 'Text', 'Set Session Time [min]', 'Position', [4 128 170 15]);
            app.SessionTimeSpinnerLabel = uilabel(app.DataAcquisitionTab, 'Text', 'Session Time', 'Position', [4 106 70 15]);
            app.SessionV = uispinner(app.DataAcquisitionTab, 'Limits', [0 Inf], 'Position', [88 102 60 22]);
            app.SamplingRateSpinnerLabel = uilabel(app.DataAcquisitionTab, 'Text', 'Sample Period [s]', 'Position', [4 82 86 15]);
            app.Add_time = uispinner(app.DataAcquisitionTab, 'Step', 0.1, 'Limits', [0.001 Inf], 'ValueDisplayFormat', '%.3f', 'Position', [96 78 52 22], 'Value', 0.1);
            app.SampleCountSpinnerLabel = uilabel(app.DataAcquisitionTab, 'Text', 'Samples', 'Position', [4 58 44 15]);
            app.Nsamples = uispinner(app.DataAcquisitionTab, 'Limits', [1 750], 'ValueDisplayFormat', '%.0f', 'Position', [84 54 60 22], 'Value', 750);

            app.PressureCtrlTab = uitab(app.TabGroup);
            app.PressureCtrlTab.Title = 'Pressure Ctrl';
            app.EnablePressureControl = uicheckbox(app.PressureCtrlTab, 'Text', 'PID servo during Get Data', 'Position', [4 158 178 15], 'Value', false);
            app.KpLabel = uilabel(app.PressureCtrlTab, 'Text', 'Kp  [%/kPa]', 'Position', [4 136 70 15]);
            app.Kp = uieditfield(app.PressureCtrlTab, 'numeric', 'Limits', [0 Inf], 'ValueDisplayFormat', '%.4g', 'Position', [96 132 50 22], 'Value', 2);
            app.KiLabel = uilabel(app.PressureCtrlTab, 'Text', 'Ki  [%/kPa/s]', 'Position', [4 112 74 15]);
            app.Ki = uieditfield(app.PressureCtrlTab, 'numeric', 'Limits', [0 Inf], 'ValueDisplayFormat', '%.4g', 'Position', [96 108 50 22], 'Value', 0.5);
            app.KdLabel = uilabel(app.PressureCtrlTab, 'Text', 'Kd  [%-s/kPa]', 'Position', [4 88 74 15]);
            app.Kd = uieditfield(app.PressureCtrlTab, 'numeric', 'Limits', [0 Inf], 'ValueDisplayFormat', '%.4g', 'Position', [96 84 50 22], 'Value', 0);
            app.SetpointkPaLabel = uilabel(app.PressureCtrlTab, 'Text', 'Set / +-dB [kPa]', 'Position', [4 64 84 15]);
            app.DesiredPressure = uispinner(app.PressureCtrlTab, 'Limits', [0 700], 'ValueDisplayFormat', '%.0f', 'Position', [96 60 40 22], 'Value', 400);
            app.PressureDeadband = uispinner(app.PressureCtrlTab, 'Limits', [0 100], 'ValueDisplayFormat', '%.0f', 'Position', [140 60 40 22], 'Value', 5);
            app.CtrlPeriodLabel = uilabel(app.PressureCtrlTab, 'Text', 'Ctrl Period [s]', 'Position', [4 40 80 15]);
            app.CtrlPeriod = uispinner(app.PressureCtrlTab, 'Limits', [0.02 5], 'Step', 0.01, 'ValueDisplayFormat', '%.2f', 'Position', [96 36 50 22], 'Value', 0.10);
            app.RunStepTestButton = uibutton(app.PressureCtrlTab, 'push', 'Text', 'Run Step Test (Dynamic Cal)', 'Position', [4 8 176 26]);
            app.RunStepTestButton.ButtonPushedFcn = createCallbackFcn(app, @RunStepTestButtonPushed, true);
            app.RunStepTestButton.FontWeight = 'bold';

            app.SaveDataTab = uitab(app.TabGroup);
            app.SaveDataTab.Title = 'Save Data';
            app.SaveFolderLabel = uilabel(app.SaveDataTab, 'Text', 'Folder', 'Position', [7 146 40 15]);
            app.SaveFolder = uieditfield(app.SaveDataTab, 'text', 'Position', [52 142 129 22]);
            app.SeriesSpinnerLabel = uilabel(app.SaveDataTab, 'Text', 'Series $', 'Position', [7 112 52 15]);
            app.SeriesSpinner = uispinner(app.SaveDataTab, 'Limits', [0 999], 'ValueDisplayFormat', '%.0f', 'Position', [72 108 45 22], 'Value', 1);
            app.RunSpinnerLabel = uilabel(app.SaveDataTab, 'Text', 'Run #', 'Position', [7 80 42 15]);
            app.RunSpinner = uispinner(app.SaveDataTab, 'Limits', [0 99], 'ValueDisplayFormat', '%02.0f', 'Position', [72 76 45 22], 'Value', 0);
            app.PrefixLabel = uilabel(app.SaveDataTab, 'Text', 'Prefix', 'Position', [7 48 40 15]);
            app.Prefix = uieditfield(app.SaveDataTab, 'text', 'Position', [72 44 109 22], 'Value', 'FlxTest');
            app.NameEditFieldLabel = uilabel(app.SaveDataTab, 'Text', 'Last file', 'Position', [7 16 50 15]);
            app.Name = uieditfield(app.SaveDataTab, 'text', 'Editable', 'off', 'Position', [60 12 121 22]);

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
            app.Axes2.XGrid = 'on';
            app.Axes2.YGrid = 'on';
            app.Axes2.Box = 'on';
            app.Axes2.Position = [8 8 350 265];
            app.TabGroup2 = uitabgroup(app.CalibrationResultPanel);
            app.TabGroup2.Position = [367 8 344 257];

            app.ValueTab = uitab(app.TabGroup2);
            app.ValueTab.Title = 'Value';
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

            app.SettingsTab = uitab(app.TabGroup2);
            app.SettingsTab.Title = 'Settings';
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

            app.KnownCalTab = uitab(app.TabGroup2);
            app.KnownCalTab.Title = 'Known LC Cal';
            app.KnownTareLabel = uilabel(app.KnownCalTab, 'Text', 'Zero offset / tare [counts]', 'Position', [10 190 180 15]);
            app.KnownTare = uieditfield(app.KnownCalTab, 'text', 'Position', [200 186 120 22]);
            app.KnownScaleLabel = uilabel(app.KnownCalTab, 'Text', 'Calibration slope [counts/g] (blank = keep)', 'Position', [10 150 250 15]);
            app.KnownScale = uieditfield(app.KnownCalTab, 'text', 'Position', [200 146 120 22]);
            app.ApplyKnownLoadCellButton = uibutton(app.KnownCalTab, 'push', 'Text', 'Apply Known Load-Cell Cal', 'Position', [70 104 200 28]);
            app.ApplyKnownLoadCellButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyKnownLoadCellButtonPushed, true);
            app.KnownCalHint = uilabel(app.KnownCalTab);
            app.KnownCalHint.Position = [10 16 324 78];
            app.KnownCalHint.FontSize = 9;
            app.KnownCalHint.FontColor = [0.35 0.35 0.35];
            app.KnownCalHint.WordWrap = 'on';
            app.KnownCalHint.Text = ['Normal procedure: hang the cell, Tare; tie on a known weight, Scale Factor; ', ...
                'remove the weight; mount horizontally tied to the tibia and Tare again. ', ...
                'Taring only updates the zero offset - the scale factor is kept. ', ...
                'A field left blank keeps its current value when you click Apply. ', ...
                'Factors persist across sessions.'];

            app.PressureCalTab = uitab(app.TabGroup2);
            app.PressureCalTab.Title = 'Pressure Cal';
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
            app.PressureCalHint = uilabel(app.PressureCalTab);
            app.PressureCalHint.Position = [10 4 324 36];
            app.PressureCalHint.FontSize = 9;
            app.PressureCalHint.FontColor = [0.35 0.35 0.35];
            app.PressureCalHint.WordWrap = 'on';
            app.PressureCalHint.Text = ['Pressure_kPa = a*Voltage_V + b. Type known a and b and Apply, ', ...
                'or run the guided 7-point calibration (regulator 0-620 kPa).'];

            app.LicenseTab = uitab(app.TabGroup2);
            app.LicenseTab.Title = '*License*';
            app.TextArea = uitextarea(app.LicenseTab, 'Editable', 'off', 'Position', [10 9 320 215]);
            app.TextArea.Value = {'Original HX711 app: copyright 2018, Nicholas Giacoboni (BSD).', ...
                'Pressure/encoder/valve customization + known-factor calibration entry:', ...
                'AARL / Bipedal_Robot project. See README_HX711_BPA.md.'};

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
        end
    end

    % App creation and deletion
    methods (Access = public)

        function app = HX711_BPA(appRoot, keepHidden)
            % HX711_BPA(appRoot, keepHidden) - both arguments optional.
            % appRoot anchors the default save folder and the calibration
            % cache; it defaults to the folder containing this file.
            if nargin < 1 || isempty(appRoot)
                appRoot = fileparts(mfilename('fullpath'));
            end
            if nargin < 2
                keepHidden = false;
            end
            app.appRoot = char(appRoot);
            createComponents(app);
            app.SaveFolder.Value = app.resolveDefaultSaveDir();
            app.loadCalCache();
            registerApp(app, app.MatlabArduinoHX711UIFigure);
            if ~keepHidden
                app.MatlabArduinoHX711UIFigure.Visible = 'on';
            end
            if nargout == 0
                clear app
            end
        end

        function delete(app)
            try
                if ~isempty(app.a) && isvalid(app.a)
                    app.a.disconnect();
                end
            catch
                % never block app teardown on disconnect problems
            end
            delete(app.MatlabArduinoHX711UIFigure);
        end
    end
end

function rt = rise10to90(t, p, p0, sp)
    % 10-90 % rise time of the pressure step response; NaN when the step
    % size is 0 or the target was not reached. Local function of this
    % class file (callable from the class methods by plain name).
    stepSize = sp - p0;
    if abs(stepSize) < eps
        rt = NaN;
        return;
    end
    i10 = find((p - p0)*sign(stepSize) >= 0.10*abs(stepSize), 1);
    i90 = find((p - p0)*sign(stepSize) >= 0.90*abs(stepSize), 1);
    if isempty(i10) || isempty(i90)
        rt = NaN;
    else
        rt = t(i90) - t(i10);
    end
end
