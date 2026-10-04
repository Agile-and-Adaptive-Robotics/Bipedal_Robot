% Verify live-gauge + valve-panel changes.
here = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Code/Matlab/HX711 v3.0/HX711 v3.0';
addpath(here);
app = HX711_BPA(here, true);
fig = app.MatlabArduinoHX711UIFigure;

% live timer: created, not running, 5 Hz
assert(~isempty(app.PressureTimer) && ~strcmp(app.PressureTimer.Running,'on'), 'timer state');
assert(abs(app.PressureTimer.Period - 0.2) < 1e-9, 'timer period');

% valve panel: at bottom-right, between readouts and gauge, correct labels
vp = app.ValvePanel.Position;
g = app.PressureGauge.Position;
f = app.ForceEdit.Position;
assert(vp(1) > f(1)+f(3) && vp(1)+vp(3) < g(1), 'valve panel not between readouts and gauge');
assert(strcmp(app.IncreasePressureButton.Text,'Open (Fill)'), 'fill label');
assert(strcmp(app.MaintainPressureButton.Text,'Hold'), 'hold label');
assert(strcmp(app.DecreasePressureButton.Text,'Deflate (Vent)'), 'vent label');

% panels pairwise no-overlap (incl. moved ValvePanel)
P = [app.GlobalSettingsPanel.Position; app.CalibrationPanel.Position; ...
     app.CleanPanel.Position; app.ValvePanel.Position; ...
     app.StatusPanel.Position; app.ArduinoHX711Panel.Position];
names = {'GlobalSettings','Calibration','Clean','Valves','Status','Arduino'};
for i = 1:size(P,1)
    for j = i+1:size(P,1)
        ox = max(0, min(P(i,1)+P(i,3), P(j,1)+P(j,3)) - max(P(i,1), P(j,1)));
        oy = max(0, min(P(i,2)+P(i,4), P(j,2)+P(j,4)) - max(P(i,2), P(j,2)));
        assert(ox == 0 || oy == 0, 'panels overlap: %s vs %s', names{i}, names{j});
    end
end
fprintf('live timer + valve panel OK (timer %.1f Hz, panel [%s])\n', ...
    1/app.PressureTimer.Period, mat2str(round(vp)));
delete(app);

cd(here);
test_HX711_BPA_offline;
fprintf('V3 TEST PASS\n');
