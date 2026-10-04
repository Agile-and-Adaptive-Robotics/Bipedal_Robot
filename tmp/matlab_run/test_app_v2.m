% Combined verification for the HX711_BPA UI changes + Arduino diagnostics.
here = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Code/Matlab/HX711 v3.0/HX711 v3.0';
addpath(here);

% 1) syntax check
r = checkcode(fullfile(here,'HX711_BPA.m'));
fprintf('checkcode: %d messages (constructor below is the hard gate)\n', numel(r));

% 2) construct + layout assertions
app = HX711_BPA(here, true);
fig = app.MatlabArduinoHX711UIFigure;
ip = fig.InnerPosition;   % client size AFTER shrink-to-fit scaling
g = app.PressureGauge.Position;
assert(g(1)+g(3) >= 0.9*ip(3) && g(2) < 0.3*ip(4), 'gauge not bottom-right: %s', mat2str(g));
f = app.ForceEdit.Position;
assert(f(1) >= 0.3*ip(3) && f(1)+f(3) < g(1), 'Force readout not left of gauge');
gs = app.GlobalSettingsPanel.Position;
assert(gs(3) > 300 && gs(4) > 250, 'GlobalSettings not bigger: %s', mat2str(gs));
cp = app.CalibrationPanel.Position;
cl = app.CleanPanel.Position;
% pairwise overlap among figure-level panels
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
assert(strcmp(app.GoButton.Text, 'GO'), 'Go button missing');
assert(strcmp(app.TakeCalPointButton.Text, 'Take Reading'), 'Take button missing');
fprintf('layout: gauge[%s] readouts<%d GSettings[%s] Cal[%s] Clean[%s] -- OK\n', ...
    mat2str(g), g(1), mat2str(gs), mat2str(cp), mat2str(cl));
delete(app);

% 3) offline regression
cd(here);
test_HX711_BPA_offline;

% 4) Arduino diagnostics (why Connect fails)
fprintf('\n--- Arduino diagnostics ---\n');
try
    ports = serialportlist('available');
    if isempty(ports)
        fprintf('COM ports visible to MATLAB: NONE (board not enumerated / not plugged in / driver missing)\n');
    else
        fprintf('COM ports visible to MATLAB: %s\n', strjoin(string(ports), ', '));
    end
catch e
    fprintf('serialportlist failed: %s\n', e.message);
end
try
    addons = matlab.addons.installedAddons;
    hit = addons(contains(lower({addons.Name}), 'arduino'), :);
    if isempty(hit)
        fprintf('Arduino support package: NOT INSTALLED in this MATLAB (Connect cannot work)\n');
    else
        for k = 1:height(hit)
            fprintf('Arduino support package: %s v%s\n', hit.Name{k}, hit.Version{k});
        end
    end
catch e
    fprintf('addon check failed: %s\n', e.message);
end
fprintf('DIAG DONE\n');
