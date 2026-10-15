% Reproduce readserialnumbers2 display + parse locally (no hardware).
td = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Testing_Data/2026_06_Festo';

% 1) parse check: both firmware line formats -> same 5 numbers
lineOld = "123456,35.250,12.345,80.500,1,0";   % 6-col encoder sketch
lineNew = "123456,12.345,80.500,1,0";          % 5-col new sketch
p6 = split(lineOld, ","); p5 = split(lineNew, ",");
v6 = str2double(p6([1 3 4 5 6])).';
v5 = str2double(p5).';
assert(isequal(v6, v5), 'parse mismatch');
fprintf('PARSE CHECK PASS: 6-col and 5-col lines give the same 5 values\n');

% 2) figure reproduction: empty state + filled state, rendered to PNG
ss = get(groot, "ScreenSize");
figW = 780; figH = 660;
figPos = [max(20, floor((ss(3)-figW)/2)), max(40, floor((ss(4)-figH)/2)), figW, figH];
fig = figure("Name", "Arduino Live Data", "NumberTitle", "off", ...
    "Position", figPos, "Visible", "off");
movegui(fig, "center");

layout = tiledlayout(fig, 3, 1, "TileSpacing", "compact", "Padding", "compact");
axForce = nexttile(layout);
forceLine = plot(axForce, NaN, NaN); ylabel(axForce, "Force (N)"); grid(axForce, "on");
axPressure = nexttile(layout);
pressureLine = plot(axPressure, NaN, NaN); ylabel(axPressure, "Pressure (kPa)"); grid(axPressure, "on");
axTorque = nexttile(layout);
torqueLine = plot(axTorque, NaN, NaN, "LineWidth", 1.5); hold(axTorque, "on");
targetLine = plot(axTorque, NaN, NaN, "r--", "LineWidth", 1.5); hold(axTorque, "off");
ylabel(axTorque, "Knee torque (N\cdotm)");
xlabel(axTorque, "Time relative to latest sample (s)"); grid(axTorque, "on");
statusTitle = title(layout, "V = valves on/start | S = save | O = valves off | Q = quit");
uKneeEdit = uicontrol(fig, "Style", "edit", "String", "-30", "Units", "normalized", "Position", [0.01 0.010 0.07 0.040]);
uLCEdit = uicontrol(fig, "Style", "edit", "String", "30", "Units", "normalized", "Position", [0.15 0.010 0.07 0.040]);
uVerdict = uicontrol(fig, "Style", "text", "String", "Enter knee + load-cell angles for the torque guard", ...
    "Units", "normalized", "Position", [0.35 0.010 0.63 0.045], "HorizontalAlignment", "left", "FontWeight", "bold");

exportgraphics(fig, 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/tmp/matlab_run/display_empty.png', 'Resolution', 110);

% synthetic 60 s of rows: force ramps 0..60 N, torque guard mock
t = 0:0.5:60;
F = 60*t/max(t);
liveData = [t'*1000, F', 400+20*sin(t'), ones(size(t')), zeros(size(t')), ...
    repmat(-30, size(t')), repmat(30, size(t')), 0.287.*F', repmat(49.84, size(t'))];
relativeTime = (liveData(:,1) - liveData(end,1))/1000;
set(forceLine, "XData", relativeTime, "YData", liveData(:,2));
set(pressureLine, "XData", relativeTime, "YData", liveData(:,3));
set(torqueLine, "XData", relativeTime, "YData", liveData(:,7));
set(targetLine, "XData", relativeTime([1 end]), "YData", [49.84 49.84]);
statusTitle.String = "NOT RECORDING | Fill = 1 | Exhaust = 0\nForce = 60.000 N | Pressure = 400.123 kPa";
uVerdict.String = "MET: 17.20 N*m vs human 49.84 N*m (+0.00 margin)";
uVerdict.String = sprintf('BELOW: %.2f N*m vs human %.2f N*m (%.2f short)', 0.287*60, 49.84, 0.287*60-49.84);
uVerdict.ForegroundColor = [0.8 0 0];
exportgraphics(fig, 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/tmp/matlab_run/display_filled.png', 'Resolution', 110);
disp('DISPLAY PROBE PASS');
