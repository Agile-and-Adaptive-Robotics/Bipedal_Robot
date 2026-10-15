% Verify Global Settings tab layout: tabgroup inside panel, contents visible.
here = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Code/Matlab/HX711 v3.0/HX711 v3.0';
addpath(here);
app = HX711_BPA(here, true);

tg = app.TabGroup.Position;            % panel-relative
assert(tg(2) >= 0 && tg(2)+tg(4) <= app.GlobalSettingsPanel.Position(4)+1, ...
    'TabGroup not fully inside panel: %s', mat2str(tg));
fprintf('TabGroup panel-relative: [%s] -- inside panel OK\n', mat2str(round(tg)));

tabs = [app.ConnectionTab, app.DataAcquisitionTab, app.PressureCtrlTab, ...
        app.SaveDataTab, app.MetadataTab];
names = {'Connection','DataAcq','PressureCtrl','SaveData','Metadata'};
s = 0.87;  % approx scale; read real scale from a known pair instead:
sReal = app.GoButton.Position(4)/26;
for t = 1:numel(tabs)
    kids = tabs(t).Children;
    ymin = inf; ymax = 0;
    for k = 1:numel(kids)
        p = kids(k).Position;
        if isprop(kids(k),'Position') && ~isempty(p)
            ymin = min(ymin, p(2)/sReal);
            ymax = max(ymax, (p(2)+p(4))/sReal);
        end
    end
    fprintf('%-13s contents y-range (design px): %.0f..%.0f (%d controls)\n', ...
        names{t}, ymin, ymax, numel(kids));
    assert(ymin >= 3, '%s has a control too low: y=%.1f', names{t}, ymin);
    assert(ymax <= 240, '%s has a control too high: %.1f', names{t}, ymax);
end
delete(app);
cd(here);
test_HX711_BPA_offline;
fprintf('TAB LAYOUT TEST PASS\n');
