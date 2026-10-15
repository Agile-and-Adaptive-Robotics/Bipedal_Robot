% Verify HX711_BPA shrink-to-fit: window on-screen, components inside.
here = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Code/Matlab/HX711 v3.0/HX711 v3.0';
addpath(here);

app = HX711_BPA(here, true);   % keepHidden
fig = app.MatlabArduinoHX711UIFigure;
ss = get(groot,'ScreenSize');
p = fig.Position;
fprintf('screen: [%d %d %d %d]\n', ss);
fprintf('fig outer: [%d %d %d %d]\n', p);
assert(p(1) >= 1 && p(2) >= 1, 'window off-screen at bottom/left');
assert(p(1)+p(3) <= ss(3)+1, 'window overflows right');
assert(p(2)+p(4) <= ss(4)+1, 'window overflows top');
assert(p(3) <= ss(3) && p(4) <= ss(4), 'window larger than screen');

ip = fig.InnerPosition;
allH = gatherKids(fig.Children);
bad = 0; minFont = inf;
for k = 1:numel(allH)
    h = allH(k);
    q = h.Position;
    % -30 y-tolerance: the original design hangs the GlobalSettings
    % TabGroup 25 px below its panel edge on purpose (tab strip flush
    % with the panel bottom); scaled it becomes -22.
    if q(1) < -1 || q(2) < -30 || q(1)+q(3) > ip(3)+2 || q(2)+q(4) > ip(4)+2
        bad = bad + 1;
        fprintf('OUT OF CLIENT: %s at [%.0f %.0f %.0f %.0f] (client %.0fx%.0f)\n', ...
            class(h), q, ip(3), ip(4));
    end
    if isprop(h,'FontSize') && isnumeric(h.FontSize)
        minFont = min(minFont, h.FontSize);
    end
end
fprintf('components checked: %d, out-of-client: %d, min font: %g, client: [%.0f %.0f %.0f %.0f]\n', ...
    numel(allH), bad, minFont, ip);
assert(bad == 0, 'components outside client area');
assert(minFont >= 8, 'font below floor');
delete(app);
fprintf('FIT TEST PASS\n');

    function kids = gatherKids(hs)
        kids = hs(:).';
        for k = 1:numel(hs)
            if isprop(hs(k),'Children') && ~isempty(hs(k).Children)
                kids = [kids, gatherKids(hs(k).Children)];
            end
        end
    end
