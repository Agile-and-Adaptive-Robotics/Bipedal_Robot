%dump_FlxPin_group_fig.m — read Ben's FlxPin_group.fig annotation geometry
%so the extensor 6-panel letters can copy his exact placement convention.
f = hgload('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo/FlxPin_group.fig');
fprintf('FIG Units=%s Position=%s\n', f.Units, mat2str(f.Position));
ax = reshape(findall(f, 'Type', 'axes'), 1, []);
for k = 1:numel(ax)
    p = ax(k).Position; op = ax(k).OuterPosition;
    ts = ax(k).Title.String; if iscell(ts), ts = strjoin(ts, '|'); elseif ~ischar(ts), ts = ''; end
    fprintf('AX %d: pos=[%.4f %.4f %.4f %.4f] outer=[%.4f %.4f %.4f %.4f] title=%s\n', ...
        k, p, op, ts);
end
tb = reshape(findall(f, 'Type', 'textbox'), 1, []);
for m = 1:numel(tb)
    s = tb(m).String; if iscell(s), s = strjoin(s, '|'); end
    fprintf('TEXTBOX: pos=[%.4f %.4f %.4f %.4f] str=%s VA=%s HA=%s rot=%g\n', ...
        tb(m).Position, s, tb(m).VerticalAlignment, tb(m).HorizontalAlignment, tb(m).Rotation);
end
tx = reshape(findall(f, 'Type', 'text'), 1, []);
for m = 1:numel(tx)
    s = tx(m).String; if iscell(s), s = strjoin(s, '|'); end
    fprintf('TEXT: units=%s pos=%s str=%s\n', tx(m).Units, mat2str(tx(m).Position, 5), s);
end
disp('DONE_dump');
