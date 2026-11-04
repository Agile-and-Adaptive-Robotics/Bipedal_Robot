function route = fr_muscleroute(osimPath, muscleName)
%FR_MUSCLEROUTE Parse one muscle's GeometryPath straight from a .osim file.
%
% Returns the declared path points exactly as the XML carries them:
%   .points(i): .name .body .loc (1x3 m) .type ('PathPoint' |
%     'ConditionalPathPoint' | 'MovingPathPoint')
%     conditional points add .coord (coordinate NAME) and .range [min max] rad
%     moving points add .coordXYZ (coordinate NAME per component) and
%       .fnXYZ (struct with .x .y .y nodes for the SimmSpline per component)
% The declared behaviour is implemented at evaluation time in fr_evalRoute;
% nothing is frozen at the default pose (master directive item 1).
doc = xmlread(osimPath);
root = doc.getDocumentElement();
% gait2392 muscles are Thelen2003Muscle elements (a plain 'Muscle'
% tag does not exist in the model)
muscles = root.getElementsByTagName('Thelen2003Muscle');
target = [];
for i = 0:muscles.getLength()-1
    m = muscles.item(i);
    if strcmp(char(m.getAttribute('name')), muscleName)
        target = m;
        break
    end
end
if isempty(target)
    error('fr_muscleroute:NoMuscle', 'muscle %s not found in %s', ...
        muscleName, osimPath);
end
gp = target.getElementsByTagName('GeometryPath');
pps = gp.item(0).getElementsByTagName('PathPointSet');
objs = pps.item(0).getElementsByTagName('objects');
kids = objs.item(0).getChildNodes();
route = struct();
route.muscle = muscleName;
route.points = {};   % NOT struct('points', {}): a cell value would make a 0x0 struct array
for i = 0:kids.getLength()-1
    c = kids.item(i);
    if c.getNodeType() ~= c.ELEMENT_NODE; continue; end
    P = struct();
    P.name = char(c.getAttribute('name'));
    P.type = char(c.getNodeName());
    P.body = fr_stripFrame(fr_xmlChild(c, 'socket_parent_frame'));
    loc = fr_xmlChild(c, 'location');
    if ~isempty(loc)
        P.loc = fr_xmlNums(loc);
    else
        P.loc = [NaN NaN NaN];   % moving point: location comes from functions
    end
    sc = fr_xmlChild(c, 'socket_coordinate');
    if ~isempty(sc)   % conditional
        P.coord = fr_stripFrame(sc);
        P.range = fr_xmlNums(fr_xmlChild(c, 'range'));
    end
    sx = fr_xmlChild(c, 'socket_x_coordinate');
    if ~isempty(sx)   % moving
        P.coordXYZ = {fr_stripFrame(sx), ...
            fr_stripFrame(fr_xmlChild(c, 'socket_y_coordinate')), ...
            fr_stripFrame(fr_xmlChild(c, 'socket_z_coordinate'))};
        P.fnXYZ = struct();
        P.fnXYZ.x = fr_parseFn(fr_xmlChild(c, 'x_location'));
        P.fnXYZ.y = fr_parseFn(fr_xmlChild(c, 'y_location'));
        P.fnXYZ.z = fr_parseFn(fr_xmlChild(c, 'z_location'));
    end
    route.points{end+1} = P;
end
end

function F = fr_parseFn(node)
F = struct('type', 'constant', 'value', 0);
if isempty(node); return; end
fnl = fr_xmlChild(node, 'LinearFunction');
fsp = fr_xmlChild(node, 'SimmSpline');
fcn = fr_xmlChild(node, 'Constant');
if ~isempty(fnl)
    F.type = 'linear'; F.coeffs = fr_xmlNums(fnl);
elseif ~isempty(fsp)
    F.type = 'spline';
    F.x = fr_xmlNums(fr_xmlChild(fsp, 'x'));
    F.y = fr_xmlNums(fr_xmlChild(fsp, 'y'));
elseif ~isempty(fcn)
    F.type = 'constant'; F.value = fr_xmlNum(fcn);
end
end

function s = fr_stripFrame(n)
s = strrep(strrep(strtrim(char(n.getTextContent())), '/bodyset/', ''), ...
    '/jointset/', '');
s = strrep(s, '/ground', 'ground');
% coordinates are addressed as /jointset/<joint>/<coord>: keep the LAST
% path segment (the bare coordinate name used as the q field)
sl = strfind(s, '/');
if ~isempty(sl)
    s = s(sl(end) + 1:end);
end
end
