function export_robotbody_survey_reference(outputFile)
% Export the saved Appendix C design's length and compliance at the starting pose.
% The source grids use the original route and bracket geometry.
here = fileparts(mfilename('fullpath'));
matlabRoot = fullfile(fileparts(here), 'Matlab');
if nargin < 1
    outputFile = fullfile(here, 'robotbody_survey_reference.csv');
end
addpath(fullfile(matlabRoot, 'Robot_Data'), fullfile(matlabRoot, 'Functions'), ...
    fullfile(matlabRoot, 'Functions','ModernRobotics'), fullfile(matlabRoot,'Mesh_Optimization'));
results = fullfile(matlabRoot, 'Mesh_Optimization','Results');
rows = readtable(fullfile(results,'KompareXXVII','route_configurations.csv'), 'TextType','string');
fid = fopen(outputFile, 'w');
closeFile = onCleanup(@() fclose(fid));
fprintf(fid, 'site,design,output_route_m,input_route_m,bracket_single_m_N,shared_tendon_N_m\n');
for i = 1:height(rows)
    row = rows(i,:);
    [~,name,ext] = fileparts(row.mat_file);
    s = load(fullfile(results,'BiPulley_surveys',name+ext),'g','spec');
    o = s.g.obj;
    assert(o.Diameter == row.diameter_saved_mm && abs(o.RestingL-row.resting_length_m)<1e-12);
    rigidLength = o.MuscleLength;
    [n1,n2] = size(rigidLength);
    bracket = zeros(n1,n2);
    ranges = s.spec.kin.thetaRangesDeg;
    axis1 = linspace(ranges(1,1),ranges(1,2),n1);
    axis2 = linspace(ranges(2,1),ranges(2,2),n2);
    pose = -30 * double(string(s.spec.kin.jointNames) == "knee");
    inputLength = rigidLength + (1-1/row.G_survey) * o.deltaL(:,:,o.ActiveCrossing);
    % Match the local bracket frame used by the source LokBi2 solver.
    switch o.Diameter
        case 20, origin = [9.48,-36.21,30.86]/1000;
        case 10, origin = [-19,22,27.6]/1000;
        otherwise, origin = [0,0,0];
    end
    delta = o.Location(1,:)-origin;
    z = atan2(delta(2),delta(1));
    rz = [cos(z),-sin(z),0;sin(z),cos(z),0;0,0,1];
    deltaZ = rz'*delta';
    y = atan2(deltaZ(3),deltaZ(1));
    ry = [cos(y),0,sin(y);0,1,0;-sin(y),0,cos(y)];
    rotation = rz*ry';
    extension = o.RestingL-rigidLength+o.TendonL+2*o.FittingLength+o.Xi0;
    force = festo4(o.Diameter,extension/(o.RestingL-o.Kmax),o.Pressure).*o.Fmax;
    for j = 1:numel(bracket)
        direction = o.Funit(j,:);
        % Source RowVecTrans applies this translation to the incoming force.
        u = (force(j)*direction-origin)*rotation;
        if norm(u) > 1e-3
            u = u/norm(u);
        else
            u = direction*rotation;
        end
        bracket(j) = sum(u.^2 ./ [o.Xi1,o.Xi2,o.Xi2]);
    end
    get = @(grid) atPose(axis1,axis2,grid,pose);
    fprintf(fid,'%s,%s,%.17g,%.17g,%.17g,%.17g\n',row.site,row.side, ...
        get(rigidLength),get(inputLength),get(bracket),o.kSpr);
end
end

function y = atPose(a,b,grid,pose)
% Bilinear interpolation of the saved five-by-five survey, including rigid axes.
[a,ia] = unique(a,'stable');
[b,ib] = unique(b,'stable');
g = grid(ia,ib);
if isscalar(a)
    y = interp1(b,g,pose(2),'linear');
elseif isscalar(b)
    y = interp1(a,g,pose(1),'linear');
else
    interpolator = griddedInterpolant({a,b},g,'linear','none');
    y = interpolator(pose(1),pose(2));
end
assert(isfinite(y),'Starting pose lies outside the saved survey grid.');
end
