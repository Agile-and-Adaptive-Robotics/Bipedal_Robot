function geo = buildGeoExclusion()
%BUILDGEOEXCLUSION Geometry for p2 and inflated-BPA exclusion.
%
% p2 = [px py pz]
%
% Case 1:
%   py >= -0.00887
%   use outer X-Z profile from SolidWorks sketch
%
% Case 2:
%   py < -0.00887
%   use inner circle in X-Z plane:
%       diameter = 0.0381 m
%       center   = [0.01586, 0, 0]

geo.clearance = 0.001;

% Collision envelope for the nominal 20 mm BPA.  The measured/selected
% maximum inflated diameter is 38.5 mm.  The tendon is 2.0 mm diameter.
% Collision is checked at +5 deg, the physical robot extension limit,
% where the flexor has the least joint clearance.  The centerline sampling
% allowance in nonlconExclusion makes the discrete test conservative.
geo.inflatedDiameter = 0.0385;
geo.bpaRadius = geo.inflatedDiameter/2;
geo.tendonDiameter = 0.002;
geo.tendonRadius = geo.tendonDiameter/2;
geo.collisionAngleD = 5;
geo.centerlineSampleStep = 0.001;

geo.yCut = -0.00887;

% Case 2 circle
geo.circleCenter = [0.01586, 0, 0];
geo.circleRadius = 0.0381/2;

% Case 1 outer profile
geo.xBack = -0.02064;
geo.zBackTop = -0.0025;
geo.zBackBot =  0.0025;
geo.R = 0.035;

% Horizontal bottom is 37.5 mm below z = 0 in the SolidWorks top view.
geo.zBottom = 0.0375;

% Femoral-condyle ellipse from the extensor geometry, expressed in the
% femur-frame X-Y plane.  Keep the physical ellipse here; the inflated BPA
% radius is applied as a true Euclidean distance offset by the constraint.
geo.femurProfileCenter = [0.01817, -0.41031];
geo.femurEllipsePhysicalA = 0.06640/2;
geo.femurEllipsePhysicalB = 0.05439/2;
geo.femurEllipseTheta = deg2rad(25.11);

% The condyles have the same Z width as the proximal-tibia outer profile.
geo.femurZMin = geo.zBackTop;
geo.femurZMax = geo.zBottom;

ellipseParameter = linspace(0, 2*pi, 721).';
ellipseLocal = [geo.femurEllipsePhysicalA*cos(ellipseParameter), ...
                geo.femurEllipsePhysicalB*sin(ellipseParameter)];
Rell = [cos(geo.femurEllipseTheta), -sin(geo.femurEllipseTheta); ...
        sin(geo.femurEllipseTheta),  cos(geo.femurEllipseTheta)];
ellipseFemur = ellipseLocal*Rell.' + geo.femurProfileCenter;
geo.femurEllipseX = ellipseFemur(:,1);
geo.femurEllipseY = ellipseFemur(:,2);

% Pick this as far as you need the horizontal segment to exist.
% Use your p2 upper x bound here if that is larger.
geo.xAfter = 0.0375;

xc = geo.xBack + geo.R;
zc = geo.zBackBot;

% Back vertical edge
zEdge = linspace(geo.zBackTop, geo.zBackBot, 20).';
xEdge = geo.xBack*ones(size(zEdge));

% 90-degree radius corner.
% Start at phi = pi:
%   x = xBack, z = zBackBot
% End at phi = pi/2:
%   x = xc, z = zBackBot + R = 0.0375
phi = linspace(pi, pi/2, 100).';
xArc = xc + geo.R*cos(phi);
zArc = zc + geo.R*sin(phi);

% Bottom horizontal continuation in positive x
xBot = linspace(xc, geo.xAfter, 40).';
zBot = geo.zBottom*ones(size(xBot));

xProf = [xEdge; xArc(2:end); xBot(2:end)];
zProf = [zEdge; zArc(2:end); zBot(2:end)];

% Convert to polar, using 0 to 2*pi so the small negative-z back-edge
% does not wrap across -pi/pi.
thetaProf = atan2(zProf, xProf);
thetaProf(thetaProf < 0) = thetaProf(thetaProf < 0) + 2*pi;

rProf = hypot(xProf, zProf);

[thetaProf, I] = sort(thetaProf);
rProf = rProf(I);
xProf = xProf(I);
zProf = zProf(I);

[thetaProf, I] = unique(thetaProf, 'stable');
rProf = rProf(I);
xProf = xProf(I);
zProf = zProf(I);

geo.thetaVec = thetaProf(:).';
geo.rVec = rProf(:).';

geo.xProf = xProf;
geo.zProf = zProf;

% -------------------------------------------------------------------------
% DEBUG PLOTS
% Set doDebug = true to sanity-check Cartesian profile vs polar lookup.
% -------------------------------------------------------------------------
doDebug = false;

if doDebug
    figure; hold on; grid on; axis equal;

    plot(geo.xProf, geo.zProf, 'o-', 'LineWidth', 1.5);

    thetaFine = linspace(min(geo.thetaVec), max(geo.thetaVec), 500);
    rFine = interp1(geo.thetaVec, geo.rVec, thetaFine, 'pchip');

    xFine = rFine.*cos(thetaFine);
    zFine = rFine.*sin(thetaFine);

    plot(xFine, zFine, '-', 'LineWidth', 1.5);

    th = linspace(0, 2*pi, 300);
    xCircle = geo.circleCenter(1) + geo.circleRadius*cos(th);
    zCircle = geo.circleCenter(3) + geo.circleRadius*sin(th);

    plot(xCircle, zCircle, '--', 'LineWidth', 1.5);

    xlabel('x [m]');
    ylabel('z [m]');
    title('X-Z Exclusion Boundary Sanity Check');
    legend('Outer Cartesian profile', ...
           'Outer polar lookup converted to Cartesian', ...
           'Inner circle for py < -0.00887', ...
           'Location', 'best');

    % Make the MATLAB plot look like the SolidWorks top-plane view,
    % where positive z points downward.
    set(gca, 'YDir', 'reverse');

    figure; hold on; grid on; axis equal;

    yOuter = 0.010*ones(size(xFine));
    plot3(xFine, yOuter, zFine, '-', 'LineWidth', 1.5);

    yCircle = geo.yCut*ones(size(xCircle));
    plot3(xCircle, yCircle, zCircle, '--', 'LineWidth', 1.5);

    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    title('3D Case-Switch Exclusion Sanity Check');
    legend('Outer profile case: py >= -0.00887', ...
           'Inner circle case: py < -0.00887', ...
           'Location', 'best');
    view(3);
end

end
