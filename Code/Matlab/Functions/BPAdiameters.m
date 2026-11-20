clear; clc; close all
%Based on Martens and Boblan and Chou and Hannaford
%% BPA geometry
L0 = [0.6, 0.5, 0.4, 0.3];    % Resting lengths, m

D0 = [0.010 0.020 0.040].';     % Initial inner diameter, m
H0 = 0.0018;                    % Initial membrane thickness, m
KMAX = [0.17 0.255 0.3]';       % Maximum contraction fraction

%% Initial braid angle
theta0a = 28.6;            % deg, Chou and Hannaford / Martens and Boblan

dtheta0b = [-0.1 -0.073 0].';   % deg, Martens and Boblan correction

theta0 = deg2rad(theta0a + dtheta0b);

%% Limit contraction to physical BPA range
relstrain = linspace(0,1,101);      %relstrain = k/KMAX, solve for k
k = KMAX*relstrain;         %Contraction

%% Instantaneous BPA diameter
% From:
%
%   Lfiber = L0/cos(theta0)
%   n      = L0*tan(theta0)/(pi*D0)
%   L      = L0*(1-k)
%
% Substitution eliminates L0, Lfiber, and n.

Din = D0 .* ...
    sqrt(1 - (1-k).^2 .* cos(theta0).^2) ./ ...
    sin(theta0);

D = Din + 2*H0;

%% Outer BPA radius
R = D/2;

%% Figure for showing effect initial diameter on clearance radius
figure('Color','w','Name','Radius')
hold on
    for i = 1:length(D0)

      pR(i) = plot(k(i,:), R(i,:)*1000, ...
            'LineWidth', 2, ...
            'DisplayName', sprintf('D_0 = %.0f mm',1000*D0(i)));

    end
hold off
xlabel('Contraction, (L_0-L)/L_0')
ylabel('Radius, outer, mm')
title('Predicted BPA Radius vs. Contraction')
legend('Location','best')

%% Figure for showing effect initial diameter on clearance diameter
figure('Color','w','Name','Diameter')
hold on
    for i = 1:length(D0)

       pD(i) = plot(k(i,:), D(i,:)*1000, ...
            'LineWidth', 2, ...
            'DisplayName', sprintf('D_0 = %.0f mm',1000*D0(i)));

    end
hold off
xlabel('Contraction, (L_0-L)/L_0')
ylabel('Diameter, outer, mm')
title('Predicted BPA Diameter vs. Contraction')
legend('Location','best')

%% Figure for showing effect of resting length on diameter.
M = numel(L0);
N = numel(KMAX);

% Minimum BPA length
Lmin = (1-KMAX)*L0; % 3xM

% Fiber geometry
Lfiber = L0./cos(theta0);  % 3xM

n = L0.*tan(theta0)./(pi*D0); %3xM

% Length range from resting length to maximum contraction
% permute k from 3x101 to 3x1x101:
% rows = diameter, columns = length, pages = BPA contraction increment
L = pagemtimes(1-permute(k,[1 3 2]),L0);     % 3x4x101

% Inner diameter from braid geometry
Din = sqrt(Lfiber.^2 - L.^2)./(n*pi);        % 3x4x101

% Approximate outer diameter
D = Din + 2*H0;                              % 3x4x101

figure
hold on
c = lines(M);                         % Color identifies resting length
mk = {'o','s','^'};                   % Marker identifies initial diameter
for ii = 1:N
    for i = 1:M

        % Stagger markers so coincident curves remain visible
        markerIdx = i:12:size(D,3);

        plot(squeeze((L0(i)-L(ii,i,:))./L0(i)), ...
             squeeze(D(ii,i,:))*1000, ...
            'LineWidth',2, ...
            'Color',c(i,:), ...
            'Marker',mk{ii}, ...
            'MarkerIndices',markerIdx, ...
            'MarkerSize',6, ...
            'HandleVisibility','off');

    end
end


xlabel('Contraction, (L_0-L)/L_0')
ylabel('Outer Diameter, mm')
title('Predicted BPA Diameter vs. Contraction')
%% Custom two-column legend

grey = [0.5 0.5 0.5];

% Blank handles for column headings
hDtitle = plot(nan,nan,'LineStyle','none','Marker','none');
hLtitle = plot(nan,nan,'LineStyle','none','Marker','none');

% D0 column: marker shape only, grey
hD1 = plot(nan,nan,'LineStyle','none', ...
    'Marker',mk{1},'MarkerSize',6, ...
    'MarkerEdgeColor',grey);

hD2 = plot(nan,nan,'LineStyle','none', ...
    'Marker',mk{2},'MarkerSize',6, ...
    'MarkerEdgeColor',grey);

hD3 = plot(nan,nan,'LineStyle','none', ...
    'Marker',mk{3},'MarkerSize',6, ...
    'MarkerEdgeColor',grey);

% Blank fourth row under D0
hBlank = plot(nan,nan,'LineStyle','none','Marker','none');

% L0 column: colored line only, no markers
hL1 = plot(nan,nan,'-','Color',c(1,:),'LineWidth',2);
hL2 = plot(nan,nan,'-','Color',c(2,:),'LineWidth',2);
hL3 = plot(nan,nan,'-','Color',c(3,:),'LineWidth',2);
hL4 = plot(nan,nan,'-','Color',c(4,:),'LineWidth',2);

% Interleave entries so legend fills row-by-row:
%       D0          L0
%       10 mm       0.6 m
%       20 mm       0.5 m
%       40 mm       0.4 m
%                   0.3 m

hLeg = [ ...
    hDtitle, ...
    hD1, ...
    hD2, ...
    hD3, ...
    hBlank, ...
    hLtitle, ...
    hL1, ...
    hL2, ...
    hL3, ...
    hL4];

legText = { ...
    '{\bf D_0}', ...
    '10 mm', ...
    '20 mm', ...
    '40 mm', ...
    ' ', ...
    '{\bf L_0}', ...
    '0.6 m', ...
    '0.5 m', ...
    '0.4 m', ...
    '0.3 m'};

lgd = legend(hLeg,legText, ...
    'NumColumns',2, ...
    'Location','best', ...
    'Interpreter','tex');

drawnow

% Get the actual legend graphics object
lgd = findobj(gcf,'Type','Legend');

% Legend border
set(lgd,'Box','on','LineWidth',2.5);

% Get legend position in normalized figure coordinates
set(lgd,'Units','normalized');
p = get(lgd,'Position');      % [x y width height]

% Vertical separator between D0 and L0 columns
xSep = p(1) + 0.5*p(3);

annotation(gcf,'line', ...
    [xSep xSep], ...
    [p(2) p(2)+p(4)], ...
    'LineWidth',1.5);

% Horizontal separator below header row
ySep = p(2) + 4/5*p(4);

annotation(gcf,'line', ...
    [p(1) p(1)+p(3)], ...
    [ySep ySep], ...
    'LineWidth',1.5);