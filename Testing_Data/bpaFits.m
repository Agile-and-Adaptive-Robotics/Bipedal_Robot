function [fitresult, gof, output, valid] = bpaFits(XX40, YY40, ZZ40, XX, YY, ZZ, XX10, YY10, ZZ10, Ax, Ay, Az, XX20, YY20, ZZ20)
%CREATEFITS(XX20,YY20,ZZ20,XX10,YY10,ZZ10,XX,YY,ZZ,XX40,YY40,ZZ40)
%  Create fits for 10mm, 20mm, and 40mm BPAs using experimental and Festo
%  data. Run after running "BenLawrenceData.m" in the Testing_Data folder.
%
%  Data for 'Exponential 20, Festo' fit:
%      X Input : XX20
%      Y Input : YY20
%      Z Output: ZZ20
%      Validation X: XX20
%      Validation Y: YY20
%      Validation Z: ZZ20
%  Data for 'Exp10, Festo vs Experiment' fit:
%      X Input : XX10
%      Y Input : YY10
%      Z Output: ZZ10
%      Validation X: XX
%      Validation Y: YY
%      Validation Z: ZZ
%  Data for 'Exp10, Experiment vs Festo' fit:
%      X Input : XX
%      Y Input : YY
%      Z Output: ZZ
%      Validation X: XX10
%      Validation Y: YY10
%      Validation Z: ZZ10
%  Data for 'Exponential 40 normalized' fit:
%      X Input : XX40
%      Y Input : YY40
%      Z Output: ZZ40
%      Validation X: XX40
%      Validation Y: YY40
%      Validation Z: ZZ40
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 17-Dec-2023 17:29:42

%% Initialization.

%Use these commands if variables are not loaded into workspace
% load FestoData.mat XX40 YY40 ZZ40 XX10 YY10 ZZ10 XX20 YY20 ZZ20
% load allData.mat XX YY ZZ
% load data20mm.mat Ax Ay Az

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell( 3, 1 );
gof = struct( 'sse', cell( 3, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
output = cell( 3,1);

%% Fit: 'Exp10, Experiment vs Festo'.
[xData, yData, zData] = prepareSurfaceData( XX, YY, ZZ );

% Set up fittype and options.
ft = fittype( 'a0*(exp(-a1.*x)-1)+y*exp(-a3*((x).^2))', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0];
opts.MaxFunEvals = 6000;
opts.MaxIter = 2000;
opts.Robust = 'LAR';
opts.StartPoint = [0.5822 4.142 0.5368];

% Fit model to data.
[fitresult{1}, gof(1), output{1}] = fit( [xData, yData], zData, ft, opts );

% Compare against validation data.
[xValidation, yValidation, zValidation] = prepareSurfaceData( XX10, YY10, ZZ10 );
residual = zValidation - fitresult{1}( xValidation, yValidation );
nNaN = nnz( isnan( residual ) );
residual(isnan( residual )) = [];
sse = norm( residual )^2;
rmse = sqrt( sse/length( residual ) );
sst = norm(zValidation - mean(zValidation))^2;
Rsquare = 1-(sse*(length(residual)-1))/(sst*(length(residual)));
fprintf( 'Goodness-of-validation for ''%s'' fit:\n', 'Exp10, Experiment vs Festo' );
fprintf( '    SSE : %f\n', sse );
fprintf( '    RMSE : %f\n', rmse );
fprintf( '    Adj. R^2 : %f\n', Rsquare );
fprintf( '    Largest residual: %f\n',max(residual));
fprintf( '    %i points outside domain of data.\n', nNaN );

valid(1,:) = table(sse, rmse, sst, Rsquare);
% Create a figure for the plots.
figure( 'Name', 'Exp10, Experiment vs Festo' );

% Compute limits for axes.
xlim = [min( [xData; xValidation] ), max( [xData; xValidation] )];
ylim = [min( [yData; yValidation] ), max( [yData; yValidation] )];

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult{1}, [xData, yData], zData, 'XLim', xlim, 'YLim', ylim );
% Add validation data to plot.
hold on
h(end+1) = plot3( xValidation, yValidation, zValidation, 'bo', 'MarkerFaceColor', 'w' );
hold off
legend( h, 'Exp10, Experiment vs Festo', 'ZZ vs. XX, YY', 'ZZ10 vs. XX10, YY10', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'XX', 'Interpreter', 'none' );
ylabel( 'YY', 'Interpreter', 'none' );
zlabel( 'ZZ', 'Interpreter', 'none' );
grid on
view( 90.0, 0.0 );

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult{1}, [xData, yData], zData, 'Style', 'Residual', 'XLim', xlim, 'YLim', ylim );
% Add validation data to plot.
hold on
h(end+1) = plot3( xValidation, yValidation, zValidation - fitresult{1}( xValidation, yValidation ), 'bo', 'MarkerFaceColor', 'w' );
hold off
legend( h, 'Exp10, Experiment vs Festo - residuals', 'Exp10, Experiment vs Festo - validation residuals', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'XX', 'Interpreter', 'none' );
ylabel( 'YY', 'Interpreter', 'none' );
zlabel( 'ZZ', 'Interpreter', 'none' );
grid on
view( 90.0, 0.0 );

%% Fit: 'Exponential 20, Festo'.
[xData, yData, zData] = prepareSurfaceData(  Ax, Ay, Az);

% Set up fittype and options.
ft = fittype( 'a0*(exp(-a1.*x)-1)+y*exp(-a3*((x).^2))', 'independent', {'x', 'y'}, 'dependent', 'z' );
excludedPoints = zData < 0;
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0];
opts.MaxFunEvals = 6000;
opts.MaxIter = 4000;
opts.Robust = 'LAR';
opts.StartPoint = [0.2607 6.398 1.303];
opts.TolFun = 1e-07;
opts.TolX = 1e-07;
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult{2}, gof(2), output{2}] = fit( [xData, yData], zData, ft, opts );

% Compare against validation data.
[xValidation, yValidation, zValidation] = prepareSurfaceData( XX20, YY20, ZZ20 );
residual = zValidation - fitresult{2}( xValidation, yValidation );
nNaN = nnz( isnan( residual ) );
residual(isnan( residual )) = [];
sse = norm( residual )^2;
rmse = sqrt( sse/length( residual ) );
sst = norm(zValidation - mean(zValidation))^2;
Rsquare = 1-(sse*(length(residual)-1))/(sst*(length(residual)));
fprintf( 'Goodness-of-validation for ''%s'' fit:\n', 'Exponential 20, Festo' );
fprintf( '    SSE : %f\n', sse );
fprintf( '    RMSE : %f\n', rmse );
fprintf( '    Adj. R^2 : %f\n', Rsquare );
fprintf( '    Largest residual: %f\n',max(residual));
fprintf( '    %i points outside domain of data.\n', nNaN );

valid(2,:) = table(sse, rmse, sst, Rsquare);

% Create a figure for the plots.
figure( 'Name', 'Exponential 20, Festo' );

% Compute limits for axes.
xlim = [min( [xData; xValidation] ), max( [xData; xValidation] )];
ylim = [min( [yData; yValidation] ), max( [yData; yValidation] )];

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult{2}, [xData, yData], zData, 'Exclude', excludedPoints, 'XLim', xlim, 'YLim', ylim );
% Add validation data to plot.
hold on
h(end+1) = plot3( xValidation, yValidation, zValidation, 'bo', 'MarkerFaceColor', 'w' );
hold off
legend( h, 'Exp20, Data vs Festo', 'Az vs. Ax, Ay', 'Excluded Az vs. Ax, Ay', 'ZZ20 vs. XX20, YY20', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Ax', 'Interpreter', 'none' );
ylabel( 'Ay', 'Interpreter', 'none' );
zlabel( 'Az', 'Interpreter', 'none' );
grid on
view( 0.0, 0.0 );

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult{2}, [xData, yData], zData, 'Style', 'Residual', 'Exclude', excludedPoints, 'XLim', xlim, 'YLim', ylim );
% Add validation data to plot.
hold on
h(end+1) = plot3( xValidation, yValidation, zValidation - fitresult{2}( xValidation, yValidation ), 'bo', 'MarkerFaceColor', 'w' );
hold off
legend( h, 'Exp20, Data vs Festo - residuals', 'Excluded Az vs. Ax, Ay', 'Exp20, Data vs Festo - validation residuals', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Ax', 'Interpreter', 'none' );
ylabel( 'Ay', 'Interpreter', 'none' );
zlabel( 'Az', 'Interpreter', 'none' );
grid on
view( 0.0, 0.0 );


%% Fit: 'Exponential 40 normalized'.
[xData, yData, zData] = prepareSurfaceData( XX40, YY40, ZZ40 );

% Set up fittype and options.
ft = fittype( 'a0*(exp(-a1.*x)-1)+y*exp(-a3*((x).^2))', 'independent', {'x', 'y'}, 'dependent', 'z' );
excludedPoints = zData < 0;
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 -Inf 0];
opts.MaxFunEvals = 7000;
opts.MaxIter = 5000;
opts.Robust = 'LAR';
opts.StartPoint = [0.1224 10.47 2.023];
opts.Exclude = excludedPoints;


% Fit model to data.
[fitresult{3}, gof(3), output{3}] = fit( [xData, yData], zData, ft, opts );

% Compare against validation data.
[xValidation, yValidation, zValidation] = prepareSurfaceData( XX40, YY40, ZZ40 );
residual = zValidation - fitresult{3}( xValidation, yValidation );
nNaN = nnz( isnan( residual ) );
residual(isnan( residual )) = [];
sse = norm( residual )^2;
rmse = sqrt( sse/length( residual ) );
sst = norm(zValidation - mean(zValidation))^2;
Rsquare = 1-(sse*(length(residual)-1))/(sst*(length(residual)));
fprintf( 'Goodness-of-validation for ''%s'' fit:\n', 'Exponential 40 normalized' );
fprintf( '    SSE : %f\n', sse );
fprintf( '    RMSE : %f\n', rmse );
fprintf( '    Adj. R^2 : %f\n', Rsquare );
fprintf( '    Largest residual: %f\n',max(residual));
fprintf( '    %i points outside domain of data.\n', nNaN );

valid(3,:) = table(sse, rmse, sst, Rsquare);
% Create a figure for the plots.
figure( 'Name', 'Exponential 40 normalized' );

% Compute limits for axes.
xlim = [min( [xData; xValidation] ), max( [xData; xValidation] )];
ylim = [min( [yData; yValidation] ), max( [yData; yValidation] )];

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult{3}, [xData, yData], zData, 'XLim', xlim, 'YLim', ylim );
% Add validation data to plot.
hold on
h(end+1) = plot3( xValidation, yValidation, zValidation, 'bo', 'MarkerFaceColor', 'w' );
hold off
legend( h, 'Exponential 40 normalized', 'ZZ40 vs. XX40, YY40', 'ZZ40 vs. XX40, YY40', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'XX40', 'Interpreter', 'none' );
ylabel( 'YY40', 'Interpreter', 'none' );
zlabel( 'ZZ40', 'Interpreter', 'none' );
grid on
view( 0.0, 0.0 );

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult{3}, [xData, yData], zData, 'Style', 'Residual', 'XLim', xlim, 'YLim', ylim );
% Add validation data to plot.
hold on
h(end+1) = plot3( xValidation, yValidation, zValidation - fitresult{3}( xValidation, yValidation ), 'bo', 'MarkerFaceColor', 'w' );
hold off
legend( h, 'Exponential 40 normalized - residuals', 'Exponential 40 normalized - validation residuals', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'XX40', 'Interpreter', 'none' );
ylabel( 'YY40', 'Interpreter', 'none' );
zlabel( 'ZZ40', 'Interpreter', 'none' );
grid on
view( 0.0, 0.0 );


%% save it
save bpaFitsResult.mat fitresult gof output valid XX40 YY40 ZZ40 XX YY ZZ XX10 YY10 ZZ10 Ax Ay Az XX20 YY20 ZZ20