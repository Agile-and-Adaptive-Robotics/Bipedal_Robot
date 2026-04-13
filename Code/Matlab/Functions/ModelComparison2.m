function Goof = ModelComparison2(D, Lrest, Lmin, Data, varargin)
%% Compare Martens and Sarosi models to ours
% Input convention:
%   D     = {10; 20};
%   Lrest = {[0.257 0.233]; [0.300 0.450 0.509]};
%   Data  = {{W1, W2, []}; {M30, M45, M50}};
%
% Each measured-data entry must be an N-by-3 matrix:
%   col 1 -> relative contraction coordinate x
%   col 2 -> pressure [kPa]
%   col 3 -> force [N]
%
% Empty entries [] are allowed and skipped.
%
% Notes:
%   - Your x is your relative contraction coordinate
%   - Sarosi uses x*kmax in your convention
%   - Sarosi pressure is passed in bar
%   - Our model pressure is normalized by 620 kPa
%   - Martens pressure is passed in Pa
%
% Output:
%   Goof is a struct containing cases, tables, and figure handles

%% Options
p = inputParser;
addParameter(p,'SarosiBest','off');
addParameter(p,'MartensBest','off');
parse(p,varargin{:});

useSarosiBest  = strcmpi(p.Results.SarosiBest,'on');
useMartensBest = strcmpi(p.Results.MartensBest,'on');

MartensDOffset_mm = 0;
if useMartensBest
    MartensDOffset_mm = 1.8;
end

%% Load our lookup models
load FestoLookup.mat f_10 f20

% Set global options for figures
set(groot, ...
    'defaultAxesLineWidth', 3, ...
    'defaultAxesFontSize', 12, ...
    'defaultAxesFontWeight', 'bold')

ms = 60;

%% Colors that are accessible
% c = cell(7,1);
% c{1} = '#FFD700'; % gold
% c{2} = '#FFB14E'; % orange
% c{3} = '#FA8775'; % light orange
% c{4} = '#EA5F94'; % pink
% c{5} = '#CD34B5'; % magenta
% c{6} = '#9D02D7'; % magenta 2
% c{7} = '#0000FF'; % indigo
% 
% d = hex2rgb(c); 


% Output will be 
d =[1.0000    0.8431         0;
    1.0000    0.6941    0.3059;
    0.9804    0.5294    0.4588;
    0.9176    0.3725    0.5804;
    0.8039    0.2039    0.7098;
    0.6157    0.0078    0.8431;
         0         0    1.0000];

% Assign meaning 
colExp     = d(7,:); 
colBolen   = d(6,:); 
colSarosiA = d(3,:); 
colSarosiB = d(2,:); 
colMartens = d(4,:); 

cmap = interp1(1:7, d, linspace(1,7,256), 'linear');
colormap(cmap)

%% Coefficients
a10  = [-9.2194029, 203.7012413, -0.34221042, -3.2255991, 109.2038216, -208.372034];
a20a = [-4.00180705, 292.4620246, -0.32930845, -9.33564098, 294.0538256, -280.498151];
a20b = [-4.35573, 281.22370, -0.32866, -9.27035, 302.20107, -263.69185];

% Converted to SI units from reported coefficients
c10 = [74.085e6, -689.20e6, 1.837e9,  -848.79e6, -0.1];
c20 = [93.232e6, -715.29e6, 1.5483e9, -502.95e6, -0.073];

%% Published data from Martens paper
% Appendix A1: DMSP-20-300 measured force map
L_A1 = [0.296 0.282 0.267 0.250 0.237 0.222];   % [m]

P_A1 = [ ...
     35000
     65000
     85000
    110000
    135000
    146250
    157500
    168750
    180000
    197500
    215000
    232500
    250000
    281250
    312500
    343750
    375000
    450000
    525000
    600000
    675000];                                    % [Pa]

F_A1 = [ ...
      0      0      0      0      0      0
    54.9     0      0      0      0      0
   113.4     0      0      0      0      0
   179.0     0      0      0      0      0
   241.8     0      0      0      0      0
   273.3   21.1     0      0      0      0
   300.2   42.9     0      0      0      0
   331.0   70.3     0      0      0      0
   355.6   93.7     0      0      0      0
   402.5  131.0   27.7     0      0      0
   446.5  165.2   61.6     0      0      0
   495.9  204.1   92.7     0      0      0
   543.4  244.7  122.0     0      0      0
   623.5  309.7  170.2   42.8     0      0
   708.8  380.4  221.0   80.6     0      0
   787.9  443.7  277.4  118.7     0      0
   870.9  518.9  328.0  157.6     0      0
  1065.2  676.4  462.3  252.1   65.5     0
  1268.9  837.5  587.8  352.3  124.7     0
  1460.1 1004.0  714.8  443.5  188.1     0
  1653.5 1164.3  840.7  535.8  253.3     0];

% Appendix A2: DMSP-10-250 measured force map
L_A2 = [0.250 0.234 0.219 0.204 0.198];   % [m]

P_A2 = [ ...
     50000
     68750
    132500
    196250
    260000
    285000
    310000
    335000
    360000
    402187
    441375
    486562
    528750
    559062
    589375
    619688
    650000];   % [Pa]

F_A2 = [ ...
      0        0        0        0       0
    32.1667    0        0        0       0
    84.8333    0        0        0       0
   139.5       0        0        0       0
   195.6667    0        0        0       0
   217.6667   17        0        0       0
   239.6667   35.1667   0        0       0
   261.6667   52.6667   0        0       0
   283.8333   71.3333   0        0       0
   321.5     101.5     22.8333   0       0
   359.0     132.3333  46.1667   0       0
   396.6667  162.5     69.3333   0       0
   434.3333  193.6667  92.1667   0       0
   461.6667  216.3333 109.3333  11.5     0
   488.3333  238.1667 126.3333  23.3333  0
   515.8333  261.0    142.6667  35.0     0
   542.5     282.8333 158.6667  46.0     0];

%% Build cases from arbitrary input convention
cases = struct([]);
kk = 0;

for iD = 1:numel(D)
    Dia_mm = D{iD};
    Lrest_i = Lrest{iD};
    Lmin_i = Lmin{iD};
    Data_i = Data{iD};

    if numel(Lrest_i) ~= numel(Data_i) || numel(Lrest_i) ~= numel(Lmin_i)
        error('For D{%d}, resting lengths, Lmin values, and data cells must have matching counts.', iD);
    end

    for iL = 1:numel(Lrest_i)
        thisData = Data_i{iL};

        if isempty(thisData)
            continue
        end

        if size(thisData,2) ~= 3
            error('Data{%d}{%d} must be N-by-3 with columns [x p z].', iD, iL);
        end

        kk = kk + 1;

        cases(kk).Diameter_mm = Dia_mm;
        cases(kk).RestLength_m = Lrest_i(iL);
        cases(kk).RestLength_mm = 1000 * Lrest_i(iL);
        cases(kk).Lmin_m = Lmin_i(iL);
        cases(kk).Lmin_mm = 1000 * Lmin_i(iL);

        cases(kk).x = thisData(:,1);
        cases(kk).p_kPa = 620*thisData(:,2);    %Convert normalized pressure to kPa
        cases(kk).z = thisData(:,3);

        cases(kk).CaseIndexWithinDiameter = iL;
        cases(kk).CaseLabel = sprintf('%d mm, %d mm', Dia_mm, round(1000 * Lrest_i(iL)));
        cases(kk).XAxisLabel = 'Relative contraction';
        cases(kk).YAxisLabel = 'Pressure (kPa)';
        cases(kk).ZAxisLabel = 'Force (N)';
        cases(kk).FigureTitle = cases(kk).CaseLabel;

        if Dia_mm == 10
            cases(kk).OurModelFun = f_10;
            cases(kk).MartensCoeff = c10;
            cases(kk).MartensDiameter_mm = 10.0 + MartensDOffset_mm;
            cases(kk).SarosiCoeffPrimary = a10;
            cases(kk).SarosiCoeffAlt = [];
            cases(kk).SarosiPrimaryLabel = 'Sarosi';
            cases(kk).SarosiAltLabel = '';
            cases(kk).HasSarosiComparison = false;
            cases(kk).Models = {'Experimental data','Bolen prediction','Sarosi prediction','MB prediction'};
        elseif Dia_mm == 20
            cases(kk).OurModelFun = f20;
            cases(kk).MartensCoeff = c20;
            cases(kk).MartensDiameter_mm = 20.0 + MartensDOffset_mm;

            if abs(Lrest_i(iL) - 0.300) < 1e-12
                cases(kk).SarosiCoeffPrimary = a20a;
                cases(kk).SarosiCoeffAlt = a20b;
                cases(kk).SarosiPrimaryLabel = 'DMSP-20-200N';
                cases(kk).SarosiAltLabel = 'DMSP-20-400N';
                cases(kk).HasSarosiComparison = true;
                cases(kk).Models = {'Experimental data','Bolen prediction','Sarosi_A prediction','Sarosi_B prediction','MB prediction'};
            else
                cases(kk).SarosiCoeffPrimary = a20b;
                cases(kk).SarosiCoeffAlt = [];
                cases(kk).SarosiPrimaryLabel = 'DMSP-20-400N';
                cases(kk).SarosiAltLabel = '';
                cases(kk).HasSarosiComparison = false;
                cases(kk).Models = {'Experimental data','Bolen prediction','Sarosi_B prediction','MB prediction'};
            end
        else
            error('Unsupported diameter: %g mm', Dia_mm);
        end
    end
end

nCase = numel(cases);

if nCase == 0
    error('No nonempty measured-data cases were provided.');
end

%% Evaluate models for each case
for i = 1:nCase
    x = cases(i).x(:);      %relative strain
    p_raw = cases(i).p_kPa(:)/620;
    % z = cases(i).z(:);      %Experimental data
    
    p_norm = p_raw;         %of course this is a little inefficient but it works
    p_bar  = p_raw * 6.2;
    p_Pa   = p_raw * 6.2E5;

    cases(i).OurModel = cases(i).OurModelFun(x, p_norm) .* ...
        maxBPAforce(cases(i).RestLength_m, num2str(cases(i).Diameter_mm), 620);

    [cases(i).Martens, cases(i).Lmin_used] = Martens( ...
        cases(i).RestLength_m, ...
        cases(i).Lmin_m, ...
        cases(i).MartensDiameter_mm, ...
        cases(i).MartensCoeff, ...
        x, p_Pa);

        [cases(i).SarosiPrimary, cases(i).SarosiPrimary_kmax] = sar( ...
        cases(i).RestLength_m, ...
        cases(i).Diameter_mm, ...
        cases(i).SarosiCoeffPrimary, ...
        x, p_bar);

    if cases(i).HasSarosiComparison
        [cases(i).SarosiAlt, cases(i).SarosiAlt_kmax] = sar( ...
            cases(i).RestLength_m, ...
            cases(i).Diameter_mm, ...
            cases(i).SarosiCoeffAlt, ...
            x, p_bar);
    else
        cases(i).SarosiAlt = [];
        cases(i).SarosiAlt_kmax = [];
    end

    [cases(i).OurModel_RMSE, cases(i).OurModel_FVU, cases(i).OurModel_MaxResidual] = ...
        Go_OfF(cases(i).z, cases(i).OurModel(:));

    [cases(i).Martens_RMSE, cases(i).Martens_FVU, cases(i).Martens_MaxResidual] = ...
        Go_OfF(cases(i).z, cases(i).Martens(:));

    [cases(i).SarosiPrimary_RMSE, cases(i).SarosiPrimary_FVU, cases(i).SarosiPrimary_MaxResidual] = ...
        Go_OfF(cases(i).z, cases(i).SarosiPrimary(:));

    if cases(i).HasSarosiComparison
        [cases(i).SarosiAlt_RMSE, cases(i).SarosiAlt_FVU, cases(i).SarosiAlt_MaxResidual] = ...
            Go_OfF(cases(i).z, cases(i).SarosiAlt(:));
    else
        cases(i).SarosiAlt_RMSE = [];
        cases(i).SarosiAlt_FVU = [];
        cases(i).SarosiAlt_MaxResidual = [];
    end

    if useSarosiBest && cases(i).HasSarosiComparison
        if cases(i).SarosiPrimary_RMSE <= cases(i).SarosiAlt_RMSE
            cases(i).SarosiBest = cases(i).SarosiPrimary;
            cases(i).SarosiBest_kmax = cases(i).SarosiPrimary_kmax;
            cases(i).SarosiBestLabel = 'Sarosi';
            cases(i).SarosiBest_RMSE = cases(i).SarosiPrimary_RMSE;
            cases(i).SarosiBest_FVU = cases(i).SarosiPrimary_FVU;
            cases(i).SarosiBest_MaxResidual = cases(i).SarosiPrimary_MaxResidual;
        else
            cases(i).SarosiBest = cases(i).SarosiAlt;
            cases(i).SarosiBest_kmax = cases(i).SarosiAlt_kmax;
            cases(i).SarosiBestLabel = 'Sarosi';
            cases(i).SarosiBest_RMSE = cases(i).SarosiAlt_RMSE;
            cases(i).SarosiBest_FVU = cases(i).SarosiAlt_FVU;
            cases(i).SarosiBest_MaxResidual = cases(i).SarosiAlt_MaxResidual;
        end
    else
        cases(i).SarosiBest = cases(i).SarosiPrimary;
        cases(i).SarosiBest_kmax = cases(i).SarosiPrimary_kmax;
        cases(i).SarosiBestLabel = cases(i).SarosiPrimaryLabel;
        cases(i).SarosiBest_RMSE = cases(i).SarosiPrimary_RMSE;
        cases(i).SarosiBest_FVU = cases(i).SarosiPrimary_FVU;
        cases(i).SarosiBest_MaxResidual = cases(i).SarosiPrimary_MaxResidual;
    end
end

%% Compare MB inputs for 20 mm, 300 mm and 450 mm
idx300 = find([cases.Diameter_mm] == 20 & abs([cases.RestLength_m] - 0.300) < 1e-12, 1, 'first');
idx450 = find([cases.Diameter_mm] == 20 & abs([cases.RestLength_m] - 0.450) < 1e-12, 1, 'first');

if ~isempty(idx300)
    x300 = cases(idx300).x(:);
    p300_kPa = cases(idx300).p_kPa(:);
    L300 = cases(idx300).RestLength_m - x300 .* (cases(idx300).RestLength_m - cases(idx300).Lmin_m);

    fprintf('\n--- MB input summary: 20 mm, 300 mm ---\n')
    fprintf('d used = %.3f mm\n', cases(idx300).MartensDiameter_mm)
    fprintf('Lrest = %.3f m, Lmin = %.3f m\n', cases(idx300).RestLength_m, cases(idx300).Lmin_m)
    fprintf('x range    = [%.4f, %.4f]\n', min(x300), max(x300))
    fprintf('p range    = [%.1f, %.1f] kPa\n', min(p300_kPa), max(p300_kPa))
    fprintf('Lm range   = [%.4f, %.4f] m\n', min(L300), max(L300))
end

if ~isempty(idx450)
    x450 = cases(idx450).x(:);
    p450_kPa = cases(idx450).p_kPa(:);
    L450 = cases(idx450).RestLength_m - x450 .* (cases(idx450).RestLength_m - cases(idx450).Lmin_m);

    fprintf('\n--- MB input summary: 20 mm, 450 mm ---\n')
    fprintf('d used = %.3f mm\n', cases(idx450).MartensDiameter_mm)
    fprintf('Lrest = %.3f m, Lmin = %.3f m\n', cases(idx450).RestLength_m, cases(idx450).Lmin_m)
    fprintf('x range    = [%.4f, %.4f]\n', min(x450), max(x450))
    fprintf('p range    = [%.1f, %.1f] kPa\n', min(p450_kPa), max(p450_kPa))
    fprintf('Lm range   = [%.4f, %.4f] m\n', min(L450), max(L450))
end

%% Debug MB on 20 mm, 450 mm case
idx450_dbg = find([cases.Diameter_mm] == 20 & abs([cases.RestLength_m] - 0.450) < 1e-12, 1, 'first');

if ~isempty(idx450_dbg)
    x_dbg = cases(idx450_dbg).x(:);
    p_dbg = cases(idx450_dbg).p_kPa(:) / 620;
    p_dbg_Pa = p_dbg * 6.2e5;
    L_dbg = cases(idx450_dbg).RestLength_m - x_dbg .* (cases(idx450_dbg).RestLength_m - cases(idx450_dbg).Lmin_m);

    [F_dbg, Fp_dbg, Fpe_dbg, Fl_dbg, D_dbg, rad_dbg] = Martens_force_from_PL( ...
        cases(idx450_dbg).RestLength_m, ...
        cases(idx450_dbg).MartensDiameter_mm, ...
        cases(idx450_dbg).MartensCoeff, ...
        p_dbg_Pa, ...
        L_dbg);

    disp('--- MB debug: 20 mm, 450 mm ---')
    disp(table(x_dbg, cases(idx450_dbg).p_kPa(:), L_dbg, F_dbg, Fp_dbg, Fpe_dbg, Fl_dbg, D_dbg, rad_dbg, ...
        'VariableNames', {'x','p_kPa','L_m','F_total','F_p','F_pe','F_l','D_m','rad'}))
end
%% Debug MB on 20 mm, 300 mm case
idx300_dbg = find([cases.Diameter_mm] == 20 & abs([cases.RestLength_m] - 0.300) < 1e-12, 1, 'first');

if ~isempty(idx300_dbg)
    x_dbg = cases(idx300_dbg).x(:);
    p_dbg = cases(idx300_dbg).p_kPa(:) / 620;
    p_dbg_Pa = p_dbg * 6.2e5;
    L_dbg = cases(idx300_dbg).RestLength_m - x_dbg .* (cases(idx300_dbg).RestLength_m - cases(idx300_dbg).Lmin_m);

    [F_dbg, Fp_dbg, Fpe_dbg, Fl_dbg, D_dbg, rad_dbg] = Martens_force_from_PL( ...
        cases(idx300_dbg).RestLength_m, ...
        cases(idx300_dbg).MartensDiameter_mm, ...
        cases(idx300_dbg).MartensCoeff, ...
        p_dbg_Pa, ...
        L_dbg);

    disp('--- MB debug: 20 mm, 300 mm ---')
    disp(table(x_dbg, cases(idx300_dbg).p_kPa(:), L_dbg, F_dbg, Fp_dbg, Fpe_dbg, Fl_dbg, D_dbg, rad_dbg, ...
        'VariableNames', {'x','p_kPa','L_m','F_total','F_p','F_pe','F_l','D_m','rad'}))
end
%% Main summary table
CaseLabel = strings(nCase,1);
Diameter_mm = zeros(nCase,1);
RestLength_mm = zeros(nCase,1);

OurModel_RMSE = zeros(nCase,1);
OurModel_FVU = zeros(nCase,1);
OurModel_MaxResidual = zeros(nCase,1);

Martens_RMSE = zeros(nCase,1);
Martens_FVU = zeros(nCase,1);
Martens_MaxResidual = zeros(nCase,1);

SarosiLabel = strings(nCase,1);
Sarosi_RMSE = zeros(nCase,1);
Sarosi_FVU = zeros(nCase,1);
Sarosi_MaxResidual = zeros(nCase,1);

for i = 1:nCase
    CaseLabel(i) = string(cases(i).CaseLabel);
    Diameter_mm(i) = cases(i).Diameter_mm;
    RestLength_mm(i) = cases(i).RestLength_mm;

    OurModel_RMSE(i) = cases(i).OurModel_RMSE;
    OurModel_FVU(i) = cases(i).OurModel_FVU;
    OurModel_MaxResidual(i) = cases(i).OurModel_MaxResidual;

    Martens_RMSE(i) = cases(i).Martens_RMSE;
    Martens_FVU(i) = cases(i).Martens_FVU;
    Martens_MaxResidual(i) = cases(i).Martens_MaxResidual;

    SarosiLabel(i) = string(cases(i).SarosiBestLabel);
    Sarosi_RMSE(i) = cases(i).SarosiBest_RMSE;
    Sarosi_FVU(i) = cases(i).SarosiBest_FVU;
    Sarosi_MaxResidual(i) = cases(i).SarosiBest_MaxResidual;
end

ResultsTable = table( ...
    CaseLabel, Diameter_mm, RestLength_mm, ...
    OurModel_RMSE, OurModel_FVU, OurModel_MaxResidual, ...
    SarosiLabel, Sarosi_RMSE, Sarosi_FVU, Sarosi_MaxResidual, ...
    Martens_RMSE, Martens_FVU, Martens_MaxResidual, ...
    'VariableNames', { ...
    'Case','Diameter_mm','RestLength_mm', ...
    'Bolen_RMSE','Bolen_FVU','Bolen_MaxResidual', ...
    'SarosiModel','Sarosi_RMSE','Sarosi_FVU','Sarosi_MaxResidual', ...
    'Martens_RMSE','Martens_FVU','Martens_MaxResidual'});

disp(' ')
disp('=== Main comparison table ===')
disp(ResultsTable)

%% Focused Sarosi repeatability comparison at 20 mm, 300 mm
idx300 = find([cases.Diameter_mm] == 20 & abs([cases.RestLength_m] - 0.300) < 1e-12, 1, 'first');

if ~isempty(idx300) && cases(idx300).HasSarosiComparison
    Sarosi300Table = table( ...
        ["DMSP-20-200N"; "DMSP-20-400N"], ...
        [cases(idx300).SarosiPrimary_kmax; cases(idx300).SarosiAlt_kmax], ...
        [cases(idx300).SarosiPrimary_RMSE; cases(idx300).SarosiAlt_RMSE], ...
        [cases(idx300).SarosiPrimary_FVU; cases(idx300).SarosiAlt_FVU], ...
        [cases(idx300).SarosiPrimary_MaxResidual; cases(idx300).SarosiAlt_MaxResidual], ...
        'VariableNames', {'Model','kmax','RMSE','FVU','MaxResidual'});

    disp(' ')
    disp('=== Sarosi comparison at 20 mm, 300 mm ===')
    disp(Sarosi300Table)
else
    Sarosi300Table = table;
end

%% Long-format results table
rows = {};
for i = 1:nCase
    rows(end+1,:) = {cases(i).CaseLabel, 'Bolen', cases(i).Diameter_mm, cases(i).RestLength_mm, ...
        cases(i).OurModel_RMSE, cases(i).OurModel_FVU, cases(i).OurModel_MaxResidual};

    rows(end+1,:) = {cases(i).CaseLabel, char(cases(i).SarosiPrimaryLabel), cases(i).Diameter_mm, cases(i).RestLength_mm, ...
        cases(i).SarosiPrimary_RMSE, cases(i).SarosiPrimary_FVU, cases(i).SarosiPrimary_MaxResidual};

    rows(end+1,:) = {cases(i).CaseLabel, 'Martens-Boblan', cases(i).Diameter_mm, cases(i).RestLength_mm, ...
        cases(i).Martens_RMSE, cases(i).Martens_FVU, cases(i).Martens_MaxResidual};

    if cases(i).HasSarosiComparison
        rows(end+1,:) = {cases(i).CaseLabel, char(cases(i).SarosiAltLabel), cases(i).Diameter_mm, cases(i).RestLength_mm, ...
            cases(i).SarosiAlt_RMSE, cases(i).SarosiAlt_FVU, cases(i).SarosiAlt_MaxResidual};
    end
end

ResultsLong = cell2table(rows, ...
    'VariableNames', {'Case','Model','Diameter_mm','RestLength_mm','RMSE','FVU','MaxResidual'});

disp(' ')
disp('=== Long comparison table ===')
disp(ResultsLong)

%% Sarosi kmax and Martens Lmin table
CaseLabel_k = strings(nCase,1);
Diameter_k = zeros(nCase,1);
RestLength_k = zeros(nCase,1);
SarosiPrimaryLabel_k = strings(nCase,1);
SarosiPrimary_kmax = NaN(nCase,1);
SarosiAltLabel_k = strings(nCase,1);
SarosiAlt_kmax = NaN(nCase,1);
Martens_Lmin = NaN(nCase,1);

for i = 1:nCase
    CaseLabel_k(i) = string(cases(i).CaseLabel);
    Diameter_k(i) = cases(i).Diameter_mm;
    RestLength_k(i) = cases(i).RestLength_m;
    SarosiPrimaryLabel_k(i) = string(cases(i).SarosiPrimaryLabel);
    SarosiPrimary_kmax(i) = cases(i).SarosiPrimary_kmax;
    Martens_Lmin(i) = cases(i).Lmin_used;

    if cases(i).HasSarosiComparison
        SarosiAltLabel_k(i) = string(cases(i).SarosiAltLabel);
        SarosiAlt_kmax(i) = cases(i).SarosiAlt_kmax;
    else
        SarosiAltLabel_k(i) = "";
    end
end

KmaxLminTable = table( ...
    CaseLabel_k, Diameter_k, RestLength_k, ...
    SarosiPrimaryLabel_k, SarosiPrimary_kmax, ...
    SarosiAltLabel_k, SarosiAlt_kmax, ...
    Martens_Lmin, ...
    'VariableNames', {'Case','Diameter_mm','RestLength_m', ...
    'SarosiPrimaryModel','SarosiPrimary_kmax', ...
    'SarosiAltModel','SarosiAlt_kmax', ...
    'Martens_Lmin_m'});

disp(' ')
disp('=== Sarosi kmax and Martens Lmin ===')
disp(KmaxLminTable)

%% Native Martens validation on the exact published grids
if MartensDOffset_mm == 0
    MBLegend = 'MB prediction';
    MBTitle = 'MB';
else
    MBLegend = sprintf('MB prediction (D_0 + %.1f mm)', MartensDOffset_mm);
    MBTitle = sprintf('MB (D_0 + %.1f mm)', MartensDOffset_mm);
end

d10_native = 10.0 + MartensDOffset_mm;
d20_native = 20.0 + MartensDOffset_mm;

F_M10_A2 = Martens_PL(0.250,d10_native,c10,P_A2,L_A2);
F_M20_A1 = Martens_PL(0.296,d20_native,c20,P_A1,L_A1);

F_M10_A2 = max(F_M10_A2,0);
F_M20_A1 = max(F_M20_A1,0);

errPct_A2 = 100 * max(abs(F_A2(:) - F_M10_A2(:))) / max(F_A2(:));
errPct_A1 = 100 * max(abs(F_A1(:) - F_M20_A1(:))) / max(F_A1(:));

alpha_A2 = (F_A2(:)' * F_M10_A2(:)) / (F_M10_A2(:)' * F_M10_A2(:));
alpha_A1 = (F_A1(:)' * F_M20_A1(:)) / (F_M20_A1(:)' * F_M20_A1(:));

F_M10_A2_scaled = alpha_A2 * F_M10_A2;
F_M20_A1_scaled = alpha_A1 * F_M20_A1;

errPct_A2_scaled = 100 * max(abs(F_A2(:) - F_M10_A2_scaled(:))) / max(F_A2(:));
errPct_A1_scaled = 100 * max(abs(F_A1(:) - F_M20_A1_scaled(:))) / max(F_A1(:));

disp('--- Martens scaling diagnostic ---')
disp([alpha_A2, errPct_A2, errPct_A2_scaled])
disp([alpha_A1, errPct_A1, errPct_A1_scaled])

[PP_A2,LL_A2g] = ndgrid(P_A2/1000, L_A2);
[PP_A1,LL_A1g] = ndgrid(P_A1/1000, L_A1);

figMartensNative = figure('Color','w','Name','Martens native validation');
tiledlayout(2,2)

nexttile
scatter3(LL_A2g(:), PP_A2(:), F_A2(:), ms, 'filled','MarkerFaceColor',colExp); hold on
scatter3(LL_A2g(:), PP_A2(:), F_M10_A2(:), ms, 'filled','MarkerFaceColor',colMartens);
xlabel('Length (m)')
ylabel('Pressure (kPa)')
zlabel('Force (N)')
title(sprintf('DMSP-10-250: MB Experiment vs MB prediction, d = %.1f mm', d10_native))
legend('MB Experiment','MB prediction','Location','best')

nexttile
resid = F_A2(:)-F_M10_A2(:);
scatter3(LL_A2g(:), PP_A2(:), resid, ms, 'filled')
xlabel('Length (m)')
ylabel('Pressure (kPa)')
zlabel('Residual (N)')
title(sprintf('DMSP-10-250 residual, err = %.2f%%', errPct_A2))

nexttile
scatter3(LL_A1g(:), PP_A1(:), F_A1(:), ms, 'filled','MarkerFaceColor',colExp); hold on
scatter3(LL_A1g(:), PP_A1(:), F_M20_A1(:), ms, 'filled','MarkerFaceColor',colMartens);
xlabel('Length (m)')
ylabel('Pressure (kPa)')
zlabel('Force (N)')
title(sprintf('DMSP-10-300: MB Experiment vs MB prediction, d = %.1f mm', d20_native))
legend('MB Experiment','MB prediction','Location','best')

nexttile
scatter3(LL_A1g(:), PP_A1(:), F_A1(:)-F_M20_A1(:), ms, 'filled')
xlabel('Length (m)')
ylabel('Pressure (kPa)')
zlabel('Residual (N)')
title(sprintf('DMSP-20-300 residual, err = %.2f%%', errPct_A1))

%% 20 mm debug sweep
dSweep_1 = 20.0 + MartensDOffset_mm;
dSweep_2 = 21.8 + MartensDOffset_mm;

F_300_20   = max(Martens_PL(0.300, dSweep_1, c20, P_A1, L_A1), 0);
F_300_21p8 = max(Martens_PL(0.300, dSweep_2, c20, P_A1, L_A1), 0);
F_296_20   = max(Martens_PL(0.296, dSweep_1, c20, P_A1, L_A1), 0);
F_296_21p8 = max(Martens_PL(0.296, dSweep_2, c20, P_A1, L_A1), 0);

figMartensSweep = figure('Color','w','Name','MB A1 debug sweep');
tiledlayout(2,2)

Fset = {F_300_20, F_300_21p8, F_296_20, F_296_21p8};
titles = { ...
    sprintf('L0=0.300, d=%.1f', dSweep_1), ...
    sprintf('L0=0.300, d=%.1f', dSweep_2), ...
    sprintf('L0=0.296, d=%.1f', dSweep_1), ...
    sprintf('L0=0.296, d=%.1f', dSweep_2)};

[PP_A1dbg,LL_A1dbg] = ndgrid(P_A1/1000, L_A1);

for k = 1:4
    nexttile
    scatter3(LL_A1dbg(:), PP_A1dbg(:), F_A1(:), ms, 'filled','MarkerFaceColor',colExp); hold on
    scatter3(LL_A1dbg(:), PP_A1dbg(:), Fset{k}(:), ms, 'filled','MarkerFaceColor',colMartens);
    xlabel('Length (m)')
    ylabel('Pressure (kPa)')
    zlabel('Force (N)')
    title(titles{k})
    legend('MB Experimental','MB prediction','Location','best')
end

%% RMSE figure
figRMSE = figure('Color','w','Name','RMSE summary by case');
tiledlayout(ceil(nCase/2),2)

for i = 1:nCase
    nexttile

    vals = [cases(i).OurModel_RMSE, cases(i).SarosiPrimary_RMSE, cases(i).Martens_RMSE];
    labels = {'Bolen', char(cases(i).SarosiPrimaryLabel), 'MB'};

    if cases(i).HasSarosiComparison
        vals = [cases(i).OurModel_RMSE, cases(i).SarosiPrimary_RMSE, cases(i).SarosiAlt_RMSE, cases(i).Martens_RMSE];
        labels = {'Bolen', char(cases(i).SarosiPrimaryLabel), char(cases(i).SarosiAltLabel), 'MB'};
    end

    bar(vals)
    set(gca,'XTick',1:numel(labels),'XTickLabel',labels)
    ylabel('RMSE')
    title(cases(i).CaseLabel)
end

%% Measured vs predicted 2D scatter
figScatter2D = figure('Color','w','Name','Measured vs predicted force by case','Renderer','painters');
set(figScatter2D,'Units','pixels','Position',[100 100 980 780],'PaperPositionMode','auto')
tiledlayout(2,2,'TileSpacing','loose','Padding','loose')

letters = {'(a)','(b)','(c)','(d)'};

for i = 1:nCase
    ax = nexttile;
    hold(ax,'on')

    scatter(cases(i).z, cases(i).OurModel, ms, 'filled','MarkerFaceColor',colBolen)

    if cases(i).HasSarosiComparison && ~useSarosiBest
        scatter(cases(i).z, cases(i).SarosiPrimary, ms, 'filled','MarkerFaceColor',colSarosiA)
        scatter(cases(i).z, cases(i).SarosiAlt, ms, 'filled','MarkerFaceColor',colSarosiB)
    else
        scatter(cases(i).z, cases(i).SarosiBest, ms, 'filled','MarkerFaceColor',colSarosiA)
    end

    scatter(cases(i).z, cases(i).Martens, ms, 'filled','MarkerFaceColor',colMartens)

    allPred = [cases(i).OurModel(:); cases(i).Martens(:)];
    if cases(i).HasSarosiComparison && ~useSarosiBest
        allPred = [allPred; cases(i).SarosiPrimary(:); cases(i).SarosiAlt(:)];
    else
        allPred = [allPred; cases(i).SarosiBest(:)];
    end

    limmax = max([cases(i).z(:); allPred(:)]);

    if i == 4
        xlim([0 limmax])
        ylim([-limmax limmax])
        plot([0 limmax],[0 limmax],'LineStyle','--', ...
            'Color',[128 128 128]/255,'LineWidth',2.5)
    else
        xlim([0 limmax])
        ylim([0 limmax])
        plot([0 limmax],[0 limmax],'LineStyle','--', ...
            'Color',[128 128 128]/255,'LineWidth',2.5)
    end

    axis square

    if i <= 2
        xticks([0 100 200 300 400 500 600])
        yticks([0 100 200 300 400 500 600])
        xlim([0 600])
        ylim([0 600])
    end

    hxl = xlabel('Experimental data force (N)');
    hxl.Units = 'normalized';
    hxl.Position(2) = hxl.Position(2) - 0.07;

    hyl = ylabel('Predicted force (N)');
    hyl.Units = 'normalized';
    hyl.Position(1) = hyl.Position(1) - 0.12;

    title(sprintf('\\phi%d mm, l_{rest} = %.3f m', ...
        cases(i).Diameter_mm, cases(i).RestLength_m), ...
        'Interpreter','tex')

    if cases(i).HasSarosiComparison && ~useSarosiBest
        legend('Bolen prediction', ...
            [char(cases(i).SarosiPrimaryLabel) ' prediction'], ...
            [char(cases(i).SarosiAltLabel) ' prediction'], ...
            'MB prediction', ...
            'Measured = predicted', ...
            'Location','southeast')
    else
        legend('Bolen prediction', ...
            'Sarosi prediction', ...
            'MB prediction', ...
            'Measured = predicted', ...
            'Location','southeast')
    end

    text(ax,0.5,-0.22,letters{i}, ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontWeight','bold', ...
        'Clipping','off')
end



%% 3D comparison plots by diameter
idx10 = find([cases.Diameter_mm] == 10);
idx20 = find([cases.Diameter_mm] == 20);

if ~isempty(idx10)
    fig10 = figure('Color','w','Name','10 mm BPA comparison');
    tiledlayout(2,3)

    for k = 1:min(2,numel(idx10))
        i = idx10(k);

        nexttile
        scatter3(cases(i).x, cases(i).p_kPa, cases(i).z, ms, 'filled','MarkerFaceColor',colExp); hold on
        scatter3(cases(i).x, cases(i).p_kPa, cases(i).OurModel, ms, 'filled','MarkerFaceColor',colBolen)
        scatter3(cases(i).x, cases(i).p_kPa, cases(i).SarosiPrimary, ms, 'filled','MarkerFaceColor',colSarosiA)
        scatter3(cases(i).x, cases(i).p_kPa, cases(i).Martens, ms, 'filled','MarkerFaceColor',colMartens)
        xlabel('Relative contraction')
        ylabel('Pressure (kPa)')
        zlabel('Force (N)')
        title(cases(i).CaseLabel)
        legend('Experimental data','Bolen prediction','Sarosi prediction','MB prediction','Location','bestoutside')
    end

    nexttile
    hold on
    for k = 1:numel(idx10)
        i = idx10(k);
        scatter3(cases(i).x, cases(i).p_kPa, cases(i).z - cases(i).OurModel, ms, 'filled')
    end
    view(0,0);
    xlabel('Relative contraction')
    ylabel('Pressure (kPa)')
    zlabel('Residual (N)')
    title('10 mm residuals: Bolen')
    legend(string({cases(idx10).CaseLabel}),'Location','best')

    nexttile
    hold on
    for k = 1:numel(idx10)
        i = idx10(k);
        scatter3(cases(i).x, cases(i).p_kPa, cases(i).z - cases(i).SarosiPrimary, ms, 'filled')
    end
    view(0,0)
    xlabel('Relative contraction')
    ylabel('Pressure (kPa)')
    zlabel('Residual (N)')
    title('10 mm residuals: Sarosi')
    legend(string({cases(idx10).CaseLabel}),'Location','best')

    nexttile
    hold on
    for k = 1:numel(idx10)
        i = idx10(k);
        scatter3(cases(i).x, cases(i).p_kPa, cases(i).z - cases(i).Martens, ms, 'filled')
    end
    view(0,0)
    xlabel('Relative contraction')
    ylabel('Pressure (kPa)')
    zlabel('Residual (N)')
    title('10 mm residuals: Martens')
    legend(string({cases(idx10).CaseLabel}),'Location','best')
else
    fig10 = [];
end

if ~isempty(idx20)
    fig20 = figure('Color','w','Name','20 mm BPA comparison');
    tiledlayout(2,3)

    shown = min(2,numel(idx20));

    for k = 1:shown
        i = idx20(k);

        nexttile
        scatter3(cases(i).x, cases(i).p_kPa, cases(i).z, ms, 'filled','MarkerFaceColor',colExp); hold on
        scatter3(cases(i).x, cases(i).p_kPa, cases(i).OurModel, ms, 'filled','MarkerFaceColor',colBolen)

        if cases(i).HasSarosiComparison
            scatter3(cases(i).x, cases(i).p_kPa, cases(i).SarosiPrimary, ms, 'filled','MarkerFaceColor',colSarosiA)
            scatter3(cases(i).x, cases(i).p_kPa, cases(i).SarosiAlt, ms, 'filled','MarkerFaceColor',colSarosiB)
            scatter3(cases(i).x, cases(i).p_kPa, cases(i).Martens, ms, 'filled','MarkerFaceColor',colMartens)
            legendEntries = {'Experimental data','Bolen prediction', ...
                [char(cases(i).SarosiPrimaryLabel) ' prediction'], ...
                [char(cases(i).SarosiAltLabel) ' prediction'], ...
                'MB prediction'};
        else
            scatter3(cases(i).x, cases(i).p_kPa, cases(i).SarosiPrimary, ms, 'filled','MarkerFaceColor',colSarosiA)
            scatter3(cases(i).x, cases(i).p_kPa, cases(i).Martens, ms, 'filled','MarkerFaceColor',colMartens)
            legendEntries = {'Experimental data','Bolen prediction', ...
                [char(cases(i).SarosiPrimaryLabel) ' prediction'], ...
                'MB prediction'};
        end
        xlabel('Relative contraction')
        ylabel('Pressure (kPa)')
        zlabel('Force (N)')
        title(cases(i).CaseLabel)
        legend(legendEntries,'Location','bestoutside')
    end

    nexttile
    hold on
    for k = 1:numel(idx20)
        i = idx20(k);
        scatter3(cases(i).x, cases(i).p_kPa, cases(i).z - cases(i).OurModel, ms, 'filled')
    end
    view(0,0)
    xlabel('Relative contraction')
    ylabel('Pressure (kPa)')
    zlabel('Residual (N)')
    title('20 mm residuals: Bolen')
    legend(string({cases(idx20).CaseLabel}),'Location','best')

    nexttile
    hold on
    sarLegendAdded = false;
    for k = 1:numel(idx20)
        i = idx20(k);
        if cases(i).HasSarosiComparison
            scatter3(cases(i).x, cases(i).p_kPa, cases(i).z - cases(i).SarosiPrimary, ms, 'filled');
            scatter3(cases(i).x, cases(i).p_kPa, cases(i).z - cases(i).SarosiAlt, ms, 'filled');
            sarLegendAdded = true;
        else
            scatter3(cases(i).x, cases(i).p_kPa, cases(i).z - cases(i).SarosiPrimary, ms, 'filled')
        end
    end
    view(0,0)
    xlabel('Relative contraction')
    ylabel('Pressure (kPa)')
    zlabel('Residual (N)')
    title('20 mm residuals: Sarosi')
    if sarLegendAdded
        legend('Sarosi_A residuals','Sarosi_B residuals','Location','best')
    end

    nexttile
    hold on
    for k = 1:numel(idx20)
        i = idx20(k);
        scatter3(cases(i).x, cases(i).p_kPa, cases(i).z - cases(i).Martens, ms, 'filled')
    end
    view(0,0)
    xlabel('Relative contraction')
    ylabel('Pressure (kPa)')
    zlabel('Residual (N)')
    title('20 mm residuals: Martens')
    legend(string({cases(idx20).CaseLabel}),'Location','best')
else
    fig20 = [];
end

%% Focused repeatability figure for 20 mm, 300 mm Sarosi coefficients
if ~isempty(idx300) && cases(idx300).HasSarosiComparison
    figSarosi300 = figure('Color','w','Name','Sarosi repeatability at 20 mm, 300 mm');
    tiledlayout(1,2)

    nexttile
    vals = [cases(idx300).SarosiPrimary_RMSE, cases(idx300).SarosiAlt_RMSE];
    bar(vals)
    set(gca,'XTick',1:2,'XTickLabel',{'a20a','a20b'})
    ylabel('RMSE')
    title('20 mm, 300 mm Sarosi comparison')

    nexttile
    scatter(cases(idx300).z, cases(idx300).SarosiPrimary, ms, 'filled','MarkerFaceColor',colSarosiA); hold on
    scatter(cases(idx300).z, cases(idx300).SarosiAlt, ms, 'filled','MarkerFaceColor',colSarosiB)
    xlabel('Experimental data force (N)')
    ylabel('Predicted force (N)')
    title('Sarosi at 300 mm: L0 = 200 vs. 400')
    legend('a20a prediction','a20b prediction','Location','best')
else
    figSarosi300 = [];
end

%% Output
Goof = struct;
Goof.cases = cases;
Goof.ResultsTable = ResultsTable;
Goof.ResultsLong = ResultsLong;
Goof.Sarosi300Table = Sarosi300Table;
Goof.KmaxLminTable = KmaxLminTable;

Goof.figMartensNative = figMartensNative;
Goof.figMartensSweep = figMartensSweep;
Goof.figRMSE = figRMSE;
Goof.figScatter2D = figScatter2D;
Goof.fig10mm3D = fig10;
Goof.fig20mm3D = fig20;
Goof.figSarosi300 = figSarosi300;

Goof.Published.Martens_A2 = F_M10_A2;
Goof.Published.Martens_A1 = F_M20_A1;
Goof.Published.errPct_A2 = errPct_A2;
Goof.Published.errPct_A1 = errPct_A1;
Goof.Published.alpha_A2 = alpha_A2;
Goof.Published.alpha_A1 = alpha_A1;
Goof.Published.errPct_A2_scaled = errPct_A2_scaled;
Goof.Published.errPct_A1_scaled = errPct_A1_scaled;

end

%% ---------------------- LOCAL FUNCTIONS ----------------------

function [FF, kmax] = sar(Lrest, Dia, u, x, pbar)
% Sarosi model
% x is your contraction coordinate, used as k = x*kmax
% pbar is already in bar

    Lrest = Lrest; %Resting length is not used
    Dia = Dia; %Diameter is not used

    a = u(1);
    b = u(2);
    c = u(3);
    d = u(4);
    e = u(5);
    f = u(6);

    myfun = @(K) (a*6.2 + b).*exp(c*K) + d*6.2.*K + e*6.2 + f;
    kmax = fzero(myfun,[14 35]);

    x = x(:);
    pbar = pbar(:);

    k = x * kmax;

    FF = (a*pbar + b).*exp(c*k) + d*pbar.*k + e*pbar + f;
end

function [FF, Lmin_used] = Martens(rest, Lmin, d, c, x, pPa)
% Martens force model using your x convention and pressure already in Pa
%
% Inputs
%   rest : resting length [m]
%   d    : initial inner diameter [mm]
%   c    : coefficient vector [c0 c1 c2 c3 d0]
%   x    : relative contraction coordinate
%   pPa  : pressure [Pa]

    x = x(:);
    pPa = pPa(:);

    if numel(x) ~= numel(pPa)
        error('Martens: x and pPa must have the same number of elements.');
    end

    Lm = rest - x .* (rest - Lmin);    % Current muscle length

    FF = Martens_force_from_PL(rest, d, c, pPa, Lm);
    Lmin_used = Lmin;
end

function [F, Fp, FpeTerm, FlTerm, D, rad] = Martens_force_from_PL(L0, d, c, P, L)

    c0 = c(1);
    c1 = c(2);
    c2 = c(3);
    c3 = c(4);
    d0 = c(5);

    t0    = 28.6;           %initial angle
    tcorr = t0 + d0;        %corrected angle
    H0    = 0.0018;         %material thickness

    D0 = d / 1000;

    P = P(:);
    L = L(:);

    if numel(P) ~= numel(L)
        error('P and L must have the same number of elements.');
    end

    Lfiber = L0 / cosd(tcorr);
    n      = L0 * tand(tcorr) / (pi * D0);

    rad = Lfiber.^2 - L.^2;
    D    = sqrt(rad) ./ (n*pi);
    dVdL = (Lfiber.^2 - 3*L.^2) ./ (4*pi*n.^2);
    dDdL = -(L ./ sqrt(rad)) ./ (n*pi);

    Eru = c3.*L.^3 + c2.*L.^2 + c1.*L + c0;

    sigL  = Eru .* (L - L0) ./ L0;
    sigPE = Eru .* (D - D0) ./ D0;

    FL  = sigL  .* H0 .* pi .* D;
    FPE = sigPE .* H0 .* L  .* pi;

    Fp      = -P .* dVdL;
    FpeTerm = FPE .* dDdL;
    FlTerm  = -FL;

    F = Fp + FpeTerm + FlTerm;
end

function F = Martens_PL(L0, d, c, P, L)

    c0 = c(1);
    c1 = c(2);
    c2 = c(3);
    c3 = c(4);
    d0 = c(5);

    t0    = 28.6;
    tcorr = t0 + d0;
    H0    = 0.0018;
    D0    = d / 1000;

    Lfiber = L0 / cosd(tcorr);
    n      = L0 * tand(tcorr) / (pi * D0);

    [PP, LL] = ndgrid(P, L);

    rad = Lfiber.^2 - LL.^2;
    D    = sqrt(rad) ./ (n*pi);
    dVdL = (Lfiber.^2 - 3*LL.^2) ./ (4*pi*n.^2);
    dDdL = -(LL ./ sqrt(rad)) ./ (n*pi);

    Eru = c3.*LL.^3 + c2.*LL.^2 + c1.*LL + c0;

    sigL  = Eru .* (LL - L0) ./ L0;
    sigPE = Eru .* (D - D0) ./ D0;

    FL  = sigL  .* H0 .* pi .* D;
    FPE = sigPE .* H0 .* LL .* pi;

    F = -PP .* dVdL + FPE .* dDdL - FL;
end

function applyAxisStyle(a, ft, fw)
    set(gca,'LineWidth',a,'FontSize',ft,'FontWeight',fw)
end