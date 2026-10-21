% Rebuild the three GoF xlsx companions in the l0/l_ten column format
tblDir = 'D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/Tables';

% PinnedFlexor_Preopt_GoF.xlsx
H = {'l0_cm', 'lten_mm', 'Hunt_RMSE', 'Hunt_FVU', 'Hunt_Max', ...
     'CharAlone_RMSE', 'CharAlone_FVU', 'CharAlone_Max', ...
     'Hybrid_RMSE', 'Hybrid_FVU', 'Hybrid_Max'};
D = [ ...
 48.5 NaN 5.779 1.1761 9.798 4.871 0.8354 9.865 1.635 0.0991 3.842; ...
 41.0 NaN 6.633 3.2313 11.409 6.074 2.7090 9.262 2.187 0.1184 4.243; ...
 45.7 NaN 5.315 0.9408 9.991 4.465 0.6640 9.536 1.246 0.0409 2.275; ...
 40.6 35 3.945 0.6750 7.527 3.226 0.4512 7.376 2.085 0.0755 3.073];
writetable(array2table(D, 'VariableNames', H), ...
    fullfile(tblDir, 'PinnedFlexor_Preopt_GoF.xlsx'));

% PinnedFlexor_GoF.xlsx (Training/Validation blocks stacked with a Group col)
G = [repmat("Training", 3, 1); repmat("Validation", 3, 1)];
l0 = [48.5 41.0 NaN 45.7 40.6 NaN]';
lt = [NaN NaN NaN NaN 35 NaN]';
lbl = [""; ""; "Mean"; ""; ""; "Mean"];
c1 = [7.617 9.694 8.656 8.803 6.373 7.588]';
c2 = [2.1278 2.1149 2.1213 1.9380 0.5935 1.2657]';
c3 = [25.531 21.076 23.303 26.264 16.859 21.562]';
c4 = [2.165 1.084 1.624 1.368 1.611 1.490]';
c5 = [0.1718 0.0277 0.0998 0.0468 0.0427 0.0447]';
c6 = [5.467 1.766 3.617 3.614 5.891 4.752]';
writetable(table(G, l0, lt, lbl, c1, c2, c3, c4, c5, c6, ...
    'VariableNames', {'Group', 'l0_cm', 'lten_mm', 'RowLabel', ...
    'CharAlone_RMSE', 'CharAlone_FVU', 'CharAlone_Max', ...
    'Corrected_RMSE', 'Corrected_FVU', 'Corrected_Max'}), ...
    fullfile(tblDir, 'PinnedFlexor_GoF.xlsx'));

% PinnedExtensor_GoF.xlsx
G = [repmat("Training", 4, 1); repmat("Validation", 4, 1)];
l0 = [43.6 43.2 45.7 NaN 40.0 40.0 46.5 NaN]';
lt = [NaN 22 NaN NaN NaN 40 NaN NaN]';
lbl = [""; ""; ""; "Mean"; ""; ""; ""; "Mean"];
c1 = [2.437 2.878 2.178 2.498 3.092 3.622 0.949 2.554]';
c2 = [0.9618 4.7304 1.0586 2.2503 0.5638 6.3919 0.3152 2.4236]';
c3 = [4.284 4.817 4.974 4.692 4.842 7.845 2.272 4.986]';
c4 = [0.978 1.525 0.490 0.997 1.244 1.350 0.699 1.098]';
c5 = [0.1634 1.3270 0.0535 0.5146 0.1303 0.8878 0.1714 0.3965]';
c6 = [1.818 2.481 1.206 1.835 2.196 1.890 1.559 1.882]';
writetable(table(G, l0, lt, lbl, c1, c2, c3, c4, c5, c6, ...
    'VariableNames', {'Group', 'l0_cm', 'lten_mm', 'RowLabel', ...
    'CharAlone_RMSE', 'CharAlone_FVU', 'CharAlone_Max', ...
    'Corrected_RMSE', 'Corrected_FVU', 'Corrected_Max'}), ...
    fullfile(tblDir, 'PinnedExtensor_GoF.xlsx'));

disp('ALL THREE XLSX REBUILT IN l0/lten FORMAT');
