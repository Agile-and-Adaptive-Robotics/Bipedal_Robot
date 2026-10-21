% Write companion xlsx for the supplement BPA routing tables (S3-S5).
% Values identical to frontiers_SupplementaryMaterial.tex; sources cited in
% the supplement-builder report. SCRATCH (tmp/paper_bridge), 2026-10-20.
outDir = 'D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/Tables';

%% S3: pinned-knee bracket reference points
Configuration = {'Pinned flexor';'Pinned flexor';'Pinned extensor'};
Point = {'Pbri (insertion bracket)';'Pbr2 (origin-side second bracket)';'Pbr (rib midpoint, medial side)'};
Frame = {'knee ICR';'hip';'hip'};
x_mm = [-48.11; -52.61; -3.84];
y_mm = [-107.81; 0; -46.44];
z_mm = [13.80; 75.06; 62.50];
Notes = {'flexion-negative arc R1=12 mm active -23<=thk<=27 deg; R2=40 mm active thk<=-23 deg (extensor)';
         ''; 'extensor wrap arcs R1=12 mm (-23..27 deg), R2=40 mm (thk<=-23 deg)'};
T3 = table(Configuration, Point, Frame, x_mm, y_mm, z_mm, Notes);
writetable(T3, fullfile(outDir,'Route_PinnedBrackets.xlsx'));

%% S4: biomimetic 10 mm routes
cfg    = {'Flexor 10mm';'Flexor 10mm';'Extensor 10mm';'Extensor 10mm';'Extensor 10mm';'Extensor 10mm';'Extensor 10mm';'Extensor 10mm';'Extensor 10mm';'Extensor 10mm'};
pt     = {'p1';'p2';'p1';'p2';'p3';'p4';'p5';'p6';'p7';'p8'};
descr  = {'origin';'insertion';'origin';'mounting-base contact';'femur channel contact';'femoral condyle contact';'tibia contact';'tibia contact';'tibia tendon contact';'patellar ligament ring'};
frame  = {'femur';'tibia theta1';'femur';'femur';'femur';'femur';'tibia theta1';'tibia theta1';'tibia theta1';'tibia theta1'};
x4 = [-50.00; -12.24; 40.00; 75.90; 59.82; 39.55; 58.71; 58.71; 54.00; 36.04];
y4 = [ 35.00;  -8.87; 35.00;-274.76;-374.27;-421.83; 25.00; 12.28;  -6.12;-28.44];
z4 = [ 50.00;  27.87;  0.00;   0.00;   0.00;   0.00;  0.00;  0.00;   0.00;  0.00];
act = {'all';'all';'all';'all';'thk<=-9.19 deg';'thk<=-69.8 deg';'thk<=-9.19 deg';'all';'all';'all'};
restmm = [415;415;520;520;520;520;520;520;520;520];
T4 = table(cfg, pt, descr, frame, x4, y4, z4, act, restmm, ...
    'VariableNames', {'Configuration','Point','Description','Frame','x_mm','y_mm','z_mm','Active_thk','RestingLength_mm'});
writetable(T4, fullfile(outDir,'Route_Biomimetic10mm.xlsx'));

%% S5: 20 mm redesign routes
cfg5 = [repmat({'Flexor 20mm'},4,1); repmat({'Extensor 20mm'},9,1)];
pt5  = {'p1';'p_end';'p1_BPA2';'p_end_BPA2';'p1';'p2';'p3';'p4';'p5';'p6';'p7';'p8';'p9'};
descr5 = {'origin, BPA 1';'insertion, BPA 1';'origin, BPA 2';'insertion, BPA 2'; ...
          'origin (fixed)';'mounting-base contact';'femur contact';'femur contact'; 'femoral condyle contact'; ...
          'cam contact (optimized)';'cam contact (optimized)';'cam contact (optimized)';'distal ring (insertion)'};
frame5 = [repmat({'femur'},3,1); {'tibia theta1'}; repmat({'femur'},5,1); repmat({'tibia theta1'},4,1)];
frame5(3) = {'femur'};  % p1_BPA2 femur frame
x5 = [-125.00; -18.72; -125.00; -18.72; 40.00; 83.00; 66.60; 54.70; 38.80; 116.32; 131.32; 127.04; 34.42];
y5 = [ 135.00; -32.69;  135.00; -32.69; 137.00;-271.70;-427.80;-443.60;-454.40; 15.04; -18.97; -47.87;-197.10];
z5 = [ 72.08; 38.96;  -5.85; -38.96; 0; 0; 0; 0; 0; 0; 0; 0; 0];
lift = {'end point';'end point';'end point';'end point';'end point';'-21.52 deg';'-54.34 deg';'-89.80 deg';'-110.81 deg';'-18.89 deg';'none';'none';'end point'};
consts = [repmat({'2 BPAs, rest 561.49 mm, tendon 25.00 mm'},4,1); repmat({'2 BPAs, rest 827.74 mm, tendon 92.49 mm, CrossPoint 6'},9,1)];
T5 = table(cfg5, pt5, descr5, frame5, x5, y5, z5, lift, consts, ...
    'VariableNames', {'Configuration','Point','Description','Frame','x_mm','y_mm','z_mm','Liftoff_thk','RouteConstants'});
writetable(T5, fullfile(outDir,'Route_Redesign20mm.xlsx'));

fprintf('Wrote 3 xlsx files to %s\n', outDir);
