function fr_screenAnkleExcursion(outputDir,stockDir)
% Apply the rectus excursion and force screen to full coupled ankle routes.
C=fr_setup();
kin=fr_kinematics(C.osimStock);
stage=fullfile(C.camp,'run_20261105_distal_stages');
groups={'ext_dig','ext_hal'}; folders={'ankle_guide_arrow','ankle_mounts_600'};
if ~isfolder(outputDir); mkdir(outputDir); end
for g=1:numel(groups)
    name=[groups{g} '_compliant'];
    S=load(fullfile(stage,folders{g},[name '.mat']),'result');
    c=S.result.chosen;
    H=readtable(fullfile(stockDir,[name '_coupled_torque.csv']));
    n=numel(unique(H.ankle_deg));
    G=fr_ankleCoupledSetup(C,kin,c.caseDef,n);
    assert(max(abs(H.ankle_deg-G.ankle_deg))<1e-9 && ...
        max(abs(H.mtp_angle_r_deg-G.second_deg))<1e-9,'Stock and robot poses differ');
    [m,meta]=fr_mechBuild(C,G.cs,c.design,kin,'compliant',G.qg);
    bay=sum(vecnorm(diff(meta.loc(1:c.caseDef.distalSleeveGuideRow,:)),2,2))-2*c.caseDef.fittingLength;
    stock=fr_excursionSizing(C,H.human_path_m,bay,H.human_Fmax_N(1));
    robot=fr_excursionSizing(C,m.cell.MuscleLength(:),bay,H.human_Fmax_N(1));
    stock.route=repmat("stock OpenSim",height(stock),1);
    robot.route=repmat("saved robot",height(robot),1);
    T=[stock;robot];
    T.half_Fmax10_screen=repmat(H.human_Fmax_N(1)<=stock.initial_10mm_half_Fmax_N(1),height(T),1);
    writetable(T,fullfile(outputDir,[groups{g} '_excursion_sizing.csv']));
    fprintf('%s longest stock path %.1f mm, stock excursion %.1f mm, robot excursion %.1f mm, bay %.1f mm\n', ...
        groups{g},1000*max(H.human_path_m),1000*range(H.human_path_m), ...
        1000*range(m.cell.MuscleLength(:)),1000*bay);
    disp(T(T.gain==1,{'route','diameter_mm','required_stroke_mm','minimum_rest_mm', ...
        'maximum_proxy_rest_mm','axial_screen_pass','parallel_count_screen'}));
end
end
