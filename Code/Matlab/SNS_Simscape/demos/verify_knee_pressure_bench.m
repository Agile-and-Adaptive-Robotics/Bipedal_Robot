function summary = verify_knee_pressure_bench(repoRoot, outputFolder)
% Drive the exact Simulink BPA block with the Chapter 3 measured lengths
% and pressures, then compare independently with the published equations.
td=fullfile(repoRoot,'Testing_Data'); lib=fullfile(repoRoot,'Code','Matlab','SNS_Simscape');
addpath(lib,fullfile(repoRoot,'Code','Matlab','Functions'));
data=load(fullfile(td,'datacompare20.mat'));
if ~isfolder(outputFolder),mkdir(outputFolder);end
load_system('SNS_Library'); mdl='KneePressureBench';
if bdIsLoaded(mdl),close_system(mdl,0);end
new_system(mdl); cleanup=onCleanup(@()close_system(mdl,0));
add_block('simulink/Sources/From Workspace',[mdl '/Pressure'],'VariableName','PressureBenchInput','Interpolate','off','OutputAfterFinalValue','Holding final value');
add_block('simulink/Sources/From Workspace',[mdl '/Length'],'VariableName','LengthBenchInput','Interpolate','off','OutputAfterFinalValue','Holding final value');
add_block('SNS_Library/BPA_20mm',[mdl '/BPA'],'Rest','rest','Kmax','free620','TendonL','0','FittingL','0');
add_block('simulink/Sinks/To Workspace',[mdl '/Force'],'VariableName','benchForce','SaveFormat','Timeseries');
add_line(mdl,'Pressure/1','BPA/1');add_line(mdl,'Length/1','BPA/2');add_line(mdl,'BPA/1','Force/1');
set_param(mdl,'Solver','FixedStepDiscrete','FixedStep','0.01','ReturnWorkspaceOutputs','on');
for k=1:2
    if k==1, rest=0.300;free620=0.225;id='30';else,rest=0.450;free620=0.335;id='45';end
    rhoMeasured=data.(['E' id]); pressure=data.(['P' id]); measured=data.(['F' id]);
    length=rest-rhoMeasured*(rest-free620);
    % Read each held sample after the discontinuity.
    tick=(0:numel(pressure)-1)'*.02;
    in=Simulink.SimulationInput(mdl);
    in=in.setVariable('rest',rest).setVariable('free620',free620);
    in=in.setVariable('PressureBenchInput',[tick,pressure]).setVariable('LengthBenchInput',[tick,length]);
    in=in.setModelParameter('StopTime',num2str(tick(end)+.01,'%.8f'));
    out=sim(in); native=interp1(out.benchForce.Time,out.benchForce.Data,tick+.01,'nearest');
    computedRho=(rest-length)/(rest-free620);
    f620=maxBPAforce(rest,'20',620);
    raw=f620*(.257852586017*(exp(-6.4766142989*computedRho)-1)+(pressure/620).*exp(-1.32087718059*computedRho.^2));
    expected=min(f620,max(0,raw));expected(computedRho>=1)=0;
    delta=abs(native-expected);err=native-measured;
    summary(k)=struct('rest_m',rest,'free620_m',free620,'measured_pairs',numel(measured),...
        'F620_N',f620,'max_equation_difference_N',max(delta),'empirical_rmse_N',sqrt(mean(err.^2)),...
        'max_absolute_empirical_difference_N',max(abs(err)),'unclipped_fit_rmse_N',sqrt(mean((raw-measured).^2)));
    assert(summary(k).max_equation_difference_N<1e-7,'BPA block failed analytical dimensional check');
    bench=table(rhoMeasured,length,pressure,measured,native,expected,raw,err);
    writetable(bench,fullfile(outputFolder,['force_bench_' id '.csv']));
end
fid=fopen(fullfile(outputFolder,'force_bench.json'),'w');
fprintf(fid,'%s\n',jsonencode(summary,PrettyPrint=true));fclose(fid);disp(struct2table(summary));
end
