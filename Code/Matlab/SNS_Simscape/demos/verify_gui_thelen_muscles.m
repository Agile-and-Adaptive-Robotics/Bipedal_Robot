function verify_gui_thelen_muscles(outputFolder)
% Compare the exact equations embedded in both Simulink muscle blocks with
% actual OpenSim calculations at 18 equilibrium and 36 moving-fiber states.
here=fileparts(mfilename('fullpath'));j=jsondecode(fileread(fullfile(here,'gait2392_reflex_geometry.json')));
rows=[j.opensim_equilibria(:);j.opensim_dynamic_states(:)];
df=zeros(numel(rows),1);dv=df;
for i=1:numel(rows)
    r=rows(i);x=j.muscles.(r.muscle);
    p=[x.Fmax x.l_opt x.l_slack deg2rad(x.alpha_deg) x.v_max x.FmaxTendonStrain ...
        x.FmaxMuscleStrain x.KshapePassive x.KshapeActive x.Af x.Flen x.fv_linear_extrap_threshold x.minimum_activation];
    [force,velocity]=thelen_reflex_fiber(r.activation,x.length_m(r.angle+1),r.fiber_length_m,p);
    df(i)=force-r.force_N;dv(i)=velocity-r.fiber_velocity_m_s;
end
result=struct('n_static',numel(j.opensim_equilibria),'n_dynamic',numel(j.opensim_dynamic_states), ...
    'max_force_difference_N',max(abs(df)),'max_velocity_difference_m_s',max(abs(dv)), ...
    'opensim_version',j.opensim,'source_sha256',j.source_sha256);
fid=fopen(fullfile(outputFolder,'thelen_reference_checks.json'),'w');
fwrite(fid,jsonencode(result,'PrettyPrint',true));fclose(fid);disp(result);
assert(max(abs(df))<1e-7 && max(abs(dv))<1e-8,'OpenSim/Simulink equations disagree');
end
