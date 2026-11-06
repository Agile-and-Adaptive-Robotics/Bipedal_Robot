function fr_audit_states_20261105()
C=fr_setup();kin=fr_kinematics(C.osimStock);outdir=fullfile(C.camp,'audit_20261105');if ~isfolder(outdir);mkdir(outdir);end
files=dir(fullfile(C.camp,'run_1','checkpoints','*_compliant.mat'));
ov=containers.Map('KeyType','char','ValueType','any');
for i=1:numel(files)
 S=load(fullfile(files(i).folder,files(i).name));ss=strsplit(S.cfgSig,'|');
 ov(S.csName)=struct('diaStr',ss{1},'nBPA',str2double(ss{2}),'family',ss{3},'gain',str2double(ss{4}));
end
cases=fr_cases27(C,kin,ov);rows={};
for i=1:numel(cases)
 cs=cases(i);
 for br={'rigid','compliant'}
  f=fullfile(C.camp,'run_1','checkpoints',[cs.name '_' br{1} '.mat']);if ~isfile(f);continue;end
  S=load(f);qg=struct('passIdx',1,'qList',{cs.passes{1}.qStructs});
  try
   [m,meta]=fr_mechBuild(C,cs,S.design,kin,br{1},qg);
   contraction=m.cell.sContraction(:,1)/S.design.L0;
   r=struct('group',cs.name,'branch',br{1},'nominal_diameter_mm',str2double(cs.diaStr), ...
    'nBPA',cs.seed.nBPA,'cfg',m.cfg,'L0_m',S.design.L0,'tendon_m',S.design.tendon, ...
    'physicalContraction',contraction,'legacyStrain',m.strain_p,'Fsingle_N',m.FsingleBPA,'Ftendon_N',m.Ftendon, ...
    'rowBodies',{meta.rowBodies},'loc',meta.loc,'poses',{qg.qList},'primaryDof',cs.passes{1}.sweptCoord, ...
    'qDefault',kin.qDefault,'status',S.status,'recordedFractionMet',S.rep.fractionPosesMet, ...
    'maxLocationDifference_m',max(abs(meta.loc-S.rep.loc),[],'all'), ...
    'limit',1-m.kmax/m.rest,'cross',cs.crossFor{1});
   rows{end+1}=r;
   fprintf('STATE %s %s physical %.4g..%.4g legacy %.4g..%.4g\n',cs.name,br{1},min(contraction),max(contraction),min(m.strain_p),max(m.strain_p));
  catch ME
   fprintf(2,'STATE ERROR %s %s: %s\n',cs.name,br{1},ME.message);
  end
 end
end
fid=fopen(fullfile(outdir,'mechanical_states.json'),'w');assert(fid>0);fwrite(fid,jsonencode(rows));fclose(fid);
fprintf('AUDIT COMPLETE %d states\n',numel(rows));
end
