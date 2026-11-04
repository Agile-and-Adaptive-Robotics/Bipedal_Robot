function resArr = fr_scanCheckpoints(ckDir, CASES)
%FR_SCANCHECKPOINTS Resume scan: one row per case over the checkpoint dir.
resArr = repmat(fr_resRow(), 1, numel(CASES));
for i = 1:numel(CASES)
    cs = CASES(i);
    resArr(i).name = cs.name;
    resArr(i).mapId = cs.groupMapId;
    ckR = fr_ckPath(cs.name, 'rigid', ckDir);
    ckC = fr_ckPath(cs.name, 'compliant', ckDir);
    ckF = fr_ckPath(cs.name, 'refine1', ckDir);
    resArr(i).ckRigid = ckR;
    resArr(i).ckCompliant = ckC;
    resArr(i).ckRefine = ckF;
    resArr(i).hasRigid = exist(ckR, 'file') == 2;
    resArr(i).hasCompliant = exist(ckC, 'file') == 2;
    resArr(i).hasRefine = exist(ckF, 'file') == 2;
end
end
