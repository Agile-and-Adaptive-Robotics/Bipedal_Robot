function fr_writeStatus(runDir, startedAt, stopped, layerResults, ...
    resArr, envRecord, CASES, refineOut, allIds)
%FR_WRITESTATUS run_0/STATUS.json - progress snapshot / final record.
% Called with 6 args while stopped early, 9 args after refinement.

st = struct();
st.run = 'run_0_weights_placeholder';
st.startedAt = startedAt;
st.updatedAt = datetime('now');
st.stoppedCleanly = stopped;
st.machine = getenv('COMPUTERNAME');
if nargin >= 6 && ~isempty(envRecord)
    st.placedEnvelopes = envRecord.Count;
end
st.groups = struct('mapId', {}, 'group', {}, 'rigidStatus', {}, ...
    'rigidF', {}, 'rigidViol', {}, 'compliantStatus', {}, ...
    'compliantF', {}, 'compliantViol', {}, 'refineVerdict', {}, ...
    'fPrev', {}, 'fNew', {});
if nargin >= 4 && ~isempty(layerResults)
    for L = 1:numel(layerResults)
        if isempty(layerResults{L}); continue; end
        po1 = layerResults{L};
        for i = 1:numel(po1)
            st.groups(end + 1) = struct('mapId', po1(i).mapId, ...
                'group', po1(i).name, ...
                'rigidStatus', po1(i).rigid.status, ...
                'rigidF', po1(i).rigid.fBest, ...
                'rigidViol', po1(i).rigid.violation, ...
                'compliantStatus', po1(i).compliant.status, ...
                'compliantF', po1(i).compliant.fBest, ...
                'compliantViol', po1(i).compliant.violation, ...
                'refineVerdict', 'pending', 'fPrev', NaN, 'fNew', NaN);
        end
    end
end
if nargin >= 9 && ~isempty(refineOut)
    gmap = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for i = 1:numel(st.groups)
        gmap(st.groups(i).mapId) = i;
    end
    for i = 1:numel(refineOut)
        ro = refineOut(i);
        if ~isKey(gmap, ro.mapId); continue; end
        gi = gmap(ro.mapId);
        st.groups(gi).refineVerdict = ro.verdict;
        st.groups(gi).fPrev = ro.fPrev;
        st.groups(gi).fNew = ro.fNew;
    end
end
% checkpoint presence as the failure-independent truth
if nargin >= 5 && ~isempty(resArr)
    for i = 1:numel(resArr)
        ra = resArr(i);
        idx = find([st.groups.mapId] == ra.mapId, 1);
        if isempty(idx); continue; end
        if ra.hasRigid && strcmp(st.groups(idx).rigidStatus, '')
            st.groups(idx).rigidStatus = 'checkpointed';
        end
        if ra.hasCompliant && strcmp(st.groups(idx).compliantStatus, '')
            st.groups(idx).compliantStatus = 'checkpointed';
        end
    end
end
f = fullfile(runDir, 'STATUS.json');
txt = jsonencode(st, 'PrettyPrint', true);
fid = fopen(f, 'w');
fwrite(fid, txt, 'char');
fclose(fid);
fprintf('status written: %s\n', f);
end
