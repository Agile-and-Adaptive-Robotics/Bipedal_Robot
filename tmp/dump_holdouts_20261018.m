% dump_holdouts_20261018.m - read-only: which tests each adopted mat held out,
% and the drivers' own recorded picks. No evaluator calls, no parpool.
root = 'D:\Github\Bipedal_Robot';
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

S = load('minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat', ...
    'PICK','NUMHOLD','TAG','TRANSMODE','labels','allBPA','numBPA','results_cv','xCols');
fprintf('FLEXOR: PICK=%s NUMHOLD=%s TAG=%s TRANSMODE=%s numBPA=%s\n', ...
    mat2str(S.PICK), mat2str(S.NUMHOLD), mat2str(S.TAG), mat2str(S.TRANSMODE), mat2str(S.numBPA));
fprintf('FLEXOR labels: %s\n', strjoin(string(S.labels), ' | '));
rcv = S.results_cv;
fprintf('FLEXOR results_cv fields: %s\n', strjoin(fieldnames(rcv)', ', '));
for f = 1:numel(rcv)
    try
        holdidx = rcv(f).holdTest; train = rcv(f).trainTests;
        fprintf('fold %d: holdTest=%s trainTests=%s\n', f, mat2str(holdidx), mat2str(train));
    catch
        try
            fprintf('fold %d fields: %s\n', f, strjoin(fieldnames(rcv)', ', '));
        catch, end
    end
end

E = load('minimizeExt10mmX3_results_20260920_flxr77.mat', ...
    'exthold','holdTests','trainTests','holdoutIdx','validLabels','fold', ...
    'train','val','pick','labels','TAG','numBPA','NUMHOLD');
fprintf('\nEXTENSOR flxr77: pick=%s TAG=%s exthold=%s\n', ...
    mat2str(E.pick), mat2str(E.TAG), mat2str(E.exthold));
fprintf('EXT labels: %s\n', strjoin(string(E.labels), ' | '));
fprintf('EXT holdTests=%s trainTests=%s holdoutIdx=%s fold=%s\n', ...
    mat2str(E.holdTests), mat2str(E.trainTests), mat2str(E.holdoutIdx), mat2str(E.fold));
fprintf('EXT train=%s val=%s validLabels=%s\n', mat2str(E.train), mat2str(E.val), strjoin(string(E.validLabels), ' | '));
fprintf('\nHOLDOUT DUMP DONE\n');
