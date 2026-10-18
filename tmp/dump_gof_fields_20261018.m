% dump_gof_fields_20261018.m - read-only provenance dump for torque tables 4.2/4.3.
% Prints stored per-member fields of the two adopted mats at the picks (77/107, 32)
% plus baseline (zero-correction) fields. No mats written.
root = 'D:\Github\Bipedal_Robot';
cd(fullfile(root,'Testing_Data','2022_02_Festo'));

S = load('minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat');
fprintf('\n=== FLEXOR mat: minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat ===\n');
disp('top-level fields:'); disp(fieldnames(S)');
fprintf('filtered_results size: %s\n', mat2str(size(S.filtered_results)));
if isfield(S,'xCols'), fprintf('xCols = %s\n', mat2str(S.xCols)); end
fn = fieldnames(S.filtered_results);
for f = 1:numel(fn)
  try
    v77 = S.filtered_results(77).(fn{f});
    v107 = S.filtered_results(107).(fn{f});
    if isnumeric(v77) && isnumeric(v107)
      fprintf('FLX %-24s row77=%-34s row107=%s\n', fn{f}, mat2str(v77,6), mat2str(v107,6));
    end
  catch
  end
end
if isfield(S,'baselineScores'), disp('FLEXOR baselineScores:'); disp(S.baselineScores); end

E = load('minimizeExt10mmX3_results_20260920_flxr77.mat');
fprintf('\n=== EXTENSOR mat: minimizeExt10mmX3_results_20260920_flxr77.mat ===\n');
disp('top-level fields:'); disp(fieldnames(E)');
fprintf('filtered_results size: %s\n', mat2str(size(E.filtered_results)));
if isfield(E,'xCols'), fprintf('xCols = %s\n', mat2str(E.xCols)); end
fe = fieldnames(E.filtered_results);
for f = 1:numel(fe)
  try
    v32 = E.filtered_results(32).(fe{f});
    if isnumeric(v32)
      fprintf('EXT %-24s row32=%s\n', fe{f}, mat2str(v32,6));
    end
  catch
  end
end
if isfield(E,'baselineScores'), disp('EXTENSOR baselineScores:'); disp(E.baselineScores); end
for lockf = {'LOCKSRC','xi1lock','xi2lock','flxRow'}
  if isfield(E,lockf{1}), fprintf('%s = %s\n', lockf{1}, mat2str(E.(lockf{1}),10)); end
end
fprintf('\nDUMP DONE\n');
