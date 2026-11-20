function check_extensor_replay
% Recalculate the row-32 results and match every stored value, including NaNs.
sourceDir = fileparts(mfilename('fullpath'));
stored = load(fullfile(sourceDir,'extensor_data.mat'),'tests','terms');
scratch = [tempname '.mat'];
export_extensor_data(scratch);
cleanup = onCleanup(@() delete(scratch)); %#ok<NASGU>
fresh = load(scratch,'tests','terms');
assert(isequaln(stored.tests,fresh.tests) && isequaln(stored.terms,fresh.terms), ...
    'The stored test arrays differ from the adopted MATLAB calculation.');
assert(isequal([fresh.tests.id],[5 6 7 1 2 8]));
trial = fresh.tests(3);
assert(isequal(trial.measured_angle(trial.measured_angle < -120),[-125;-123;-120.5]));
assert(isequal(trial.hybrid_angle(trial.hybrid_angle < -120),[-125;-123;-120.5]));
assert(isequal(trial.measured(trial.measured_angle == 2),-0.057265164));
count = arrayfun(@(s) numel(s.measured),fresh.tests);
hybridCount = arrayfun(@(s) numel(s.hybrid),fresh.tests);
assert(sum(count) == 113 && sum(hybridCount) == 113);
fprintf('Pinned extensor: %d measured, %d hybrid, and all six model pairs match.\n', ...
    sum(count),sum(hybridCount));
end
