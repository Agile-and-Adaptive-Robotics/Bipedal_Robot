%% sns_build_simscape_lib.m — compile the +SNS Simscape package into SNS_lib.slx
% Blocks become available as SNS_lib/NonSpikingNeuron and SNS_lib/NonSpikingSynapse.
% Requires Simscape (licensed on this machine; Simscape Multibody NOT needed).

thisDir = fileparts(mfilename('fullpath'));
cd(fullfile(thisDir, 'simscape_sources'));

if exist('SNS_lib', 'file') || exist([thisDir filesep 'simscape_sources' filesep 'SNS_lib.slx'], 'file')
    % ssc_build overwrites anyway
end
ssc_build('SNS');
fprintf('SNS_lib built at %s\n', fullfile(pwd, 'SNS_lib.slx'));

%% post-process: give blocks clean single-line names (ssc_build wraps long names)
libFile = fullfile(thisDir, 'simscape_sources', 'SNS_lib.slx');
fileattrib(libFile, '+w');           % ssc_build writes it read-only
load_system('SNS_lib');
set_param('SNS_lib', 'LibraryLock', 'off');   % ssc_build locks the library
rename = {'NonSpikingNeuron', 'NonSpikingSynapse'};
hits = find_system('SNS_lib', 'SearchDepth', 1, 'Type', 'Block');
for k = 1:numel(rename)
    for h = 1:numel(hits)
        nm = strrep(strrep(get_param(hits{h}, 'Name'), newline, ' '), ' ', '');
        if contains(nm, rename{k})
            set_param(hits{h}, 'Name', rename{k});
            fprintf('renamed block -> %s\n', rename{k});
        end
    end
end
save_system('SNS_lib');
close_system('SNS_lib', 0);
fprintf('SNS_lib post-processed.\n');
