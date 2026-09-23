%% sns_build_circuit_view.m — KneeReflexCircuit.slx (2026-09-22)
%
% The 2026-09-22 library redesign made the KneeReflexDemo top level itself
% print-clean: compact one-input synapses sit against the neurons they
% synapse onto, summation lives inside the neurons, and the mechanics are one
% masked KneeModel block. So the circuit VIEW is now a runnable 1:1 TWIN of
% the demo: same blocks, same mask values, same wiring — it simulates
% standalone AND exports as the dissertation circuit figure
% (sns_export_diagram prints either model).
%
% (The 2026-09-20 version was a print-only copy with stubbed In/Out ports —
% it could not run, which Ben flagged: "neither work".)

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                       % SNS_Simscape (SNS_Library)
cd(here);

src = 'KneeReflexDemo';
mdl = 'KneeReflexCircuit';
load_system('SNS_Library');
load_system(fullfile(here, [src '.slx']));
if bdIsLoaded(mdl), close_system(mdl, 0); end
f = fullfile(here, [mdl '.slx']);
if exist(f, 'file'), delete(f); end

% 1:1 twin: save the loaded demo under the circuit-view name (the demo file
% on disk is untouched; the in-memory model is renamed by save_system).
save_system(src, f);
try
    anno = Simulink.Annotation(mdl, ...
        'KneeReflexCircuit — runnable 1:1 circuit view of KneeReflexDemo (identical blocks, values, wiring)');
    anno.Position = [40 -115 900 -75];
catch
end
save_system(mdl);
close_system(mdl, 0);
close_system(src, 0);
fprintf('%s.slx built (runnable twin of %s)\n', mdl, src);
