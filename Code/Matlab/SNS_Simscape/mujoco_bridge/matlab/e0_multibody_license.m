function e0_multibody_license()
% E0: can THIS machine build + RUN hand-built Simscape Multibody models
% (primitive blocks, no CAD import)? Builds World Frame -> Revolute Joint ->
% Cylindrical Solid, sims 0.1 s. Verdict printed + any license error shown.
% Also reports whether a generic Force Element exists (needed to later
% implement a muscle pulling between two frames - Ben's elbow idea).

try
    load_system('sm_lib');
catch e
    fprintf(['E0 VERDICT: FAIL - sm_lib will not even load: %s\n'], e.message);
    return
end
blks = find_system('sm_lib', 'SearchDepth', 3, 'Type', 'Block');

hits = @(pat) blks(~cellfun(@isempty, regexp(blks, pat, 'once')));
wf  = hits('World Frame');
rj  = hits('Joints[/\\]Revolute Joint');
cyl = hits('Cylindrical Solid');
fe  = hits('Force Element');
if ~isempty(fe)
    fprintf('library has Force Element: %s\n', fe{1});
else
    fprintf('no generic Force Element found in top-3 library levels\n');
end
if isempty(wf) || isempty(rj) || isempty(cyl)
    fprintf('E0 VERDICT: INCONCLUSIVE - block paths not found:\n');
    fprintf('  world=%d revolute=%d cylinder=%d\n', ...
            numel(wf), numel(rj), numel(cyl));
    return
end

mdl = 'sm_license_probe';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
add_block(wf{1}, [mdl '/WF'], 'Position', [100 100 130 130]);
add_block(rj{1}, [mdl '/RJ'], 'Position', [200 100 250 150]);
add_block(cyl{1}, [mdl '/Body'], 'Position', [320 100 390 160]);
% conserving connections by port index (B/F ordering is library-dependent,
% so try both orders)
ok = false; err = '';
    for pr = {{'WF/1','RJ/1','RJ/2','Body/1'}, ...
              {'WF/1','RJ/2','RJ/1','Body/1'}}
        for prl = {{'WF/1','RJ/1'}, {'WF/1','RJ/2'}, ...
                   {'RJ/1','Body/1'}, {'RJ/2','Body/1'}}
            try
                delete_line(mdl, prl{1}, prl{2});
            catch
            end
        end
        add_line(mdl, pr{1}, pr{2}, 'autorouting', 'on');
        add_line(mdl, pr{3}, pr{4}, 'autorouting', 'on');
        set_param(mdl, 'StopTime', '0.1');
        try
            sim(mdl);
            ok = true; break
        catch e
            err = e.message;
        end
    end
if ok
    fprintf(['E0 VERDICT: PASS - hand-built Simscape Multibody RUNS on this ' ...
             'machine (World+Revolute+Cylinder, 0.1 s clean)\n']);
else
    fprintf('E0 VERDICT: FAIL - sim error:\n%s\n', err);
end
try
    close_system(mdl, 0);
catch
end
bdclose('sm_lib');
end
