function calibrate_muscles_20260926()
% measure actual origin->insertion vector at t=0, correct origT/insT affinely,
% re-verify; target relative vectors (world axes, rig CAD pose):
%   EXT: (0.005, -0.130, 0) | FLX: (-0.005, -0.130, 0)
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

target = struct('EXT', {[0.005 -0.130 0]}, 'FLX', {[-0.005 -0.130 0]});

for pass = 1:3
    % enable xyz sensing
    for pfx = {'EXT', 'FLX'}
        set_param([mdl '/' pfx{1} '_TS'], 'SenseX', 'on', 'SenseY', 'on', 'SenseZ', 'on');
    end
    % tap converters: reuse c4 pattern - add temp converters + ToWorkspace
    taps = {};
    for pfx = {'EXT', 'FLX'}
        for c = 1:3
            cn = [pfx{1} '_t' num2str(c)];
            for suf = {'', 'w'}
                if getSimulinkBlockHandle([mdl '/' cn suf{1}]) > 0
                    delete_block([mdl '/' cn suf{1}]);
                end
            end
            if true
                add_block('nesl_utility/PS-Simulink Converter', [mdl '/' cn], 'Position', [260 400+40*c+100*strcmp(pfx{1},'FLX') 310 430+40*c+100*strcmp(pfx{1},'FLX')]);
                tsph = get_param([mdl '/' pfx{1} '_TS'], 'PortHandles');
                tsP = [tsph.LConn tsph.RConn];
                cph = get_param([mdl '/' cn], 'PortHandles');
                try
                    add_line(mdl, tsP(2+c), cph.LConn(1));
                catch
                    fprintf('tap %s already wired\n', cn);
                end
                add_block('simulink/Sinks/To Workspace', [mdl '/' cn 'w'], 'Position', [340 400+40*c+100*strcmp(pfx{1},'FLX') 400 430+40*c+100*strcmp(pfx{1},'FLX')]);
                set_param([mdl '/' cn 'w'], 'VariableName', ['t_' cn], 'SaveFormat', 'Structure With Time');
                cph2 = get_param(getSimulinkBlockHandle([mdl '/' cn]), 'PortHandles');
                wph2 = get_param(getSimulinkBlockHandle([mdl '/' cn 'w']), 'PortHandles');
                fprintf('conv ports: out=%d lconn=%d | w in=%d\n', numel(cph2.Outport), numel(cph2.LConn), numel(wph2.Inport));
                try
                    add_line(mdl, cph2.Outport(1), wph2.Inport(1));
                catch ME
                    % Simulink sometimes throws AFTER connecting: verify
                    ln2 = get_param(cph2.Outport(1), 'Line');
                    if ~(isscalar(ln2) && ln2 > 0)
                        rethrow(ME);
                    end
                    fprintf('wire threw but connected OK\n');
                end
            end
            taps{end+1} = ['t_' cn]; %#ok<SAGROW>
        end
    end
    try
        set_param(mdl, 'StopTime', '0.005');
        out = sim(mdl);
    catch ME
        fprintf('sim failed: %s\n', ME.message);
        break;
    end
    maxErr = 0;
    for pfx = {'EXT', 'FLX'}
        r = zeros(1, 3);
        for c = 1:3
            s = out.(['t_' pfx{1} '_t' num2str(c)]);
            r(c) = s.signals(1).values(1);
        end
        tgt = target.(pfx{1});
        err = r - tgt;
        maxErr = max(maxErr, norm(err));
        fprintf('pass %d %s: measured [%.4f %.4f %.4f] target [%.4f %.4f %.4f] err %.4f\n', ...
            pass, pfx{1}, r, tgt, norm(err));
        % split the correction between origin and insertion transforms
        dO = str2num(get_param([mdl '/' pfx{1} '_origT'], 'TranslationCartesianOffset')); %#ok<ST2NM>
        dI = str2num(get_param([mdl '/' pfx{1} '_insT'], 'TranslationCartesianOffset')); %#ok<ST2NM>
        half = err / 2;
        % measured vector is expressed in the ORIGIN frame axes; the origin
        % transform shifts the origin, the insertion transform shifts the
        % insertion -> correcting both by -/+ half err closes the gap
        set_param([mdl '/' pfx{1} '_origT'], 'TranslationCartesianOffset', mat2str(dO + half, 8));
        set_param([mdl '/' pfx{1} '_insT'], 'TranslationCartesianOffset', mat2str(dI - half, 8));
    end
    if maxErr < 2e-3
        fprintf('converged after %d pass(es)\n', pass);
        break;
    end
end
% restore: turn xyz sensing back off, remove temp taps
for pfx = {'EXT', 'FLX'}
    set_param([mdl '/' pfx{1} '_TS'], 'SenseX', 'off', 'SenseY', 'off', 'SenseZ', 'off', 'SenseDist', 'on');
    for c = 1:3
        for suf = {'', 'w'}
            try, delete_block([mdl '/' pfx{1} '_t' num2str(c) suf{1}]); catch, end
        end
    end
end
try
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('UPDATE OK\n');
catch ME
    fprintf('UPDATE FAILED: %s\n', ME.message);
end
save_system(mdl);
close_system(mdl, 0);
fprintf('=== calibrate_muscles DONE ===\n');
end
