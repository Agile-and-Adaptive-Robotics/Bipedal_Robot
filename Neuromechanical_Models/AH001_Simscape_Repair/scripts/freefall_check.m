% freefall_check.m  (AH001_Simscape_Repair, 2026-11-01)
% 0.5 s free-fall validation of the repaired closed-loop model:
%  - integrates cleanly (no solver failure / NaN)
%  - SixDOF base Py decreases ~ -1.229 m (y-up gravity check)
%  - knee four-bar pins respond smoothly (closed-loop consistency)

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'freefall_check_20261101.txt');
fid = fopen(logFile, 'w');

mdl = 'ah001_repair';
load_system(fullfile(repairRoot, 'working', [mdl '.slx']));

try
    out = sim(mdl, 'StopTime', '0.5', 'SimscapeLogType', 'all', ...
              'ZeroCrossControl', 'UseLocalSettings');
    p(fid, 'sim(0.5 s, free fall): OK\n');
    sl = out.simlog;

    try
        ts = sl.SixDOF.Py.p.series;
        v = ts.values('m');
        p(fid, 'SixDOF Py: %g -> %g m (delta %g m; ideal free fall 0.5 s = -1.229 m)\n', ...
            v(1), v(end), v(end) - v(1));
    catch ME
        p(fid, 'SixDOF Py read: %s\n', ME.message);
    end

    pins = struct('name', {'Revolute3', 'Revolute5'}, ...
                  'pair',  {'KT_L-FL_2', 'KT_R-FL_1'});
    for k = 1:numel(pins)
        try
            ts = sl.x09_BA_001_1.(pins(k).name).Rz.q.series;
            qv = ts.values('deg');
            p(fid, 'pin %-9s (%s): %g -> %g deg, range %g deg\n', ...
                pins(k).name, pins(k).pair, qv(1), qv(end), max(qv) - min(qv));
        catch ME
            p(fid, 'pin %s read fail: %s\n', pins(k).name, ME.message);
        end
    end
catch ME
    p(fid, 'sim FAIL: %s\n', getReport(ME, 'basic', 'hyperlinks', 'off'));
end

close_system(mdl, 0);
p(fid, 'FREEFALL CHECK DONE.\n');
fclose(fid);

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
