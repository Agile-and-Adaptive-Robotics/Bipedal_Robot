function cut_knee_welds_20260926()
% cut the leftover import welds KT<->KB (both legs), then compile
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));

% ports recorded in the surgery1 plain-line map:
%   KT_R port8 <-> KB_R port3 ; KT_L port4 <-> KB_L port6
cuts = { {{'x04_01_KT_R_001_1_RIGID', 8}, {'x04_02_KB_R_001_1_RIGID', 3}}; ...
         {{'x04_03_KT_L_001_1_RIGID', 4}, {'x04_04_KB_L_001_1_RIGID', 6}} };

for c = 1:2
    pa = cuts{c}{1}; pb = cuts{c}{2};
    phA = get_param([sub '/' pa{1}], 'PortHandles');
    portsA = [phA.RConn phA.LConn];
    phB = get_param([sub '/' pb{1}], 'PortHandles');
    portsB = [phB.RConn phB.LConn];
    pA = portsA(pa{2});
    pB = portsB(pb{2});
    % verify they share a line (possibly branched)
    lnA = get_param(pA, 'Line');
    lnB = get_param(pB, 'Line');
    fprintf('%s: lineA=%d lineB=%d\n', pa{1}, lnA, lnB);
    if lnA > 0
        % collect ALL port handles on this connection (branches incl.)
        phs = [];
        try, phs = [get_param(lnA, 'SrcPortHandle') get_param(lnA, 'DstPortHandle')]; catch, end
        br = -1;
        try, br = get_param(pA, 'Branch'); catch, end
        if isscalar(br) && br > 0
            try, phs = [phs get_param(br, 'BranchHandles')]; catch, end
        end
        touchesB = false;
        for q = 1:numel(phs)
            if isscalar(phs(q)) && phs(q) > 0 && phs(q) == pB
                touchesB = true;
            end
        end
        if touchesB
            delete_line(lnA);
            fprintf('cut %s <-> %s\n', pa{1}, pb{1});
        else
            % walk one more level: other ports on the same line
            for q = 1:numel(phs)
                h = phs(q);
                if isscalar(h) && h > 0 && h ~= pA
                    ln2 = get_param(h, 'Line');
                    if isscalar(ln2) && ln2 > 0
                        phs2 = [];
                        try, phs2 = [get_param(ln2, 'SrcPortHandle') get_param(ln2, 'DstPortHandle')]; catch, end
                        for r = 1:numel(phs2)
                            if isscalar(phs2(r)) && phs2(r) > 0 && phs2(r) == pB
                                delete_line(ln2);
                                fprintf('cut (branch) %s <-> %s\n', pa{1}, pb{1});
                                touchesB = true;
                                break;
                            end
                        end
                    end
                end
                if touchesB, break; end
            end
            if ~touchesB
                fprintf('WARNING: %s port%d does not reach %s port%d directly\n', pa{1}, pa{2}, pb{1}, pb{2});
            end
        end
    end
end

try
    set_param(mdl, 'StopTime', '3');
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('UPDATE OK\n');
    save_system(mdl);
    out = sim(mdl);
    LR = out.L_EXT_R.signals(1).values; LL = out.L_EXT_L.signals(1).values;
    fprintf('R knee EXT length: %.4f -> %.4f..%.4f m\n', LR(1), min(LR), max(LR));
    fprintf('L knee EXT length: %.4f -> %.4f..%.4f m\n', LL(1), min(LL), max(LL));
    save(fullfile(here, '..', 'results', 'humanoid_bpa_sns_20260926.mat'), 'out');
catch ME
    fprintf('UPDATE/SIM FAILED: %s\n', ME.message);
    for c = 1:min(numel(ME.cause), 8)
        fprintf('  CAUSE %d: %s\n', c, ME.cause{c}.message(1:min(end, 170)));
    end
end
save_system(mdl);
close_system(mdl, 0);
fprintf('=== cut_knee_welds DONE ===\n');
end
