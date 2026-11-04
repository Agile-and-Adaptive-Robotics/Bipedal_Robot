function fr_probe27()
%FR_PROBE27 Structural probe of all 27 contract groups' seed routes.
% Prints, per group: seed muscle (argmax per_muscle_Lmax_m), raw declared
% row bodies + types, conditional ranges, and kin joint topology.
% Read-only; no case building.

C = fr_setup();
kin = fr_kinematics(C.osimStock);
tc = jsondecode(fileread(fullfile(C.targetsDir, 'target_contract.json')));

fprintf('=== kin joints ===\n');
jn = fieldnames(kin.joints);
for i = 1:numel(jn)
    J = kin.joints.(jn{i});
    coords = {};
    for a = 1:6
        if ~isempty(J.axes{a}.coords)
            coords = [coords, J.axes{a}.coords]; %#ok<AGROW>
        end
    end
    fprintf('  %-10s parent %s <- child %s | coords: %s\n', J.name, ...
        J.parentFrame, J.childFrame, strjoin(coords, ','));
end
fprintf('=== chain (child body -> parent body, joint) ===\n');
cn = fieldnames(kin.chain);
for i = 1:numel(cn)
    fprintf('  %-10s <- %-10s via %s\n', cn{i}, ...
        kin.chain.(cn{i}).parentBody, kin.chain.(cn{i}).jointName);
end

fprintf('\n=== 27 groups ===\n');
for g = 1:numel(tc.groups)
    grp = tc.groups(g);
    [mx, im] = max(structfun(@(x) x, grp.sizing_data.per_muscle_Lmax_m));
    seed = grp.stockMuscles{im};
    rt = fr_muscleroute(C.osimStock, seed);
    fprintf('\n[%2d] %-9s seed %s (Lmax %.4f m) primary %s nBPA %d dia %s\n', ...
        grp.map_id, grp.group, seed, mx, grp.primaryDof, ...
        grp.sizing_data.n20_ceil_MIF_over_Fmax20, ...
        ternary(grp.sizing_data.ten_mm_screen_pass, '10', '20'));
    for r = 1:numel(rt.points)
        P = rt.points{r};
        extra = '';
        if isfield(P, 'range') && ~isempty(P.range)
            extra = sprintf(' COND %s [%.1f %.1f]deg', P.coord, ...
                P.range(1) * 180 / pi, P.range(2) * 180 / pi);
        end
        if isfield(P, 'fnXYZ')
            extra = ' MOVING';
        end
        fprintf('     row %d %-22s body %-10s%s\n', r, P.name, P.body, extra);
    end
end
end

function s = ternary(b, a, c)
if b
    s = a;
else
    s = c;
end
end
