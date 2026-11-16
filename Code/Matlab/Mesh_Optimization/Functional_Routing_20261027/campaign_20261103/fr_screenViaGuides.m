function [bestCase, bestF] = fr_screenViaGuides(C, cs, kin, x, branch, opts, bestF)
% Screen one new body-fixed guide at a time around existing sleeve routes.
% The no-guide solution remains in the comparison.
bestCase = cs;
clouds = struct('femur_r', 'Femur', 'tibia_r', 'Tibia', 'pelvis', 'Pelvis_R');
for r = 1:numel(cs.rows)-1
    body = cs.rowOwning{r};
    if ~strcmp(body, cs.rowOwning{r+1}) || ~isfield(clouds, body)
        continue
    end
    sharedFrame = true;
    for p = 1:numel(cs.passes)
        if all(cs.rowUse{p}(r:r+1))
            sharedFrame = sharedFrame && strcmp(cs.rowFrame{p}{r}, cs.rowFrame{p}{r+1});
        end
    end
    if ~sharedFrame
        continue
    end
    midpoint = mean(cs.seedLoc(r:r+1,:), 1);
    patch = fr_patch(C, clouds.(body), body, midpoint, 30);
    offsets = [0 0; -.75 0; .75 0; 0 -.75; 0 .75];
    for j = 1:size(offsets,1)
        point = patch.project(offsets(j,1), offsets(j,2), C.standoff_mm(2));
        candidate = fr_addViaGuide(cs, r, point(:)');
        try
            f = fr_objective(C, candidate, kin, x, branch, opts);
        catch ME
            fprintf('Guide rejected for %s row %d: %s\n', cs.name, r, ME.message);
            continue
        end
        if isfinite(f) && f < bestF
            bestCase = candidate;
            bestF = f;
        end
    end
end
end
