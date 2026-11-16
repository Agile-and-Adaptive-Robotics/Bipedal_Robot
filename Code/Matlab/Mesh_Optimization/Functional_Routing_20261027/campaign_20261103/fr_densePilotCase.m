function cs = fr_densePilotCase(C,kin,cs)
% Independent pose grid. Interpolate the existing human target samples.
for p = 1:numel(cs.passes)
    pass = cs.passes{p};
    old = pass.th1(:);
    v = linspace(min(old),max(old),31);
    if strcmp(pass.sweptCoord,'knee_angle_r')
        v = linspace(max(min(old),deg2rad(-120)),min(max(old),0),31);
    end
    if numel(pass.th2)==numel(old)
        pass.th2 = interp1(old,pass.th2(:),v,'linear');
    else
        pass.th2 = repmat(pass.th2(1),size(v));
    end
    pass.th1 = v;
    pass.qStructs = cell(size(v));
    for j = 1:numel(v)
        q = kin.qDefault;
        q.(pass.sweptCoord) = v(j);
        if ~isempty(pass.sweptCoord2)
            q.(pass.sweptCoord2) = pass.th2(j);
        end
        pass.qStructs{j} = q;
    end
    cs.passes{p} = pass;
end
cs.target = fr_targets(C,kin,cs);
end
