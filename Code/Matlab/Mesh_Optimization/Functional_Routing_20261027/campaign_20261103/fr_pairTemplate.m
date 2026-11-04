function po = fr_pairTemplate()
%FR_PAIRTEMPLATE Fixed-schema per-group result (rigid + compliant outs).
po = struct('name', '', 'mapId', NaN, 'rigid', fr_emptyOut(), ...
    'compliant', fr_emptyOut());
end
