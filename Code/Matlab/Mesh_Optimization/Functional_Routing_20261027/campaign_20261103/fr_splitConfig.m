function cfg = fr_splitConfig(~, cs, x)
%FR_SPLITCONFIG The recorded transmission configuration for a candidate
% (fr_mechBuild's section bookkeeping without a full mechanics build).
design = fr_unpackX(cs, x);
cfg = struct();
cfg.seriesSections = 1 + double(design.seriesSplit < 0.999);
cfg.sectionLengths_m = [design.L0 * design.seriesSplit, ...
    design.L0 * (1 - design.seriesSplit)];
end
