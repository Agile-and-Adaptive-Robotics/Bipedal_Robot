function TestMonoPam_multi
% Compatibility entry point for the earlier test name.
% MonoPam_multi.m was the discarded wrapper. The maintained standalone
% class and test use the requested _mult spelling.
warning('TestMonoPam_multi:Renamed', ...
    'Running TestMonoPam_mult. MonoPam_multi.m is obsolete.')
TestMonoPam_mult
end
