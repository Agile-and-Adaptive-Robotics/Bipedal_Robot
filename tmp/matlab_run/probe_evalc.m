% Probe: does a local function's evalc-assigned variable return as output?
function probe_evalc()
    v = getViaEvalc();
    assert(isa(v, 'double') && v == 42, 'evalc output binding failed');
    disp('EVALC PROBE PASS');
end

function out = getViaEvalc()
    evalc("out = 6*7;");
end
