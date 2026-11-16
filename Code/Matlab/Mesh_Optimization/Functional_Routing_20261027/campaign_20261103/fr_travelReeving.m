function r = fr_travelReeving(gain)
% Ideal parallel-fall topology. Guide and cam hardware are separate counts.
switch gain
    case 1
        nodes = {'moving_input', 'output'};
        moving = [1 0];
    case 2
        nodes = {'fixed_dead_end', 'moving_sheave', 'output'};
        moving = [0 1 0];
    case 3
        nodes = {'moving_dead_end', 'fixed_return_sheave', 'moving_sheave', 'output'};
        moving = [1 0 1 0];
    otherwise
        error('fr_travelReeving:Gain', 'No defined reeving for gain %g', gain);
end
% An output stroke must absorb the length changes of all moving falls.
falls = sum(abs(diff(moving)));
assert(falls == gain);
r = struct('nodes', {nodes}, 'sheaveCount', sum(contains(nodes, 'sheave')), ...
    'movingFalls', falls, 'forceRatio', 1/falls, ...
    'physicalStatus', 'ideal topology only; sheave centers and cable tangencies unverified');
if gain == 1
    r.physicalStatus = 'direct transmission; no tackle sheaves';
end
end
