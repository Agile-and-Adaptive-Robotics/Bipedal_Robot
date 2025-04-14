function Fmn = f_festo(rel, P, D)
% Returns normalized muscle force as a function of relative contraction and pressure
% Inputs:
%   rel - relative strain (scalar or array)
%   P   - normalized pressure (scalar)
%   D   - BPA diameter (10, 20, 40)

    switch D
        case 10
            a0 = 0.5822;
            a1 = 4.142;
            a3 = 0.5368;
        case 20
            a0 = 0.2607;
            a1 = 6.398;
            a3 = 1.303;
        case 40
            a0 = 0.1224;
            a1 = 10.47;
            a3 = 2.023;
        otherwise
            error('Unsupported diameter D: %d. Use 10, 20, or 40.', D);
    end

    Fmn = a0 * (exp(-a1 .* rel) - 1) + P .* exp(-a3 .* (rel.^2));
end

