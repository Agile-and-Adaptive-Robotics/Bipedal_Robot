function dx = calculate_dx( x_typical, bnds, coarse_dx )
%CALCULATE_DX Use the current solution for the problem to adjust the finite
%difference value, dx.
    n = length(x_typical);
    zero_ind = find(x_typical == 0);
    nonzero_ind = 1:n;
    nonzero_ind(zero_ind) = [];
    
    %if at least some values in x_typical are not 0, we can use them to
    %estimate a typical x value (the mean of the nonzero values).
    if isempty(bnds)
        if isempty(nonzero_ind) %If they are all zeros
            x_typical = ones(n,1);
        else %If some are zeros
            x_typical(zero_ind) = mean(x_typical(nonzero_ind));
        end
    else
        x_typical(zero_ind) = abs(bnds(zero_ind,2));
    end
    
    %Rule of thumb for minimization problems.
    if coarse_dx
        dx = eps^(1/6)*x_typical; 
    else
        dx = eps^(1/2)*x_typical; 
    end
end