function output = objective_function(f,x_map,T_inv,x,weight)
    f_output = f(x_map(T_inv*x));
    
    f_length = length(f_output);
    if f_length == 0
        warning('Objective function returns an empty vector.')
        output = NaN;
    elseif f_length == 1
        output = f_output;
    else
        output = f_output(1) + weight*sum(f_output(2:end));
    end
    if length(output) > 1
        error('Objective function does not return a scalar.')
    end
end