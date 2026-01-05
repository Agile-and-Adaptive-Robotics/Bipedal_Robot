function [ x_final, solution_found ] = bisect_min( f, dx, x_lo, x_hi, x_tol, f_tol, max_it )
%BISECT_MIN finds a local minimum of a one dimensional function within the bounds
%given
%   Detailed explanation goes here
    
%     dx = 1e-3*(x_hi - x_lo);
    df = @(x) (f(x+dx) - f(x-dx))/(2*dx);
    [x_final, solution_found] = bisect(df,x_lo,x_hi,x_tol,f_tol,max_it);
    
    f_lo = f(x_final - dx);
    f_mid = f(x_final);
    f_hi = f(x_final + dx);
    
    concave_up = ((f_lo - 2*f_mid + f_hi)/(dx^2) > 0);
    
    if ~concave_up
        %we found a maximum, not a minimum. Repeat with -df.
        warning('A local maximum was found. Repeating bisection to find local minimum.')
        neg_df = @(x) -df(x);
        [x_final, solution_found] = bisect(neg_df,x_lo,x_hi,x_tol,f_tol,max_it);
    end
      
end

