function [ x_final, solution_found ] = bisect( f, x_lo, x_hi, x_tol, f_tol, max_it )
%BISECT finds the root of a one dimensional function within the bounds
%given
%   Detailed explanation goes here
    
    x_hi_orig = x_hi;
    x_lo_orig = x_lo;
    x_range = x_hi_orig - x_lo_orig;
    
    nan_tries = 1;
    f_lo = NaN;
    while isnan(f_lo) && nan_tries <= 2
        f_lo = f(x_lo);
        if isnan(f_lo)
            warning('f_lo evaluated to NaN.')
        end
        nan_tries = nan_tries + 1;
    end
    if isnan(f_lo)
        error('Bisect method stopped because function evaluated to NaN.')
    end
    
    nan_tries = 1;
    f_hi = NaN;
    while isnan(f_hi) && nan_tries <= 2
        f_hi = f(x_hi);
        if isnan(f_hi)
            warning('f_hi evaluated to NaN.')
        end
        nan_tries = nan_tries + 1;
    end
    if isnan(f_hi)
        error('Bisect method stopped because function evaluated to NaN.')
    end
    
    i = 1;
    solution_found = 0;
    converged = 0;
    while ~converged
%         fprintf('%i, %f \n',i,x_range)
        x_mid = x_lo + .5*x_range;
        
        nan_tries = 1;
        f_mid = NaN;
        while isnan(f_mid) && nan_tries <= 2
            f_mid = f(x_mid);
            if isnan(f_mid)
                warning('f_hi evaluated to NaN.')
            end
            nan_tries = nan_tries + 1;
        end
        if isnan(f_mid)
            error('Bisect method stopped because function evaluated to NaN.')
        end
        
        if f_hi == 0 || f_lo == 0
%             disp('At least one is 0')
            if f_hi == 0 && f_lo ~= 0
                x_mid = x_hi;
                converged = 1;
            elseif f_lo == 0 && f_hi ~= 0
                x_mid = x_lo;
                converged = 1;
            else
                %Both x_hi and x_lo give solutions
%                 disp('Multiple solutions found.')
                x_hi = x_hi + .01*x_range;
                x_lo = x_lo - .01*x_range;
                f_hi = f(x_hi);
                f_lo = f(x_lo);
            end       
        else
            if sign(f_lo*f_mid) == 1
                f_lo = f_mid;
                x_lo = x_mid;
            elseif sign(f_mid*f_hi) == 1
                f_hi = f_mid;
                x_hi = x_mid;
            else
%                 fprintf('%f, %f, %f. what?? \n',f_lo, f_hi, f_mid)
                %We may have 2 solutions, or 0.
                if sign(f_lo*f_mid) == -1 && sign(f_mid*f_hi) == -1
                    if abs(f_lo) < abs(f_hi)
                        %x_lo is closer to a root
                        f_hi = f_mid;
                        x_hi = x_mid;
                    elseif abs(f_hi) <= abs(f_lo)
                        %x_hi is closer to a root
                        f_lo = f_mid;
                        x_lo = x_mid;
                    end
                else
                    %No solution may exist
%                     disp('No solution may exist.')
                    converged = 1;
                end
            end
        end
        
        x_range = x_hi - x_lo;
            
        i = i + 1;
        if i > max_it
            converged = 1;
            x_mid = (x_lo + x_hi)/2;
        end
        if norm(x_hi-x_lo) < x_tol
            converged = 1;
            x_mid = (x_lo + x_hi)/2;
        end
        if norm(f_hi-f_lo) < f_tol
            converged = 1;
            x_mid = (x_lo + x_hi)/2;
        end
        
    end
    x_final = x_mid;
    
    %If we are really close to the high edge, then we are probably seeking a
    %point that is not within our bounds.
    if abs(x_final - x_hi_orig) > x_tol
        solution_found = true;
    else
        solution_found = false;
    end
    
end

