function [x_final,s,f_curr,f_calls] = bfgs_search(f,x_0,f_0,p,s,g,options)
%BFGS_SEARCH Perform a line search in the bfgs direction
    
    if norm(s) == 0 || isnan(norm(s))
        %The direction step may specify a very small step size. In that
        %case, skip the search.
        x_final = x_0;
        f_curr = f_0;
        s = 0;
        f_calls = 0;
    else
        to_backtrack = options(1);
        alpha = options(2);
        backtrack_ratio = options(3);
        max_it_search = options(4);
        
        if to_backtrack
            converged = false;
            i = 1;
            while ~converged && i <= max_it_search
                %Check alpha condition
                f_curr = f(x_0+s*p);
                if f_curr <= f_0+alpha*s*g'*p;
                    %Success; assign values and exit
                    converged = true;
                    x_final = x_0 + s*p;
                    f_calls = 1 + i;
                else
                    if norm(s) < eps
                        %alpha is not satisfied, but our step size is
                        %practically zero. This suggests that we are at the
                        %minimum already. Return NaN, as a code to the main
                        %routine that we are already at the mininum.
                        converged = true;
                        x_final = x_0;
                        s = NaN;
                        f_curr = NaN;
                        f_calls = 1 + i;
                    else
                        %alpha is not satisfied, and we can improve our
                        %position.
                        s = s * backtrack_ratio;
                        i = i + 1;
                    end
                end     
            end
            
            if i >= max_it_search
                x_final = NaN;
                s = NaN;
                f_curr = NaN;
                f_calls = 1 + i;
            end
            
        else %fit
            %Calculate our initial lambda (1)
            lambda_curr = 1;
            %Calculate the objective function at lambda_1
            f_curr = f(x_0+lambda_curr*s*p);
            %Calculate the directional derivative at lambda=0
            f_prime_0 = g'*p;
            
            %Test the alpha condition
            if f_curr <= f_0+alpha*s*g'*p 
                %If we pass, good, end the function
                x_final = x_0+s*p;
                f_calls = 1;
            else
                %Try a quadratic fit first
                lambda_prev = 1;
                
                %Use our experssion for quadratic fit to find our next
                %lambda to test
                lambda_curr = -f_prime_0/(2*(f_curr-f_0-f_prime_0));
                %Make sure lambda_curr is a reasonable size. Too large or
                %small defeats the purpose.
                lambda_curr = max(.1*lambda_prev,min(lambda_curr,0.5*lambda_prev));
                %Limite lambda based on alpha
                lambda_curr = min(lambda_curr, 1/(2*(1-alpha)));
                
                f_curr = f(x_0 + lambda_curr*p);
                if f_curr <= f_0+alpha*lambda_curr*g'*p
                    %alpha condition is satisfied
                    s = lambda_curr;
                    x_final = x_0 + lambda_curr*p;
                    f_calls = 2;
                else
                    %alpha condition is not satisfied; try a cubic fit
                    %(final try).
                    %Shift our f values
                    %f_prev_prev = f_prev;
                    f_prev = f_curr;
                    
                    %Shift our lambdas
                    lambda_prev_prev = lambda_prev;
                    lambda_prev = lambda_curr;
                    
                    %We will continue to discard previous guesses and
                    %improve our cubic fit until we converge or try too
                    %many times.
                    converged = false;
                    i = 1;
                    while ~converged && i <= max_it_search
                        %Solve the 2x2 system of equations to find a
                        %and b, the coefficients of the lambda^3 and
                        %lambda^2 terms, respectively. This is
                        %inverting a 2x2 matrix and multiplying it by
                        %another. The inverse has been solved
                        %analytically to improve performance.
                        temp = 1/(lambda_prev^3*lambda_prev_prev^2-...
                            lambda_prev_prev^3*lambda_prev^2)*...
                            [lambda_prev_prev^2 -lambda_prev^2; ...
                            -lambda_prev_prev^3 lambda_prev^3]*...
                            [f_prev - lambda_prev*f_prime_0-f_0; ...
                            f_prev - lambda_prev_prev*f_0 - f_0];
                        a = temp(1);
                        b = temp(2);
                    
                        %We can use the quadratic formula to find the
                        %zeros of the derivative. This is the lower
                        %guess. 
                        lambda_1 = (-2*b - sqrt((2*b)^2-4*a*f_prime_0))/(6*a);
                        
                        %Shift our lambdas
                        lambda_prev_prev = lambda_prev;
                        lambda_prev = lambda_curr;
                        
                        %Choose which lambda to use
                        if 6*a*lambda_1 + 2*b > 0
                            %This is the minimizer because M''>0
                            lambda_curr = lambda_1;
                        else
                            %Otherwise, the larger lambda must be the
                            %mimizer
                            lambda_curr = (-2*b + sqrt((2*b)^2-4*a*f_prime_0))/(6*a);
                        end
                        
                        %Limit lambda to be within reasonable ranges
                        lambda_curr = max(lambda_curr, lambda_prev/10);
                        lambda_curr = min(lambda_curr, lambda_prev/2);
                        
                        %Shift our functional values
                        %f_prev_prev = f_prev;
                        f_prev = f_curr;
                        
                        %Evaluate the function at the new lambda
                        f_curr = f(x_0 + lambda_curr*p);
                        
                        %Check the alpha condition
                        if f_curr <= f_0+alpha*lambda_curr*g'*p
                            %alpha satisfied
                            converged = true;
                            s = lambda_curr;
                            x_final = x_0 + lambda_curr*p;
                            f_calls = 2 + i;
                        end
                        
                        if lambda_curr < eps
                            %If lambda gets really small, that
                            %suggests that we are at an optimal
                            %point; no step, no matter how small,
                            %can decrease the function.
                            converged = true;
                            s = lambda_curr;
                            x_final = x_0;
                            f_calls = 2 + i;
                        end
                        i = i + 1;
                    end
                        
                    %If we run out of iterations, return NaN so we know to 
                    %try trust region.
                    if i >= max_it_search
                        x_final = NaN;
                        s = NaN;
                        f_curr = NaN;
                        f_calls = max_it_search + 2;
                    end
                end
            end
        end
    end
end
