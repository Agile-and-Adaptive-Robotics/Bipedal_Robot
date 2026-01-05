function [x,p,f_try,g,H,H_inv,full_newton_step,improved_point_found,delta,f_calls] = bfgs_trust_region(f,x_0,x_prev,f_0,g_prev,H_prev,H_inv_prev,dx,f_batch,options)
    alpha = options(1);
    eta = options(2);
    delta = options(3);
    need_updated_gradient = options(4);
    full_newton_step = false;
    improved_point_found = false;

    %Compute the gradient and Hessian at this point.
    if need_updated_gradient
        [~,~,g,H,H_inv,f_calls] = bfgs_direction(f,x_0,x_prev,f_0,g_prev,H_prev,H_inv_prev,dx,f_batch);
    else
        g = g_prev;
        H = H_prev;
        H_inv = H_inv_prev;
        f_calls = 0;
    end
    %double dogleg
    lambda_cp = norm(g,2).^2/(g'*H*g);
    if lambda_cp < 0
        warning('Cauchy step is incorrectly negative. This may be caused by a discontinuous objective function.')
    end
    
    i = 1;
    converged = false;
    
    while ~converged
        if delta < lambda_cp * norm(g,2)
            %If our leash is tighter than the Cauchy step (steepest
            %descent), we should limit our travel to the leash. Return.
            p = delta*-g/norm(g,2); %cauchy direction, leash length
            x = x_0 + p;
        else
            %Otherwise, we should compute the next leg of the path.
            %Calculate the full cauchy step.
            p_cauchy = lambda_cp * -g; %cauchy direction and length
            
            %Calculate the full Newton step
            p_newton = H_inv_prev*-g;
            
            %Calculate the point N_hat, which is along the Newton
            %direction. The double dogleg runs from x_0 to x_cauchy to
            %N_hat to x_newton.
            N_hat = x_0 + eta*p_newton;
            if delta < eta*norm(p_newton)
                %Use the leash! Limited Newton step is still too far. Now
                %we must find where the x_cauchy to N_hat vector intersects
                %with our trust region boundary (radius delta).
                r_inter = N_hat - p_cauchy;
                %To calculate how far along the first dogleg to travel, we
                %use the law of cosines. We know the hypotenuse has length
                %delta, another leg is p_cauchy, and the other leg is  
                %length lambda, which we must solve for, along r_inter.
                %delta^2 = p_cauchy'*p_cauchy + lambda^2*a'*a - 
                %2*|s_cp|*|lambda*a|*cos(beta). Beta is the angle opposite
                %delta, which is pi minus the angle between p_cauchy and
                %r_inter. Solving for lambda produces two solutions; we
                %want the positive one, which requires we add the sqrt term
                %.                
                a = 1;
                b = 2*(p_cauchy'*r_inter)/norm(r_inter);
                c = p_cauchy'*p_cauchy - delta^2;
                
                lambda = (-b + sqrt(b^2 - 4*a*c))/(2*a);
                
                if ~isreal(lambda) && delta > eps
                    warning('The trust region includes infeasible points. Taking steepest descent direction and reducing trust region.')
                    delta = 0.9*(lambda_cp * norm(g,2));
                else
                    %do nothing
                end
                p = p_cauchy + lambda*r_inter;
                x = x_0 + p;
            else
                %We will be between N_hat and p_newton
                if delta < norm(p_newton)
                    %use the leash
                    p = delta*p_newton/(norm(p_newton));
                    x = x_0 + p;
                else
                    %Take the full Newton step
                    p = p_newton;
                    x = x_0 + p;
                    full_newton_step = true;
                end
            end  
        end
        %Evaluate the new point and test the alpha condition. But first, if
        %we took a trivial step, end the program.
        if norm(p) < eps
            converged = true;
            f_try = f_0;
        else
            f_try = f(x);
            alpha_satisfied = (f_try < f_0 + alpha*g'*(x - x_0));
            if alpha_satisfied
                %If alpha is satisfied, we can either stop or increase delta
                %and continue. First, compare actual performance to predicted
                %performance. If the objective is modeled well locally,
                %increase the leash length (delta).
                d_f = f_try - f_0;
                d_f_predicted = f_0 + g'*(x - x_0);

                if norm(p) > eps && norm(d_f_predicted - d_f) <= 0.1*norm(d_f) &&...
                        ~isinf(d_f) && ~isinf(delta)
                    %Model is accurate and we are taking nontrivla steps;
                    %increase the leash
                    delta = delta * 2;
                else
                    %Model is inaccurate, we should take what we have and exit.
                    converged = true;
                    f_calls = f_calls + i;
                end
                %We have found a point that satisfies the alpha condition, so
                %we should continue searching.
                improved_point_found = true;
            else
                %We failed the alpha condition, so shorten lambda.
                lambda = -g'*(x-x_0)/(2*(f_try - f_0 - g'*(x - x_0)));
                new_delta = delta * max(0.1,min(0.5,lambda*norm(x - x_0)));

                if new_delta < eps || new_delta == delta
                    %Trust region is too small
                    converged = true;
                else
                    delta = new_delta;
                end
            end
            i = i + 1;
    %         keyboard
        end
    end
end