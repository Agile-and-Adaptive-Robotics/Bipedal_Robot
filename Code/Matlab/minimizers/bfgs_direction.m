function [p,s,g,H,H_inv,f_calls] = bfgs_direction(f,x_curr,x_prev,f_curr,g_prev,H_prev,H_inv_prev,dx,f_batch)
%BFGS_DIRECTION compute the Hessian and gradient for a quasi-newton
%optimization routine
    
    n = length(x_curr);
    %assign basic values
    s_curr = x_curr - x_prev;
    
    %compute gradient, return number of function calls
    g = NaN(n,1);
    f_calls = 0;
    
    %Generate all of the inputs we will need to evaluate to calculate the
    %gradient.
    xs = NaN(n,2*n);
    fs = NaN(2*n,1);
    for i=1:n
        x_for = x_curr;
        x_bac = x_curr;

        x_for(i) = x_for(i) + dx(i);
        x_bac(i) = x_bac(i) - dx(i);

        xs(:,2*i-1) = x_for;
        xs(:,2*i) = x_bac;
    end
    
    if f_batch
        %Calculate the centered finite difference derivative
        %Initialize a matrix in which each column is a point we will need
        %to evaluate to compute the gradient.
        fs = f(xs);
        f_calls = f_calls + 2*n;
        
    else
        %If we evaluate these one by one, then we must use a loop.
        for i=1:n
            f_for = f(xs(:,2*i-1));
            f_bac = f(xs(:,2*i));

            fs(2*i-1) = f_for;
            fs(2*i) = f_bac;
            
            f_calls = f_calls + 2;
        end
    end
    
    %Take the function evaluations and compute the gradient along each
    %direction. Make sure that if infinite or NaN slopes are calculated,
    %those are reevaluated as forward or backward differences.
    for i=1:n
        f_for = fs(2*i-1);
        f_bac = fs(2*i);
        grad_try = (f_for - f_bac)/(2*dx(i));

        if isinf(grad_try) || isnan(grad_try)
            grad_try_for = (f_for - f_curr)/dx(i);
            grad_try_bac = (f_curr - f_bac)/dx(i);

            if ~isinf(grad_try_for) && ~isnan(grad_try_for)
                g(i) = grad_try_for;
            elseif ~isinf(grad_try_bac) && ~isnan(grad_try_bac)
                g(i) = grad_try_bac;
            else
                g(i) = 0;
            end
        else
            g(i) = grad_try;
        end        
    end    
    
    %Definitions for the Broyden updates
    %g is grad_curr
    y_curr = g - g_prev;
    
    do_not_update = false;
    %Broyden updates
    if s_curr'*y_curr > 0
        H_inv = H_inv_prev + (s_curr'*y_curr + y_curr'*H_inv_prev*y_curr)*...
            (s_curr*s_curr')/(s_curr'*y_curr)^2 - ...
            (H_inv_prev*y_curr*s_curr'+s_curr*y_curr'*H_inv_prev)/...
            (s_curr'*y_curr);
        
        H = H_prev + (y_curr*y_curr')/(y_curr'*s_curr) - ...
            (H_prev*(s_curr*s_curr')*H_prev)/(s_curr'*H_prev*s_curr);
        
        try chol(H);
        catch
            warning('The Broyden update for the Hessian is not positive definite.')
            do_not_update = true;
        end
    else
        do_not_update = true;
    end
    
    if do_not_update
        %If we are not to update the hessian, use the previous value. For
        %the first step, this will be identity, resulting in a steepest
        %descend step.
        H_inv = H_inv_prev;
        H = H_prev;
    end
    
%     g = [-400*x_curr(1)*(x_curr(2)-x_curr(1)^2)-2*(1-x_curr(1));...
%         200*(x_curr(2)-x_curr(1)^2)]
%     
% 
%     H = [-400*(x_curr(2)-x_curr(1)^2)+800*x_curr(1)^2+2,-400*x_curr(1);...
%         -400*x_curr(1),200]
%         
%     H_inv = H^-1
    p = -H_inv*g;
    
      
%     keyboard
    s = norm(p);
    p = p/s;
end