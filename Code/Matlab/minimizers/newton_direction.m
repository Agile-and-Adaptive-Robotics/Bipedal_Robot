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
    %%
    

    H = zeros(n,n);

    %compute centered difference gradient along with forwrd difference Hessian
    xs = zeros(n,1);

    %Initialize our counter
    index = 1;

    for i=1:n
        for j=i:n
            x_add1 = zeros(n,1);
            x_add1(i) = dx(i);
            x_add2 = zeros(n,1);
            x_add2(j) = dx(j);

            %
            if i == j
                %We're on the diagonal, so we can use these function
                %evaluations to compute the gradient as well  
                xs(:,index) = x_curr-x_add1;
                index = index + 1;
                xs(:,index) = x_curr+x_add2;
                index = index + 1; 

            else
                %Otherwise, just compute that element in the Hessian once,
                %and mirror it across the diagonal
                xs(:,index) = x_curr-x_add1+x_add2;
                index = index + 1;
                xs(:,index) = x_curr-x_add1;
                index = index + 1;
                xs(:,index) = x_curr+x_add2;
                index = index + 1;
            end
        end
    end  

    if f_batch
        %Call f once with all of the x values and barrier terms
        fs = f(xs);
    else
        %Preallocate space for f, then call f repeatedly, for
        %each x and barrier.
        fs = zeros(1,size(xs,2));
        for i=1:size(xs,2)
            fs(i) = f(xs(:,i));
        end
    end

    index = 1;
    for i=1:n
        for j=i:n
            if i == j
                f_neg = fs(:,index);
                index = index + 1;
                f_pos = fs(:,index);
                index = index + 1;

                H(i,j) = (f_pos - 2*f_curr + f_neg)/(x_add1'*x_add1);

                f_calls = f_calls + 2;
            else
                %Otherwise, just compute that element in the Hessian once,
                %and mirror it across the diagonal
                f1 = fs(:,index);
                index = index + 1;
                f2 = fs(:,index);
                index = index + 1;
                f3 = fs(:,index);
                index = index + 1;

                H(i,j,:) = ((f1 - f2) - (f3 - f_curr))/(-dx(i)*dx(j));

                H(j,i) = H(i,j);

                f_calls = f_calls + 3;
            end
        end
    end

    try chol(H);
    catch
        warning('Hessian is not PD.')
        H = eye(n);
    end
%     H_inv = H^-1;
    
%     clc
%     g
%     g_anal = [-400*x_curr(1)*(x_curr(2)-x_curr(1)^2)-2*(1-x_curr(1));...
%         200*(x_curr(2)-x_curr(1)^2)]
%     
%     H
%     H_anal = [-400*(x_curr(2)-x_curr(1)^2)+800*x_curr(1)^2+2,-400*x_curr(1);...
%         -400*x_curr(1),200]
%         
%     H_inv
%     1/det(H)*[H(2,2),-H(1,2);-H(2,1),H(1,1)]
%     H_inv = H^-1
%     H_anal^-1
%     p_anal = -H_anal^-1*g_anal
%     p = -H_inv*g
    H_inv = H^-1;
    p = linsolve(-H,g);
    
    s = norm(p);
    p = p/s;
end