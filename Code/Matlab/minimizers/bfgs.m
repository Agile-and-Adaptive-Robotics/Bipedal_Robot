function [x_final, f_final, stats] = bfgs(f, x0, bounds, linear_constraints,...
    penalty_weight, conv_criteria, direction_options, search_options, plot_options,...
    print_options)
%BFGS Broyden-Fletcher-Goldfard-Shanno function minimizer with boundaries
%and equality constraints.
%
%   INPUTS:
%   f - function to minimize.
%
%   x0 - initial point (n-element column vector).
%
%   bounds - limits on each variable. The user should input an empty
%   vector. If bounds are desired, use bnd_bfgs.m.
%
%   linear_constraints - 2 element cell array, {A,b}. Constraints are
%   specified by the equation Ax=b. A should be m x n, b should be m x 1,
%   in which m is the number of constraints. If there are none, an empty
%   vector or empty cell should be passed.
%
%   conv_criteria - up to 1 x 3 vector holding the maximum number of
%   iterations, the gradient tolerance for convergence, and the x stalling
%   tolerance for convergence. Default values are [1000,1e-6,1e-6].
%
%   direction_options - up to 1 x 2 vector that holds a
%   boolean of whether f is vectorized, and an initial finite difference
%   value (for the scaled problem)
%
%   search_options - up to 1 x 6 cell array that holds the string
%   'backtrack' or 'fit' to specify a line search method, alpha for the
%   alpha condition (1e-4), the backtrack ratio (0.5), the maximum number
%   of linesearch iterations (20), eta for the double dogleg (0.8), and
%   initial trust region radius (0.01).
%
%   plot_options - vector of what to plot (function values, x
%   ratios, if possible, x values on a contour of the function, norm(grad),
%   relative gradient (convergence criterion), relative x (convergence
%   criterion), step length, function calls by subroutine
%
%   print_options - 0 - Plot nothing
%                   1 - Plot results at the end of routine
%                   2 - Plot updates at each iteration 

    %Parse our input arrays to make sure the default values are applied if
    %information is missing.
    n_orig = length(x0);
    if isempty(linear_constraints)
       A = [];
       b = [];
    else
        %We have either empty vectors in a cell, or actual data
        all_empty = true;
        i = 1;
        while all_empty && i <= length(linear_constraints)
            if isempty(linear_constraints{i})
                %do nothing
            else
                all_empty = false;
            end
            i = i + 1;
        end
        
        if ~all_empty
            %If we have linear equality constraints, then we can eliminate some
            %variables, reducing the dimensionality of the problem. This is
            %important to all of the other problem setup we will do, so we must do
            %it first.
            A = linear_constraints{1};
            b = linear_constraints{2};
        else
            A = [];
            b = [];
        end
    end
    %Pass our constraints to the elim_var function. It will return a
    %function to map between the old and new x0, a new x0, and the indices
    %of the nonbasic (i.e. not eliminated) variables.
    allow_infeasible_start = false;
    [x_map,x0_nb,nb_var] = elim_var(A,b,x0,allow_infeasible_start);
    x0 = x0_nb;

    %How many variables do we have now?
    n = length(x0);

    %We may need to scale the variables. We can scale them to the same
    %order of magnitude, solve the problem, and then scale them back.
    %We need to do this now, because calculating our finite difference
    %requires that the problem be scaled.
    T = scale(x0);
    
    %Save the inverse to avoid computing it often. This is not expensive 
    %since T is diagonal.
    T_inv = T^-1;
    x0 = T*x0;
    x_typical = x0;
    
    %Make sure the convergence criteria input is the proper form.
    switch length(conv_criteria)
        case 0
            max_it = 1000;
            grad_tol = 1e-6;
            x_tol = 1e-6;
        case 1
            max_it = conv_criteria(1);
            grad_tol = 1e-6;
            x_tol = 1e-6;
        case 2
            max_it = conv_criteria(1);
            grad_tol = conv_criteria(2);
            x_tol = 1e-6;
        case 3
            max_it = conv_criteria(1);
            grad_tol = conv_criteria(2);
            x_tol = conv_criteria(3);
        otherwise
            max_it = 1000;
            grad_tol = 1e-6;
            x_tol = 1e-6;
            warning('Improper convergence criteria input; defaults are being used.')
    end    
    
    %The user may enter only the first input, which is a string. Test this
    %case, and if the input is a cell array, treat it like the others
    %above.
    if ischar(search_options)
        linesearch_string = search_options;
        alpha = 1e-4;
        backtrack_ratio = 0.5;
        max_linesearch_steps = 20;
        eta = 0.8;
        init_ddl_radius = 0.01;
    elseif isempty(search_options)
        linesearch_string = 'backtrack';
        alpha = 1e-4;
        backtrack_ratio = 0.5;
        max_linesearch_steps = 20;
        eta = 0.8;
        init_ddl_radius = 0.01;
    elseif iscell(search_options)
        switch length(search_options)
            case 0
                linesearch_string = 'backtrack';
                alpha = 1e-4;
                backtrack_ratio = 0.5;
                max_linesearch_steps = 20;
                eta = 0.8;
                init_ddl_radius = 0.01;
            case 1
                linesearch_string = search_options{1};
                alpha = 1e-4;
                backtrack_ratio = 0.5;
                max_linesearch_steps = 20;
                eta = 0.8;
                init_ddl_radius = 0.01;
            case 2
                linesearch_string = search_options{1};
                alpha = search_options{2};
                backtrack_ratio = 0.5;
                max_linesearch_steps = 20;
                eta = 0.8;
                init_ddl_radius = 0.01;
            case 3
                linesearch_string = search_options{1};
                alpha = search_options{2};
                backtrack_ratio = search_options{3};
                max_linesearch_steps = 20;
                eta = 0.8;
                init_ddl_radius = 0.01;
            case 4
                linesearch_string = search_options{1};
                alpha = search_options{2};
                backtrack_ratio = search_options{3};
                max_linesearch_steps = search_options{4};
                eta = 0.8;
                init_ddl_radius = 0.01;
            case 5
                linesearch_string = search_options{1};
                alpha = search_options{2};
                backtrack_ratio = search_options{3};
                max_linesearch_steps = search_options{4};
                eta = search_options{5};
                init_ddl_radius = 0.01;
            case 6
                linesearch_string = search_options{1};
                alpha = search_options{2};
                backtrack_ratio = search_options{3};
                max_linesearch_steps = search_options{4};
                eta = search_options{5};
                init_ddl_radius = search_options{6};
            otherwise
                linesearch_string = 'backtrack';
                alpha = 1e-4;
                backtrack_ratio = 0.5;
                max_linesearch_steps = 20;
                eta = 0.8;
                init_ddl_radius = 0.01;
                warning('Improper search options input; defaults are being used.')
        end
    else
        linesearch_string = 'backtrack';
        alpha = 1e-4;
        backtrack_ratio = 0.5;
        max_linesearch_steps = 20;
        eta = 0.8;
        init_ddl_radius = 0.01;
        warning('Improper search options input; defaults are being used.')
    end
    ddl_radius = init_ddl_radius;
    
    if print_options
        fprintf('Linesearch method: %s\nalpha: %f\nbacktrack ratio: %f\nbacktrack steps: %i\nddl eta: %f\ninitial ddl radius: %f\n\n',...
            linesearch_string,alpha,backtrack_ratio,max_linesearch_steps,eta,init_ddl_radius)
    end
    
    %Now, make sure that the linesearch method entered is correct.
    if strcmp(linesearch_string,'backtrack')
        to_backtrack = true;
        %to_fit = false;
    elseif strcmp(linesearch_string,'fit')
        %to_fit = true;
        to_backtrack = false;
    else
        to_backtrack = true;
        %to_fit = false;
        warning('Improper linesearch method input; backtracking will be used.')
    end
    
    if (isempty(bounds) && ismatrix(bounds)) || iscell(bounds)
        if (isempty(bounds) && ismatrix(bounds))
            bounds = {[],[]};
        end
        %input is probably ok
    else
        error('Bounds must be an empty array or a cell array with a boundary matrix and mu value (see bnd_bfgs.m)')
    end
    
    if isempty(bounds{1})
        bnds = [];
    else
        bnds = T*bounds{1}(nb_var,:);
    end
    
    %Make sure the direction options input is in the proper form.
    coarse_dx = false;
    switch length(direction_options)
        case 0
            f_batch = false;
            finite_diff = calculate_dx(abs(x0),bnds,false);
        case 1
            f_batch = direction_options(1);
            finite_diff = calculate_dx(abs(x0),bnds,false);
        case 2
            f_batch = direction_options(1);
            if length(direction_options(2)) == 1 && (direction_options(2) == 0 || direction_options(2) == 1)
                coarse_dx = true;
                finite_diff = calculate_dx(abs(x0),bnds,coarse_dx);
            else
                warning('Finite difference boolean is not formatted properly. Using fine sampling.')
                finite_diff = calculate_dx(abs(x0),bnds,coarse_dx);
            end
        otherwise
            f_batch = false;
            finite_diff = calculate_dx(abs(x0),bnds,coarse_dx);
            warning('Improper direction options input; defaults are being used.')
    end
    
    %If the penalty function shows that the starting point is infeasible,
    %we cannot run.
%     if ~isempty(bnds) && ( any(T_inv*x0 <= bnds(:,1)) || any(T_inv*x0 >= bnds(:,2)) )
%         infeasible = ((T_inv*x0 <= bnds(:,1)) | (T_inv*x0 >= bnds(:,2)));
    if ~isempty(bnds)
        if any(bnds(:,1) > bnds(:,2))
            warning('The BFGS bounds are not monotonic. They are being rearranged.');
            new_bnds = bnds;
            new_bnds(bnds(:,1) > bnds(:,2),1) = bnds(bnds(:,1) > bnds(:,2),2);
            new_bnds(bnds(:,1) > bnds(:,2),2) = bnds(bnds(:,1) > bnds(:,2),1);
            
            bnds = new_bnds;
            bounds{1} = T_inv*bnds;
        end
        if ( any(x0 <= bnds(:,1)) || any(x0 >= bnds(:,2)) )
            infeasible = ((x0 <= bnds(:,1)) | (x0 >= bnds(:,2)));
            infeasible_index = find(infeasible,1,'first');
            error(['The BFGS starting point for variable ',num2str(infeasible_index),...
                ' is outside the provided bounds.']);
        end
    end
    
    test_eval = f(x_map(T_inv*x0));
    try num_lines = length(test_eval);
    catch
        error('Function does not evaluate with the starting point and constraints provided. Check that the start point is the proper length.');
    end
    
    if num_lines ~= 1
        warning('The function to minimize does not return a scalar. The first index will be minimized, and the others will be treated as penalties.')
    end
    
    %Generate a logarithmic penalty function for f, g, and H based on the
    %boundaries of the space. If no boundaries are proivided, this function
    %is 0 everywhere.
    if isempty(bnds) || ~isequal(size(bnds),[n,2])
        f_bnd_penalty = @(x) 0;
    else
        %mu is very important, because it keeps the effects of the log
        %barrier function minimal when the argument of the log is greater
        %than 0. The same problem should be solved repeatedly with
        %decreasing values of mu and the previous solution as the starting
        %point. However, the user may simply specify one value to use.
        mu = bounds{2};
        %Scale the constraints appropriately
        f_bnd_penalty = @(x) min(0,mu*( sum(log(max(0,x - repmat(bounds{1}(:,1),1,size(x,2))))) + sum(log(max(0,-x + repmat(bounds{1}(:,2),1,size(x,2))))) ));
    end
    
    %This is our function to minimize, including our bound penalty.
%     f_bnded = @(x)f(x_map(T_inv*x)) - [f_bnd_penalty(x_map(T_inv*x));zeros(num_lines-1,1)];
    f_bnded = @(x)objective_function(f,x_map,T_inv,x,penalty_weight) - f_bnd_penalty(x_map(T_inv*x));
    
    x = NaN(n,max_it+1);
    x(:,1) = x0;
    
    f_eval = NaN(1,max_it+1);
    f_start = f_bnded(x0);
    f_eval(1) = f_start(1);
    
    gradient = NaN(n,max_it);
    
    %Data to record for plotting
    f_calls = zeros(3,max_it);
    descent_dir = NaN(n,max_it);
    step_length = NaN(1,max_it);
    
    gradient_convergence = NaN(1,max_it);
    x_convergence = NaN(1,max_it);
    
    direction_time = NaN(1,max_it);
    search_time = NaN(1,max_it);
    check_convergence_time = NaN(1,max_it);
    
    %Initialize previous values
    x_prev = zeros(n,1);
    g_prev = zeros(n,1);
    H_prev = eye(n);
    H_inv_prev = eye(n);
    
    %Flags for termination
    loop = 1;
    grad_tol_met = false;
    x_tol_met = false;
    nan_step = false;
    
    if print_options
        fprintf('***** BEGINNING BFGS MINIMIZER ***** \n');
    end
    
    while(loop < max_it && ~grad_tol_met && ~x_tol_met && ~nan_step)
        if print_options == 2
            fprintf('Beginning loop %i.\n',loop)
        end
        
        if mod(loop,5) == 0
            if print_options == 2
                fprintf('Updating dx...')
            end
            x_typical = abs(x(:,loop));
            finite_diff = max(calculate_dx(x_typical,bnds,coarse_dx),eps^(.5));
            if print_options == 2
                fprintf('Complete.\n')
            end
        end
        
        direction_timer = tic;
        
        %Find the new search direction
        x_curr = x(:,loop);
        f_curr = f_eval(loop);
        [p,s,g,H,H_inv,f_calls_direction] = bfgs_direction(f_bnded,x_curr,x_prev,f_curr,g_prev,H_prev,H_inv_prev,finite_diff,f_batch);
        
        %record the time used to find the new direction
        direction_time(loop) = toc(direction_timer);
        
        %save the output
        descent_dir(:,loop) = p;
        f_calls(1,loop) = f_calls_direction;
        x_prev = x_curr;
        g_prev = g;
        gradient(:,loop) = g;
        H_prev = H;
        H_inv_prev = H_inv;
        
        search_timer = tic;
        
        %Search along the search direction
        [x_search,s_search,f_search,f_calls_search] = bfgs_search(f_bnded,x_curr,f_curr,p,s,g,[to_backtrack,alpha,backtrack_ratio,max_linesearch_steps]);
        
        %record the time used to search for a new point
        search_time(loop) = toc(search_timer);
        
        %save the output
        f_calls(2,loop) = f_calls_search;
        step_length(loop) = s_search;
        
        %If any or all of the variables are NaN, then we do not want to
        %accept the new point, so set our flag.
        if any(isnan(x_search))
            nan_step = true;
        else
            nan_step = false;
        end
        
        convergence_timer = tic;
        %If we took a valid step, test it for convergence. If we did not
        %(line search returned nan), then attempt trust region from the
        %previous point.
        if ~nan_step
            x(:,loop+1) = x_search;
            f_eval(loop+1) = f_search;
            
            %test for convergence
            f_typical = max(abs(f_eval(loop+1)),1e-4);
            gradient_convergence(loop) = abs(g)'*abs(x_typical)/abs(f_typical);
            x_convergence(loop) = norm((x_search - x_prev)./max(abs(x_search),abs(x_typical)),Inf);
            
            if gradient_convergence(loop) <= grad_tol
                %We may stop; we can check trust region to make sure we are
                %close enough to a minimum
                trust_region_options = [alpha,eta,ddl_radius,true];
                [x_tr,p_tr,f_tr,g_tr,H_tr,H_inv_tr,~,improved_point_found,delta,f_calls_tr] = bfgs_trust_region(f_bnded,x_search,x_prev,f_search,g,H,H_inv,finite_diff,f_batch,trust_region_options);
                f_calls(3,loop) = f_calls_tr;
                if improved_point_found
                    %Overwrite our previous point.
                    x(:,loop+1) = x_tr;
                    f_eval(loop+1) = f_tr;
                    descent_dir(:,loop) = p_tr/norm(p_tr);
                    f_calls(3,loop) = f_calls_direction;
                    x_prev = x_tr;
                    g_prev = g_tr;
                    gradient(:,loop) = g;
                    H_prev = H_tr;
                    H_inv_prev = H_inv_tr;
                    ddl_radius = delta;
                else
                   %gradient condition met, quit
                    grad_tol_met = true;
                end
            elseif x_convergence(loop) <= x_tol && ~grad_tol_met
                %We may stop; we can check trust region to make sure we are
                %close enough to a minimum
                trust_region_options = [alpha,eta,ddl_radius,true];
                [x_tr,p_tr,f_tr,g_tr,H_tr,H_inv_tr,~,improved_point_found,delta,f_calls_tr] = bfgs_trust_region(f_bnded,x_search,x_prev,f_search,g,H,H_inv,finite_diff,f_batch,trust_region_options);
                f_calls(3,loop) = f_calls_tr;
                if improved_point_found
                    %Overwrite our previous point.
                    x(:,loop+1) = x_tr;
                    f_eval(loop+1) = f_tr;
                    descent_dir(:,loop) = p_tr/norm(p_tr);
                    f_calls(3,loop) = f_calls_direction;
                    x_prev = x_tr;
                    g_prev = g_tr;
                    gradient(:,loop) = g;
                    H_prev = H_tr;
                    H_inv_prev = H_inv_tr;
                    ddl_radius = delta;
                else
                    %stall condition met, quit
                    x_tol_met = true;
                end
            end
            
           	%Record the time required to check for convergence, and
           	%increment the variable if we have not converged.
            check_convergence_time(loop) = toc(convergence_timer);
            if grad_tol_met || x_tol_met
                %do not increment the counter
            else
                loop = loop + 1;
            end
        else
            %We have a nan step. Try trust region. If it finds a better
            %point, overwrite the NaNs in the saved variables, increment
            %the counter and continue to the next loop. If we can't find an
            %improved point, we are done (via nan_step flag).
            trust_region_options = [alpha,eta,ddl_radius,false];
            [x_tr,p_tr,f_tr,g_tr,H_tr,H_inv_tr,~,improved_point_found,delta,f_calls_tr] =...
                bfgs_trust_region(f_bnded,x(:,loop),[],f_eval(loop),g,H,H_inv,finite_diff,f_batch,trust_region_options);
            f_calls(3,loop) = f_calls_tr;
            if improved_point_found
                %Save this as a new point
                x(:,loop+1) = x_tr;
                f_eval(loop+1) = f_tr;
                descent_dir(:,loop) = p_tr/norm(p_tr);
                step_length(loop) = norm(p_tr);
                f_calls(3,loop) = f_calls_direction;
                x_prev = x_tr;
                g_prev = g_tr;
                gradient(:,loop) = g;
                H_prev = H_tr;
                H_inv_prev = H_inv_tr;
                ddl_radius = delta;    
                
                check_convergence_time(loop) = toc(convergence_timer);
                loop = loop + 1;
                nan_step = false;
            else
                %trust region couldn't help us.
                x(:,loop+1) = x(:,loop);
                f_eval(loop+1) = f_eval(loop);
                gradient_convergence(loop) = 0;
                x_convergence(loop) = 0;
                check_convergence_time(loop) = toc(convergence_timer);
            end
        end
    end
    
    if loop == max_it
        %For the sake of indexing and plotting, if we reach the maximum
        %iterations, reduce the loop number to max_it-1 to properly access
        %data.
        loop = loop - 1;
    end
    
    %We only want to keep the non-NaN indices.
    x = x(:,1:loop+1);
    f_eval = f_eval(1:loop+1);
    
    step_length = step_length(1:loop);
    descent_dir = descent_dir(:,1:loop);
    gradient_convergence = gradient_convergence(1:loop);
    x_convergence = x_convergence(1:loop);    
    gradient = gradient(:,1:loop);
    direction_time = direction_time(:,1:loop);
    search_time = search_time(:,1:loop);
    check_convergence_time = check_convergence_time(:,1:loop);
    f_calls = f_calls(:,1:loop);
    x_scaled = x;
    
    %We now need to unscale x and map it back to the original variables.
    x_unscaled = NaN(size(x_map(x(:,1)),1),size(x,2));
    for j=1:size(x_unscaled,2)
        x_unscaled(:,j) = x_map(T_inv*x(:,j));
    end
    x = x_unscaled;
    
    %Assign the final output values
    x_final = x(:,loop+1);
    f_final = f_eval(:,loop+1);

    if print_options
        if grad_tol_met == 1
            fprintf('Algorithm halted because an optimal point was found. Solution reached in %i iteration(s).\n',loop-1);
        elseif x_tol_met == 1
            fprintf('Algorithm halted because x stopped changing. Solution reached in %i iteration(s).\n',loop-1);
        elseif nan_step
            fprintf('Algorithm halted because the point cannot be improved. \nThe validity of this solution may be questionable. \nSolution reached in %i iteration(s).\n',loop-1);
        else
            disp('Routine Terminated. No solution found.')
        end
        fprintf('The objective function was evaluated %i times. \n \n',sum(sum(f_calls))+loop);
    end
    %Print a summary of performance
    if print_options
        fprintf('***** BFGS MINIMIZER PERFORMANCE SUMMARY *****\n')
        %How many lines need we print?
        num_lines = max(15,loop);
        num_lines = min(num_lines,max_it);
        %Which lines to print? first ten and last 5.
        if num_lines <= 15
            its = 1:loop;
        else
            its = [1:10,num_lines-4:num_lines];
        end
%         keyboard
        %Print each line
        for k=1:length(its)
            %Print the new direction
            fprintf('p_%i: [',its(k))
            for j=1:n
                if j ~= n
                    fprintf('%2.4f, ',descent_dir(j,its(k)));
                else
                    fprintf('%2.4f]  ',descent_dir(j,its(k)));
                end
            end
            %Print the step length
            fprintf('step_%i: ',its(k))
%             fprintf('%2.4f  ',step_len(its(k)));
            fprintf('%e  ',step_length(its(k)));
            %Print the latest solution
            fprintf('x_%i: [',its(k))
            for j=1:size(x,1)
                if j ~= size(x,1)
                    fprintf('%2.4f, ',x(j,its(k)));
                else
                    fprintf('%2.4f]  ',x(j,its(k)));
                end                    
            end
            
            fprintf('\n')
            if k == 10
                fprintf('...\n')
            end
        end
    end    
    
    %Plot metrics
    to_plot = zeros(9,1);
    to_plot(plot_options) = 1;
    if to_plot(1)
        %Plot function values by iteration
        figure
        clf
        if min(f_eval) > 0
            %We cannot plot negative points on a semilog plot
            semilogy(1:loop+1,f_eval,'.-','Linewidth',2);
        else
            plot(1:loop+1,f_eval,'.-','Linewidth',2);
        end
        hold on
        grid on
        title('Function Values by Iteration')
        xlabel('Iteration')
        ylabel('Value')
        hold off
    end
    if to_plot(2)
        %Plot convergence, defined as the norm of the ratio of each guess
        %and the final answer
        figure
        clf
        hold on
        title('Convergence of x (ratio of norms)')
        xlabel('Iteration')
        ylabel('|x_k|/|x_{k-1}|')
        x_norm = NaN(1,size(x_scaled,2));
        for i=1:length(x_norm)
            x_norm(i) = norm(x_scaled(:,i) - x_scaled(:,end));
        end
        grid on
        x_norm = x_norm(2:end)./x_norm(1:end-1);
        if min(x_norm) > 0
            semilogy(1:loop,x_norm,'.-','Linewidth',2);
        else
            plot(1:loop,x_norm,'.-','Linewidth',2);
        end
        hold off
    end
    if to_plot(3)
        %If the problem is one or two dimensional, plot the contour and path taken
        %by the algorithm
        if n == 2
            %We can plot xs and fs
            figure
            clf
            hold on
            title('x and f by Iteration')
            xlabel(['x_',num2str(nb_var(1))])
            ylabel(['x_',num2str(nb_var(2))])
            zlabel('f')
            plot3(x(nb_var(1),:),x(nb_var(2),:),f_eval,'g.-','Linewidth',2);
            plot3(x(nb_var(1),end),x(nb_var(2),end),f_eval(end),'r*');
            max_x1 = 1.5*max(abs(x(nb_var(1),:)));
            max_x2 = 1.5*max(abs(x(nb_var(2),:)));
            [X,Y] = meshgrid(linspace(-max_x1,max_x1,20),linspace(-max_x2,max_x2,20));
            
            max_x1_scaled = 1.5*max(abs(x_scaled(1,:)));
            max_x2_scaled = 1.5*max(abs(x_scaled(2,:)));
            [X_scaled,Y_scaled] = meshgrid(linspace(-max_x1_scaled,max_x1_scaled,20),linspace(-max_x2_scaled,max_x2_scaled,20));
            
            F = NaN(1,400);
            
            for j=1:400
                ftemp = f_bnded([X_scaled(j);Y_scaled(j)]);
                if ~isinf(ftemp)
                    F(j) = ftemp;
                end
            end
            F_sq = reshape(F,sqrt(length(F)),sqrt(length(F)));
            contour(X,Y,F_sq,100)
            hold off
        elseif n == 1
            %We can plot xs and fs
            figure
            clf
            hold on
            title('x and f by Iteration')
            xlabel(['x_',num2str(nb_var(1))])
            ylabel('f')
            
            plot(x(nb_var(1),:),f_eval,'g.-','Linewidth',2);
            plot(x(nb_var(1),end),f_eval(end),'r*');
            max_x1 = 1.5*max(abs(x(nb_var(1),:)));
            X = linspace(-max_x1,max_x1);
            
            %We know this is one-dimensional, because n = 1.
            max_x1_scaled = 1.5*max(abs(x_scaled));
            X_scaled = linspace(-max_x1_scaled,max_x1_scaled);
            
            F = NaN(length(X),1);
            for j=1:length(X)
                ftemp = f_bnded(X_scaled(j));
                F(j) = ftemp;
            end
            plot(X,F,'Linewidth',4);
            hold off
        end   
    end
    %Norm of the gradient (optimality)
    if to_plot(4)
        norm_grad = NaN(size(gradient,2),1);
        for j=1:length(norm_grad)
            norm_grad(j) = norm(gradient(:,j));
        end
        figure
        clf
        if min(norm_grad) > 0
            semilogy((1:loop),norm_grad,'.-','Linewidth',2)
        else
            plot((1:loop),norm_grad,'.-','Linewidth',2)
        end
        hold on
        title('Norm of the Gradient')
        xlabel('iteration')
        grid on
        hold off
    end
    %stalling convergence check
    if to_plot(5)
        figure
        clf
        hold on
        grid on
        if min(x_convergence) > 0
            semilogy(1:loop,x_convergence,'.-','Linewidth',2)
        else
            plot(1:loop,x_convergence,'.-','Linewidth',2)
        end
        title('Relative x for Termination')
        xlabel('iteration')
        hold off
    end
    %Gradient convergence check
    if to_plot(6)
        figure
        clf
        hold on
        grid on
        if min(gradient_convergence) > 0
            semilogy(1:loop,gradient_convergence,'.-','Linewidth',2)
        else
            plot(1:loop,gradient_convergence,'.-','Linewidth',2)
        end
        title('Relative Gradient for Termination')
        xlabel('iteration')
        hold off
    end
    %Step length
    if to_plot(7)
        figure
        clf
        if min(step_length) > 0
            semilogy(1:loop,step_length,'.-','Linewidth',2)
        else
            plot(1:loop,step_length,'.-','Linewidth',2)
        end
        hold on
        grid on
        title('Step Length by Iteration')
        xlabel('iteration')
        ylabel('step length magnitude')
        hold off
    end
    %function calls by iteration and process
    if to_plot(8)
        
        figure
        clf
        hold on
        plot(1:loop,f_calls(1,:),'.-','Linewidth',2)
        plot(1:loop,f_calls(2,:),'g.-','Linewidth',2)
        plot(1:loop,f_calls(3,:),'magenta.-','Linewidth',2)
        plot(1:loop,sum(f_calls),'r.-','Linewidth',2)
        ylim([0 max(sum(f_calls))+1])
        total_str = strcat('Total (',num2str(sum(sum(f_calls)+1)),')');
%         legend('Direction finding','Line search','Trust region',total_str,'Location','EastOutside')
        legend('Direction finding','Line search','Trust region',total_str,'Location','EastOutside')
        title('Function Calls per Iteration');
        xlabel('iteration')
        ylabel('Function Calls');
        grid on
        hold off
    end
    %time consumed by iteration and process
    if to_plot(9)
        figure
        clf
        hold on
        plot(1:loop,direction_time,'.-','Linewidth',2)
        plot(1:loop,search_time+direction_time,'g.-','Linewidth',2)
        plot(1:loop,check_convergence_time+direction_time+search_time,'r.-','Linewidth',2)
        legend('Direction','Plus Step','Plus Conv Check (total)')
        title_str = strcat('Stacked Line Graph of Time per Iteration (s) (Total: ',num2str(sum([sum(direction_time),sum(search_time),sum(check_convergence_time)])),'s)');
        title(title_str)
        xlabel('iteration')
        ylabel('time (s)')
        grid on
        hold off        
    end
    
    if print_options
        fprintf('\n***** BFGS MINIMIZER COMPLETE *****\n\n')
        
    end
    stats = {f_eval,x_scaled,x,gradient,x_convergence,gradient_convergence,step_length,f_calls,[direction_time;search_time;check_convergence_time],loop};
    
end
