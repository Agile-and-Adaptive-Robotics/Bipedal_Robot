function [ x_final, f_final ] = bnd_con_pen_bfgs(f, x0, bounds, linear_constraints,...
    conv_criteria, direction_options, search_options, plot_options,...
    print_options)
%BND_BFGS Find roots or minima of a function with a given starting point.
%Descent direction selected by BFGS, step length can be adjusted through a 
%line search (backtrack or cubic fit). The final x vector
%and f value are returned.
%
%   INPUTS:
%   f - function to minimize.
%
%   x0 - initial point (n-element column vector).
%
%   bounds - limits on each variable (n x 2 matrix).
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
    
    try f_orig = f(x0);
    catch
        error('The function to minimize does not evaluate at the initial point.')
    end
    f_start_oom = ceil(log10(abs(f_orig(1))));
    mu = max(1e-6,10^(f_start_oom-3));
    if length(f_orig) > 1
        f_pen_mag = sum(f_orig(2:end));
        if f_pen_mag == 0
            penalty_weight = 0;
        else
            penalty_weight = 1/f_pen_mag;
        end
    else
        penalty_weight = 0;
    end
    
    n = length(x0);
    
    A_bnds = [];
    b_bnds = [];
    if isempty(bounds) && length(f_orig) == 1
        warning('If no bounds or penalties are needed, use bfgs.m, which is unbounded.')
        bnds = [];
    else
        if size(bounds,1) == length(x0) && size(bounds,2) == 2
            bnds = {bounds,mu};
            %If any of the bounds leave no room in the middle (upper bound
            %equals lower bound), then create a constraint to eliminiate
            %that variable.
            vars_to_remove = find(bounds(:,1) == bounds(:,2));
            num_vars_to_remove = length(vars_to_remove);
            A_bnds = zeros(num_vars_to_remove,n);
            b_bnds = zeros(num_vars_to_remove,1);
            for i=1:num_vars_to_remove
                A_bnds(i,vars_to_remove(i)) = 1;
            end
        elseif isempty(bounds)
            bnds = [];
        else
            error('Bounds are not the correct size.');
        end
    end
    
    %Parse the linear constraints so we can add additional constraints if
    %the bounds are too restrictive.
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
    linear_constraints{1} = [A;A_bnds];
    linear_constraints{2} = [b;b_bnds];
    %****************
    
    if length(conv_criteria) >= 4
        x_tol = conv_criteria(4);
        conv_criteria = conv_criteria(1:3);
    else
        x_tol = 1e-6;
    end
    
    max_it = 16;
    
    xs = NaN(n,max_it);
    xs(:,1) = x0;
        
    to_plot = zeros(11,1);
    to_plot(plot_options) = 1;

    f_eval = cell(max_it,1);
    x_scaled = cell(max_it,1);
    x = cell(max_it,1);
    gradient = cell(max_it,1);
    x_convergence = cell(max_it,1);
    gradient_convergence = cell(max_it,1);
    step_length = cell(max_it,1);
    f_calls = cell(max_it,1);
    time = cell(max_it,1);
    num_loops = NaN(max_it,1);
    mus = NaN(max_it,1);
    
    if print_options
        disp('***** BEGINNING INTERIOR POINT BFGS MINIMIZER *****')
    end
    loop = 1;
    converged = false;
    while loop <= max_it && ~converged && mu > eps
        [x_interim,f_interim,stats] = bfgs(f, xs(:,loop), bnds, linear_constraints,...
        penalty_weight,conv_criteria, direction_options, search_options, [], print_options);
        
        if ~iscell(stats)
            loop = max_it + 1;
        else
            %Assign the new x value
            xs(:,loop+1) = x_interim;
            
            %Compare its change from the previous. If it has not changed
            %much or mu is too small, we have converged.
            converged = (norm(xs(:,loop+1)-xs(:,loop)) < x_tol || mu < eps || penalty_weight > 1/eps);
            
            %Store the output from this trial. We will plot it all
            %together.            
            f_eval{loop} = stats{1};
            x_scaled{loop} = stats{2};
            x{loop} = stats{3};
            gradient{loop} = stats{4};
            x_convergence{loop} = stats{5};
            gradient_convergence{loop} = stats{6};
            step_length{loop} = stats{7};
            f_calls{loop} = stats{8};
            time{loop} = stats{9};
            num_loops(loop) = stats{10}; 
            mus(loop) = mu;

            %We need to update mu. It should always get smaller. The basic
            %technique is to reduce it by an order of magnitude every time.
            %However, we never want it to be larger than 3 orders of
            %magnitude less than the function's value. Therefore we should
            %first compute mu that way, and default to the basic rule if it
            %does not decrease or if f_eval is negative (mu is imaginary).
            mu_prev = mu;
            f_last = f_eval{loop}(end);
            if f_last > 0
                mu_new = 10^(ceil(log10(f_last)) - 3);
            else
                mu_new = 10^(floor(log10(abs(f_last))) - 3);
            end
            if mu_new >= mu_prev || ~isreal(mu_new)
                mu = mu/10;
            else
                mu = mu_new;
            end
            
            %Increase the penalty weight by an order of magnitude
            penalty_weight = penalty_weight*10;
%             penalty_weight_prev = penalty_weight;
%             penalty_weight_new = 

            %update bnds
            bnds = {bounds,mu};
            
            if loop < max_it && ~converged && mu > eps
                %increment the counter
                loop = loop + 1;
            end
        end
    end
    
    if print_options
        if converged
            disp('Algorithm halted because x stopped changing.')
        else
            disp('Algorithm halted because mu is less than machine epsilon.')
        end
    end
    
    %Only keep the counts that were populated
    num_loops = num_loops(1:loop);
    
    %Our most recent guesses for x and f are now our final answers
    x_final = x_interim;
    f_final = f_interim;
    
    %We want to know how each iteration is sequenced, spanning calls to
    %bfgs. This will let us plot all of the data in one graph, rather than
    %in many separate ones.
    %For "loop" based metrics (step length, gradient, things that occur
    %BETWEEN guesses for x), we come can find this relationship by
    %cumulatively summing the numer of loops.
    end_pts_loops = cumsum(num_loops);
    start_pts_loops = end_pts_loops-num_loops+1;
    
    %For "point" based metrics (x, f, things that occur for EACH guess of
    %x), we will have one additional point per call to bfgs, so account for
    %that.
    end_pts_points = cumsum(num_loops+1);
    start_pts_points = end_pts_points-(num_loops+1)+1;
    
    if to_plot(1)
        %Plot function values by iteration
        figure
        clf
        
        f_eval_min = 1;
        for i=1:loop
            f_eval_min = min(f_eval_min,min(f_eval{i}));
        end           
        
        for i=1:loop
            if f_eval_min > 0
                semilogy(start_pts_points(i):end_pts_points(i),f_eval{i},'blue','Linewidth',2)
                if i == 1
                    hold on
                end
                semilogy(start_pts_points(i),f_eval{i}(1),'go');
                semilogy(end_pts_points(i),f_eval{i}(end),'r*');
            else
                hold on
                plot(start_pts_points(i):end_pts_points(i),f_eval{i},'blue','Linewidth',2)
                plot(start_pts_points(i),f_eval{i}(1),'go');
                plot(end_pts_points(i),f_eval{i}(end),'r*');
            end
        end
        grid on
        title('Function Values by Iteration')
        xlabel('Iteration')
        ylabel('Value')
        hold off
    end
    
    if to_plot(2)
        %Plot convergence
        figure
        clf
        x_norm = cell(loop,1);
        x_norm_min = 1;
        for i=1:loop
            x_norm_min = min(x_norm_min,min(x_norm{i}));
        end
        
        for i=1:loop
            for j=1:num_loops(i)+1
                x_norm{i}(j) = norm(x_scaled{i}(:,j) - x_scaled{i}(:,end));
            end
            x_norm{i} = x_norm{i}(2:end)./x_norm{i}(1:end-1);
            x_norm{i}(isnan(x_norm{i})) = 0;
            if x_norm_min > 0
                semilogy(start_pts_loops(i):end_pts_loops(i),x_norm{i},'blue','Linewidth',2)
                semilogy(start_pts_loops(i),x_norm{i}(1),'go');
                semilogy(end_pts_loops(i),x_norm{i}(end),'r*');
            else
                plot(start_pts_loops(i):end_pts_loops(i),x_norm{i},'blue','Linewidth',2)
                plot(start_pts_loops(i),x_norm{i}(1),'go');
                plot(end_pts_loops(i),x_norm{i}(end),'r*');
            end
            if i == 1
                hold on
            end
        end
        title('Convergence of x (ratio of norms)')
        xlabel('Iteration')
        ylabel('|x_k|/|x_{k-1}|')
        grid on
        hold off
    end
    %Plot the path of the optimizer over the contour.
    if to_plot(3)
        if size(x0,1) == 2 && ~isempty(bounds)
            %We can plot xs and fs
            figure
            clf
            hold on
            title('x and f by Iteration')
            xlabel('x_1')
            ylabel('x_2')
            zlabel('f')
            for i=1:loop
                plot3(x{i}(1,:),x{i}(2,:),f_eval{i},'g','Linewidth',2);
                plot3(x{i}(1,end),x{i}(2,end),f_eval{i}(end),'r*');
            end
            [X,Y] = meshgrid(linspace(bounds(1,1),bounds(1,2),20),linspace(bounds(2,1),bounds(2,2),20));
            xtemp = X(:);
            ytemp = Y(:);
            F = NaN(400,1);
            for j=1:400
                F(j) = f([xtemp(j);ytemp(j)]);
            end
            F_sq = reshape(F,sqrt(length(F)),sqrt(length(F)));
            contour(X,Y,F_sq,100)
            grid on
            hold off
        elseif size(x0,1) == 2 && isempty(bounds)
            %We can plot xs and fs
            figure
            clf
            hold on
            title('x and f by Iteration')
            xlabel('x_1')
            ylabel('x_2')
            zlabel('f')
            max_x1 = 0;
            max_x2 = 0;
            for i=1:loop
                plot3(x{i}(1,:),x{i}(2,:),f_eval{i},'g','Linewidth',2);
                plot3(x{i}(1,end),x{i}(2,end),f_eval{i}(end),'r*');
                max_x1 = max(max_x1,max(abs(x{i}(1,:))));
                max_x2 = max(max_x2,max(abs(x{i}(2,:))));
            end
            max_x1 = 1.5*max_x1;
            max_x2 = 1.5*max_x2;
            [X,Y] = meshgrid(linspace(-max_x1,max_x1,20),linspace(-max_x2,max_x2,20));
            xtemp = X(:);
            ytemp = Y(:);
            F = NaN(400,1);
            for j=1:400
                f_obj = f([xtemp(j);ytemp(j)]);
                F(j) = f_obj(1);
            end
            F_sq = reshape(F,sqrt(length(F)),sqrt(length(F)));
            contour(X,Y,F_sq,100)
            grid on
            hold off
        end        
    end
    %Norm of the gradient (optimality)
    if to_plot(4)
        figure
        clf
        norm_grad = cell(size(gradient));
        
        norm_grad_min = 1;
        for i=1:loop
            norm_grad{i} = NaN(num_loops(i),1);
            for j=1:num_loops(i)
                norm_grad{i}(j) = norm(gradient{i}(:,j));
            end
            norm_grad_min = min(norm_grad_min,min(norm_grad{i}));
        end           
        
        for i=1:loop
            if norm_grad_min > 0
                semilogy(start_pts_loops(i):end_pts_loops(i),norm_grad{i},'blue','Linewidth',2)
                if i == 1
                    hold on
                end
                semilogy(start_pts_loops(i),norm_grad{i}(1),'go');
                semilogy(end_pts_loops(i),norm_grad{i}(end),'r*');
            else
                hold on
                plot(start_pts_loops(i):end_pts_loops(i),norm_grad{i},'blue','Linewidth',2)
                plot(start_pts_loops(i),norm_grad{i}(1),'go');
                plot(end_pts_loops(i),norm_grad{i}(end),'r*');
            end
        end
        grid on
        title('Norm of the Gradient')
        xlabel('Iteration')
        ylabel('Value')
        hold off
    end
    
    if to_plot(5)
        %Plot stalling termination metric by iteration
        figure
        clf
        x_conv_min = 1;
        for i=1:loop
            x_conv_min = min(x_conv_min,min(x_convergence{i}));
        end           
        
        for i=1:loop
            if x_conv_min > 0
                semilogy(start_pts_loops(i):end_pts_loops(i),x_convergence{i},'blue','Linewidth',2)
                if i == 1
                    hold on
                end
                semilogy(start_pts_loops(i),x_convergence{i}(1),'go');
                semilogy(end_pts_loops(i),x_convergence{i}(end),'r*');
            else
                hold on
                plot(start_pts_loops(i):end_pts_loops(i),x_convergence{i},'blue','Linewidth',2)
                plot(start_pts_loops(i),x_convergence{i}(1),'go');
                plot(end_pts_loops(i),x_convergence{i}(end),'r*');
            end
        end
        grid on
        title('Relative x for Termination')
        xlabel('Iteration')
        hold off
        
    end
    if to_plot(6)
        %Plot gradient termination metric by iteration
        figure
        clf
        grad_conv_min = 1;
        for i=1:loop
            grad_conv_min = min(grad_conv_min,min(gradient_convergence{i}));
        end           
        
        for i=1:loop
            if grad_conv_min > 0
                semilogy(start_pts_loops(i):end_pts_loops(i),gradient_convergence{i},'blue','Linewidth',2)
                if i == 1
                    hold on
                end
                semilogy(start_pts_loops(i),gradient_convergence{i}(1),'go');
                semilogy(end_pts_loops(i),gradient_convergence{i}(end),'r*');
            else
                hold on
                plot(start_pts_loops(i):end_pts_loops(i),gradient_convergence{i},'blue','Linewidth',2)
                plot(start_pts_loops(i),gradient_convergence{i}(1),'go');
                plot(end_pts_loops(i),gradient_convergence{i}(end),'r*');
            end
        end
        grid on
        title('Relative gradient for Termination')
        xlabel('Iteration')
        hold off
        
    end

    if to_plot(7)
        %Plot step length by iteration
        figure
        clf
        step_len_min = 1;
        for i=1:loop
            step_len_min = min(step_len_min,min(step_length{i}));
        end           
        for i=1:loop
            if step_len_min > 0
                semilogy(start_pts_loops(i):end_pts_loops(i),step_length{i},'blue','Linewidth',2)
                if i == 1
                    hold on
                end
                semilogy(start_pts_loops(i),step_length{i}(1),'go');
                semilogy(end_pts_loops(i),step_length{i}(end),'r*');
            else
                hold on
                plot(start_pts_loops(i):end_pts_loops(i),step_length{i},'blue','Linewidth',2)
                plot(start_pts_loops(i),step_length{i}(1),'go');
                plot(end_pts_loops(i),step_length{i}(end),'r*');
            end
        end
        grid on
        title('Scaled Step Length')
        xlabel('Iteration')
        hold off
    end
    %Calculate the total number of function evaluations
    f_call_sum = 0;
    for i=1:loop
        f_call_sum = f_call_sum + sum(sum(f_calls{i}));
    end
    %Plot the number of function calls per iteration
    if to_plot(8)
        figure
        clf
        hold on
        
        for i=1:loop
            plot(start_pts_loops(i):end_pts_loops(i),f_calls{i}(1,:),'blue','Linewidth',2)
            plot(start_pts_loops(i):end_pts_loops(i),f_calls{i}(2,:),'green','Linewidth',2)
            plot(start_pts_loops(i):end_pts_loops(i),sum(f_calls{i}),'r','Linewidth',2)
            plot(end_pts_loops(i),sum(f_calls{i}(:,end)),'r*','Linewidth',2)
        end
        
        total_str = strcat('Total (',num2str(f_call_sum),')');
        legend('Direction finding','Line search',total_str,'Location','EastOutside')
        title('Function Calls per Iteration');
        xlabel('Iteration')
        ylabel('Function Calls');
        hold off
    end
    %Plot the time spend per iteration
    if to_plot(9) 
        
        figure
        clf
        hold on
        time_sum = 0;
        for i=1:loop
            plot(start_pts_loops(i):end_pts_loops(i),time{i}(1,:),'blue','Linewidth',2)
            plot(start_pts_loops(i):end_pts_loops(i),time{i}(2,:),'green','Linewidth',2)
            plot(start_pts_loops(i):end_pts_loops(i),time{i}(3,:),'red','Linewidth',2)
            plot(start_pts_loops(i):end_pts_loops(i),sum(time{i}),'magenta','Linewidth',2)
            plot(end_pts_loops(i),sum(time{i}(:,end)),'r*','Linewidth',2)
            time_sum = time_sum + sum(sum(time{i}(~isnan(time{i}))));
        end
        total_str = strcat('Total (',num2str(time_sum),'s)');
        legend('Direction finding','Line search','Check Convergence',total_str,'Location','EastOutside')
        title('Time Spent per Iteration');
        xlabel('Iteration')
        ylabel('Time (s)');
        hold off
    end
    %Plot the number of iterations per loop (spanning bfgs() calls)
    if to_plot(10)
        figure
        clf
        hold on
        grid on
        plot(num_loops,'Linewidth',2)
        title('Number of loops per bfgs() call')
        xlabel('Iteration')
        hold off
    end
    if to_plot(11)
        figure
        clf
        semilogy(mus,'Linewidth',2)
        hold on
        grid on
        title('\mu')
        xlabel('Iteration')
        hold off
    end
    
    if print_options
        fprintf('The solution was reached in %i bfgs() calls,\namounting to %i quasi-Newton steps and %i function evaluations.\n',loop,sum(num_loops),f_call_sum)
        fprintf('***** INTERIOR POINT BFGS MINIMIZER COMPLETE *****\n\n')
    end
end

