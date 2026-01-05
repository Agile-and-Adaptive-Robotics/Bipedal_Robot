function [x_final, f_final, x_all] = bnd_nelder_mead( f, x0, bounds, linear_constraints, conv_criteria, params, plot_options, print_options )
%NELDER_MEAD_BND Summary of this function goes here
%   Detailed explanation goes here

    try f_orig = f(x0);
    catch
        error('The function to minimize does not evaluate at the initial point.')
    end
    oom_diff = 3;
    mu = max(1e-12,10^(ceil(log10(abs(f_orig))) - oom_diff));
    n = length(x0);
    
    A_bnds = [];
    b_bnds = [];
    if isempty(bounds)
        warning('If no bounds are needed, use nelder_mead.m, which is unbounded.')
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
        x_tol = 1e-12;
    end
    
    max_it = 16;
    
    xs = NaN(n,max_it);
    f_eval = NaN(max_it,1);
    xs(:,1) = x0;
    x_interim = x0;
    x_all = x0;
    f_avg = cell(max_it,1);
    f_best = cell(max_it,1);
    f_worst = cell(max_it,1);
    max_x_dist = cell(max_it,1);
    max_f_dist = cell(max_it,1);
    f_calls = cell(max_it,1);
    time = cell(max_it,1);
    centroids = cell(max_it,1);
    num_loops = NaN(max_it,1);
    mus = NaN(max_it,1);
    
    to_plot = zeros(11,1);
    to_plot(plot_options) = 1;
    penalty_weight = 0; %If we want penalties, use bnd_con_pen_nelder_mead
    
    if print_options
        disp('***** BEGINNING INTERIOR POINT NELDER-MEAD SIMPLEX MINIMIZER *****')
    end
    loop = 1;
    converged = false;
    while loop <= max_it && ~converged && mu > eps
        [x_interim,f_interim,x_all,stats] = nelder_mead(f, x_all, bnds, linear_constraints,...
        penalty_weight, conv_criteria, params, [], 2);
    
%         [x_interim,f_interim,x_all,stats] = nelder_mead(f, x_interim, bnds, linear_constraints,...
%         conv_criteria, params, 1, 2);
    
        xs(:,loop+1) = x_interim;
        f_eval(loop) = f_interim;
    
        %Compare its change from the previous. If it has not changed
        %much or mu is too small, we have converged.
        converged = (norm(xs(:,loop+1)-xs(:,loop)) < x_tol || mu < eps);

        %Store the output from this trial. We will plot it all
        %together.     
        f_avg{loop} = stats{1};
        f_best{loop} = stats{2};
        f_worst{loop} = stats{3};
        max_x_dist{loop} = stats{4};
        max_f_dist{loop} = stats{5};
        f_calls{loop} = stats{6};
        time{loop} = stats{7};
        num_loops(loop) = stats{8}; 
        mus(loop) = mu;
        centroids{loop} = stats{9};

        %We need to update mu. It should always get smaller. The basic
        %technique is to reduce it by an order of magnitude every time.
        %However, we never want it to be larger than 3 orders of
        %magnitude less than the function's value. Therefore we should
        %first compute mu that way, and default to the basic rule if it
        %does not decrease or if f_eval is negative (mu is imaginary).
        mu_prev = mu;
        f_last = f_eval(loop);
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

        %update bnds
        bnds = {bounds,mu};

        if loop < max_it && ~converged && mu > eps
            %increment the counter
            loop = loop + 1;
        end
    end
    
    if print_options
        if converged
            disp('Algorithm halted because x stopped changing.')
        else
            disp('Algorithm halted because mu is less than machine epsilon.')
        end
    end
    
    num_loops = num_loops(1:loop);
    
    %We want to know how each iteration is sequenced, spanning calls to
    %bfgs. This will let us plot all of the data in one graph, rather than
    %in many separate ones.
    %For "loop" based metrics (step length, gradient, things that occur
    %BETWEEN guesses for x), we come can find this relationship by
    %cumulatively summing the numer of loops.
    end_pts_loops = cumsum(num_loops);
    start_pts_loops = end_pts_loops-num_loops+1;
    
    if to_plot(1)
        %Plot function values by iteration
        figure
        clf
        
        f_eval_min = 1;
        for i=1:loop
            f_eval_min = min(f_eval_min,min(f_avg{i}));
            f_eval_min = min(f_eval_min,min(f_best{i}));
            f_eval_min = min(f_eval_min,min(f_worst{i}));
        end           
        
        for i=1:loop
            if f_eval_min > 0
                semilogy(start_pts_loops(i):end_pts_loops(i),f_avg{i},'blue','Linewidth',2)
                if i == 1
                    hold on
                end
                semilogy(start_pts_loops(i),f_avg{i}(1),'blueo');
                semilogy(end_pts_loops(i),f_avg{i}(end),'blue*');
                
                semilogy(start_pts_loops(i):end_pts_loops(i),f_best{i},'green','Linewidth',2)
                semilogy(start_pts_loops(i),f_best{i}(1),'greeno');
                semilogy(end_pts_loops(i),f_best{i}(end),'green*');
                
                semilogy(start_pts_loops(i):end_pts_loops(i),f_worst{i},'red','Linewidth',2)
                semilogy(start_pts_loops(i),f_worst{i}(1),'redo');
                semilogy(end_pts_loops(i),f_worst{i}(end),'red*');
            else
                hold on
                plot(start_pts_loops(i):end_pts_loops(i),f_avg{i},'blue','Linewidth',2)
                plot(start_pts_loops(i),f_avg{i}(1),'blueo');
                plot(end_pts_loops(i),f_avg{i}(end),'blue*');
                
                plot(start_pts_loops(i):end_pts_loops(i),f_best{i},'green','Linewidth',2)
                plot(start_pts_loops(i),f_best{i}(1),'greeno');
                plot(end_pts_loops(i),f_best{i}(end),'green*');
                
                plot(start_pts_loops(i):end_pts_loops(i),f_worst{i},'red','Linewidth',2)
                plot(start_pts_loops(i),f_worst{i}(1),'redo');
                plot(end_pts_loops(i),f_worst{i}(end),'red*');
            end
        end
        grid on
        title('Function Values by Iteration')
        xlabel('Iteration')
        ylabel('Value')
        hold off
    end
    
    if to_plot(2) 
        figure
        clf
        for i=1:loop
            semilogy(start_pts_loops(i):end_pts_loops(i),max_x_dist{i},'blue','Linewidth',2)
            if i == 1
                hold on
            end
            semilogy(start_pts_loops(i),max_x_dist{i}(1),'greeno');
            semilogy(end_pts_loops(i),max_x_dist{i}(end),'red*');           
        end
        grid on
        title('Maximum Simplex Chord')
        xlabel('Iteration')
        ylabel('Value')
        hold off
    end
    
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
                plot3(centroids{i}(1,:),centroids{i}(2,:),f_avg{i},'g','Linewidth',2);
                plot3(centroids{i}(1,end),centroids{i}(2,end),f_avg{i}(end),'r*');
            end
            [X,Y] = meshgrid(linspace(bounds(1,1),bounds(1,2),20),linspace(bounds(2,1),bounds(2,2),20));
            xtemp = X(:);
            ytemp = Y(:);
            for j=1:400
                F(j) = f([xtemp(j);ytemp(j)]);
            end
            F_sq = reshape(F,sqrt(length(F)),sqrt(length(F)));
            contour(X,Y,F_sq,100)
            grid on
            hold off
        end        
    end
    
    if to_plot(4) 
        figure
        clf
        for i=1:loop
            semilogy(start_pts_loops(i):end_pts_loops(i),max_f_dist{i},'blue','Linewidth',2)
            if i == 1
                hold on
            end
            semilogy(start_pts_loops(i),max_f_dist{i}(1),'greeno');
            semilogy(end_pts_loops(i),max_f_dist{i}(end),'red*');           
        end
        grid on
        title('Largest Function Value Ratio within Simplex')
        xlabel('Iteration')
        ylabel('Value')
        hold off
    end
    
    %Calculate the total number of function evaluations
    f_call_sum = 0;
    f_call_max = 1;
    for i=1:loop
        f_call_sum = f_call_sum + sum(f_calls{i});
        f_call_max = max(f_call_max,max(f_calls{i}));
    end
    %Plot the number of function calls per iteration
    if to_plot(5)
        figure
        clf
        hold on
        for i=1:loop
            stairs(start_pts_loops(i):end_pts_loops(i),f_calls{i},'blue','Linewidth',2)
            stairs(start_pts_loops(i),sum(f_calls{i}(1)),'go','Linewidth',2)
            stairs(end_pts_loops(i),sum(f_calls{i}(end)),'r*','Linewidth',2)
        end
        ylim([0,1.1*f_call_max])
        total_str = strcat('Function calls (',num2str(f_call_sum),')');
        legend(total_str,'Location','EastOutside')
        title('Function Calls per Iteration');
        xlabel('Iteration')
        ylabel('Function Calls');
        grid on
        hold off
    end
    
    %Calculate the total time taken
    time_sum = 0;
    time_max = 1e-12;
    for i=1:loop
        time_sum = time_sum + sum(sum(time{i}(:)));
        time_max = max(time_max,max(time{i}(:)));
    end
    
    %Plot the time per iteration
    if to_plot(6)        
        figure
        clf
        hold on
        for i=1:loop
            plot(start_pts_loops(i):end_pts_loops(i),time{i}(:,1),'blue','Linewidth',2)
            plot(start_pts_loops(i):end_pts_loops(i),time{i}(:,1)+time{i}(:,2),'magenta','Linewidth',2)
        end
        ylim([0,1.1*time_max])
        total_str = strcat('Total (',num2str(time_sum),')');
        legend('Evaluation',total_str,'Location','EastOutside')
        title('Time Taken per Iteration');
        xlabel('Iteration')
        ylabel('Function Calls');
        grid on
        hold off
    end
    
    if to_plot(7)
        figure
        clf
        semilogy(mus,'Linewidth',2)
        hold on
        grid on
        title('\mu')
        xlabel('Iteration')
        hold off
    end
    
    x_final = x_interim;
    f_final = f_interim;
    if print_options
        fprintf('The solution was reached in %i nelder_mead() calls,\namounting to %i steps and %i function evaluations.\n',loop,sum(num_loops),f_call_sum)
        fprintf('***** INTERIOR POINT NELDER-MEAD SIMPLEX MINIMIZER COMPLETE *****\n\n')
    end
    
end

