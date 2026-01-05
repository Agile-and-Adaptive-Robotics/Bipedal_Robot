function [ x_final, f_final, x_all, stats ] = nelder_mead( f, x0, bounds, ...
    linear_constraints, penalty_weight, conv_crit, params, to_plot, to_print )
%NELDER_MEAD Nelder-Mead Simplex function minimizer with boundaries
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
%   iterations, the x stalling tolerance for convergence, and the f stalling
%   tolerance for convergence. Default values are [1000,1e-6,1e-6].
%
%   plot_options - vector of what to plot (function values, x
%   ratios, if possible, x values on a contour of the function, norm(grad),
%   relative gradient (convergence criterion), relative x (convergence
%   criterion), step length, function calls by subroutine
%
%   print_options - 0 - Plot nothing
%                   1 - Plot results at the end of routine
%                   2 - Plot updates at each iteration 
    
    if isempty(params) || length(params) ~= 4
        alpha = 1;
        gamma = 2;
        rho = -.5;
        sigma = .5;
    else
        alpha = params(1);
        gamma = params(2);
        rho = params(3);
        sigma = params(4);
    end
    
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
    [x_map,~,nb_var] = elim_var(A,b,x0(:,1),allow_infeasible_start);
    x0 = x0(nb_var,:);
    
    %How large is our starting
    [n,num_individuals] = size(x0);
    
    %We may need to scale the variables. We can scale them to the same
    %order of magnitude, solve the problem, and then scale them back.
    %We need to do this now, because calculating our finite difference
    %requires that the problem be scaled.
    T = scale(x0);
    
    %Save the inverse to avoid computing it often. This is not expensive 
    %since T is diagonal.
    T_inv = T^-1;
    %scale the initial point
    for i=1:num_individuals
        x0(:,i) = T*x0(:,i);
    end
    
    if (isempty(bounds) && ismatrix(bounds)) || iscell(bounds)
        if (isempty(bounds) && ismatrix(bounds))
            bounds = {[],[]};
        end
        %input is probably ok
    else
        error('Bounds must be an empty array or a cell array with a boundary matrix and mu value (see nelder_mead.m)')
    end
    
    test_eval = f(x_map(T_inv*x0(:,1)));
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
    if isempty(bounds{1})
        bnds = [];
    else
        bnds = T*bounds{1}(nb_var,:);
    end
    
    if isempty(bounds) || ~isequal(size(bnds),[n,2])
        f_bnd_penalty = @(x) 0;
    else
        %mu is very important, because it keeps the effects of the log
        %barrier function minimal when the argument of the log is greater
        %than 0. The same problem should be solved repeatedly with
        %decreasing values of mu and the previous solution as the starting
        %point. However, the user may simply specify one value to use.
        mu = bounds{2};
        %Scale the constraints appropriately
%         f_bnd_penalty = @(x) mu*( sum(log(max(0,x - repmat(bounds{1}(:,1),1,size(x,2))))) + sum(log(max(0,-x + repmat(bounds{1}(:,2),1,size(x,2))))) );
        f_bnd_penalty = @(x) min(0,mu*( sum(log(max(0,x - repmat(bounds{1}(:,1),1,size(x,2))))) + sum(log(max(0,-x + repmat(bounds{1}(:,2),1,size(x,2))))) ));
    end
    
    %This is our function to minimize, including our bound penalty.
%     f_bnded = @(x)f(x_map(T_inv*x)) - f_bnd_penalty(x_map(T_inv*x));
    f_bnded = @(x)objective_function(f,x_map,T_inv,x,penalty_weight) - f_bnd_penalty(x_map(T_inv*x));
    
    if num_individuals == n+1
        %The user provided the entire starting population
        to_initialize = false;
    else
        %We need to generate the starting population
        to_initialize = true;
    end
    
    if isempty(conv_crit)
        max_it = 1000;
        x_stall = 1e-6;
        f_stall = 1e-6;
    else
        if length(conv_crit) == 1
            max_it = conv_crit(1);
            x_stall = 1e-6;
            f_stall = 1e-6;
        elseif length(conv_crit) == 2
            max_it = conv_crit(1);
            x_stall = conv_crit(2);
            f_stall = 1e-6;
        elseif length(conv_crit) == 3
            max_it = conv_crit(1);
            x_stall = conv_crit(2);
            f_stall = conv_crit(3);
        end
    end
    
    f_calls = zeros(max_it,1);
    f_avg = NaN(max_it,1);
    f_best = NaN(max_it,1);
    f_worst = NaN(max_it,1);
    max_x_dist = NaN(max_it,1);
    max_f_dist = NaN(max_it,1);
    eval_time = NaN(max_it,1);
    check_conv_time = NaN(max_it,1);
    centroids = NaN(n,max_it);
    x_best = NaN(n,max_it);
    
    %Compute the initial point, and calculate some basic system properties
    f_c = f_bnded(x0(:,1));
    cost = zeros(n+1,1);
    %We may need to create the first population.
    if to_initialize
        %Similar to the GA, make a matrix in which each column is an
        %individual
        population = zeros(n,n+1);
        
        population(:,1) = x0(:,1);
        cost(1) = f_c;
        for i=1:n
            
            if isempty(bnds)
                x_try = x0 + .01*T_inv*(rand(n,1)-.5);
            else
%                 x_try = x0 + .01*(rand(n,1)-.5).*range(bnds,2);
                x_try = x0 + .01*T_inv*(rand(n,1)-.5);
                x_try_unacceptable = any(x_try < bnds(:,1)) || any(x_try > bnds(:,2));

                if x_try_unacceptable
                    unacceptable = (x_try < bnds(:,1)) | (x_try > bnds(:,2));
                    x_try(unacceptable == true) = bnds(unacceptable == true,1) + ...
                        rand(length(x_try(unacceptable == true)),1).*range(bnds(unacceptable == true,:),2);
                    warning('out of bounds.')
                end

                if any(x_try < bnds(:,1)) || any(x_try > bnds(:,2))
                    warning('The Nelder-Mead bounds are not monotonic. They are being rearranged.');
                    new_bnds = bnds;
                    new_bnds(bnds(:,1) > bnds(:,2),1) = bnds(bnds(:,1) > bnds(:,2),2);
                    new_bnds(bnds(:,1) > bnds(:,2),2) = bnds(bnds(:,1) > bnds(:,2),1);

                    bnds = new_bnds;
                    bounds{1} = T_inv*bnds;
                end
            end
            
            population(:,i+1) = x_try;
            cost(i+1) = f_bnded(x_try);
        end 
        f_calls(1) = f_calls(1) + n+1;
    else
        %If we aren't coming up with population values, then assume that
        %the x0 is the initial population.
        population = x0;
        for i=1:n+1
            cost(i) = f_bnded(population(:,i));
        end
        f_calls(1) = f_calls(1) + 1;
    end
    x_dist = 1+x_stall;
    f_dist = 1+f_stall;
    
    loop = 1;
    if to_print
        fprintf('***** BEGINNING NELDER-MEAD SIMPLEX MINIMIZER ***** \n');
    end
    
    while loop < max_it && x_dist > x_stall && f_dist > f_stall
        if mod(loop,max_it/20)==1 && to_print == 2
            fprintf('Beginning loop %i.\n',loop)
        end
        
        eval_timer = tic();
        
        %Remove the worst point
        [~,inds] = sort(cost);
        
        %Find which point to replace and which to keep
        to_replace = inds(end);
        pop_best = inds(1);
        to_keep = inds(1:end-1);
        
        %Calculate the centroid of the points
        centroid = mean(population(:,to_keep),2);
        centroids(:,loop) = centroid;
        x_best(:,loop) = population(:,pop_best);
        
        x_try = centroid + alpha*(centroid-population(:,to_replace));
        f_try = f_bnded(x_try);
%         if isnan(f_try)
%             disp('isnan')
%             keyboard
%         end

        if f_try < cost(to_replace) && f_try > cost(pop_best)
            %Better than the worst, but worse than the best.
            %The reflection is acceptable.
            population(:,to_replace) = x_try;
            cost(to_replace) = f_try;
            f_calls(loop) = 1;
        elseif f_try <= cost(pop_best)
            %This is the best point so far, so expand the simplex
            x_expand = centroid + gamma*(centroid-population(:,to_replace));
            f_expand = f_bnded(x_expand);
            
            if f_expand < f_try
                %The expanded point is in fact better than our first
                %improved point, so keep this point now. We will not expand
                %further.
                population(:,to_replace) = x_expand;
                cost(to_replace) = f_expand;
            else
                %The expanded point was not better than our first improved
                %point, so keep the first point.
                population(:,to_replace) = x_try;
                cost(to_replace) = f_try;
            end
            f_calls(loop) = 2;
        elseif f_try >= cost(to_replace) || isnan(f_try)
            %We found a worse or equal point, so contract the simplex.
            x_contract = centroid + rho*(centroid-population(:,to_replace));
            f_contract = f_bnded(x_contract);
            
            if f_contract < cost(to_replace)
                %This point is better than our worst point, so keep it
                %because it generally improves the simplex.
                population(:,to_replace) = x_contract;
                cost(to_replace) = f_contract;
                f_calls(loop) = 2;
            else
                %If this point is still worse than all of our other points,
                %contract the entire simplex around the best point. 
                disp('contract entire pop')
                population(:,inds(2:end)) = repmat(population(:,pop_best),1,n) + sigma*(population(:,inds(2:end))-repmat(population(:,pop_best),1,n));
                for j=1:n+1
                    cost(j) = f_bnded(population(:,j));
                end
                f_calls(loop) = n+1+2;
            end
        else
            error('No action taken on simplex.')
        end
        f_avg(loop) = mean(cost);
        f_best(loop) = min(cost);
        f_worst(loop) = max(cost);
        
        eval_time(loop) = toc(eval_timer);
        
        conv_timer = tic();
        
        x_dist = 0;
        f_dist = 0;
        f_smallest = max(min(abs(cost)),1e-4);
        %Check for convergence
        for j=2:n+1
            for k=1:j-1
                x_dist = max(x_dist, norm(population(:,j)-population(:,k)));
                f_dist = max(f_dist, norm(cost(j)-cost(k)));
            end
        end
        
        x_dist = x_dist/n;
        f_dist = f_dist/f_smallest;
        
        check_conv_time(loop) = toc(conv_timer);
        
        max_x_dist(loop) = x_dist;
        max_f_dist(loop) = f_dist;
        
        loop = loop + 1;
    end
    
    %Why did we stop?
    if to_print
        if loop >= max_it
            fprintf('Algorithm halted; no solution found. %i iteration(s).\n',loop-1);
        elseif x_dist <= x_stall
            fprintf('Algorithm halted because the simplex is too small. Solution reached in %i iteration(s).\n',loop-1);
        elseif f_dist <= f_stall
            fprintf('Algorithm halted better solutions cannot be found. Solution reached in %i iteration(s).\n',loop-1);
        else
            fprintf('Algorithm halted. %i iteration(s).\n',loop-1);
        end
    end
    
    loop = loop - 1;
    
    if loop == max_it
        %For the sake of plotting, fix the index.
        loop = loop - 1;
    end
    
    f_avg = f_avg(1:loop);
    f_best = f_best(1:loop);
    f_worst = f_worst(1:loop);
    x_best = x_best(:,1:loop);
    
    centroids = centroids(:,1:loop);
    
    f_calls = f_calls(1:loop);
    max_x_dist = max_x_dist(1:loop);
    max_f_dist = max_f_dist(1:loop);  
    
    check_conv_time = check_conv_time(1:loop);
    eval_time = eval_time(1:loop);
    
        
    if to_plot
        figure 
        clf
        
        if all(f_best > 0) && all(f_avg > 0) && all(f_worst > 0)
            semilogy(f_best,'green','Linewidth',2)
            hold on
            semilogy(f_avg,'blue','Linewidth',2)
            semilogy(f_worst,'red','Linewidth',2)
            legend('Best','Mean','Worst')
        else
            hold on
            plot(f_best,'green','Linewidth',2)
            plot(f_avg,'blue','Linewidth',2)
            plot(f_worst,'red','Linewidth',2)
            legend('Best','Mean','Worst')
        end
        grid on
        
        hold on
        title('Average Population Error')
        xlabel('Iteration')
        hold off
        
        figure 
        clf
        hold on
        stairs(f_calls,'Linewidth',2)
        grid on
        title(['Function Calls by Iteration (',num2str(sum(f_calls)),' total)'])
        xlabel('Iteration')
        hold off
        
        figure 
        clf
        if all(max_x_dist > 0)
            semilogy(max_x_dist,'Linewidth',2)
            hold on
        else
            hold on
            plot(max_x_dist,'Linewidth',2)
        end
        grid on
        title('Maximum Simplex Chord')
        xlabel('Iteration')
        hold off
        
        figure 
        clf
        if all(max_f_dist > 0)
            semilogy(max_f_dist,'Linewidth',2)
            hold on
        else
            hold on
            plot(max_f_dist,'Linewidth',2)
        end
        grid on
        title('Largest Function Value Ratio within Simplex')
        xlabel('Iteration')
        hold off
        
        figure 
        clf
        hold on
        grid on
        plot(eval_time,'Linewidth',2)
        plot(eval_time+check_conv_time,'g','Linewidth',2)
        sum([sum(eval_time),sum(check_conv_time)]);
        title_str = ['Elapsed Time (total = ',num2str(sum([sum(eval_time),sum(check_conv_time)])),'s)'];
        title(title_str)
        ylabel('s')
        xlabel('Iteration')
        
        %If the problem is one or two dimensional, plot the contour and path taken
        %by the algorithm
        if n == 2
            centroids_scaled = T_inv*centroids;
            x_best_scaled = T_inv*x_best;
            %We can plot xs and fs
            figure
            clf
            hold on
            title('x and f by Iteration')
            xlabel(['x_',num2str(nb_var(1))])
            ylabel(['x_',num2str(nb_var(2))])
            zlabel('f')

            plot3(centroids_scaled(1,:),centroids_scaled(2,:),f_avg,'cyan*-','Linewidth',1);
            plot3(centroids_scaled(1,1),centroids_scaled(2,1),f_avg(1),'blue+','Linewidth',2);
            plot3(centroids_scaled(1,end),centroids_scaled(2,end),f_avg(1),'ro','Linewidth',4);

            plot3(x_best_scaled(1,:),x_best_scaled(2,:),f_avg,'g*-','Linewidth',2);
            plot3(x_best_scaled(1,1),x_best_scaled(2,1),f_avg(1),'blue+','Linewidth',2);
            plot3(x_best_scaled(1,end),x_best_scaled(2,end),f_avg(1),'ro','Linewidth',4);
            
            max_x1 = 1.5*max(abs(centroids(1,:)));
            max_x2 = 1.5*max(abs(centroids(2,:)));
            [X,Y] = meshgrid(linspace(-max_x1,max_x1,20),linspace(-max_x2,max_x2,20));
            
            max_x1_scaled = 1.5*max(abs(centroids_scaled(1,:)));
            max_x2_scaled = 1.5*max(abs(centroids_scaled(2,:)));
            [X_scaled,Y_scaled] = meshgrid(linspace(-max_x1_scaled,max_x1_scaled,20),linspace(-max_x2_scaled,max_x2_scaled,20));
            
            F = NaN(1,400);
            
            for j=1:400
                ftemp = f_bnded([X(j);Y(j)]);
                if ~isinf(ftemp)
                    F(j) = ftemp;
                end
            end
            F_sq = reshape(F,sqrt(length(F)),sqrt(length(F)));
            contour(X_scaled,Y_scaled,F_sq,100)
            hold off
        end
        
    end
    
    [~,inds] = sort(cost);
    
    %sort the entire population
    temp_pop = population;
    for i=1:n
        population(:,i) = temp_pop(:,inds(i));
    end

    %Return the lowest cost
    x_final = x_map(T_inv*population(:,1));
    f_final = min(cost);
    
    x_all = x_map(T_inv*population);
    
    stats = {f_avg,f_best,f_worst,max_x_dist,max_f_dist,f_calls,[eval_time,check_conv_time],loop,T_inv*centroids};
    
    if to_print        
        fprintf('***** NELDER-MEAD SIMPLEX MINIMIZER PERFORMANCE SUMMARY *****\n')
        fprintf('Iterations: %i\n',loop)
        fprintf('Function calls: %i\n',sum(f_calls))
        fprintf('Final point: [')
        for i=1:length(x_final)
            fprintf('%2.4f',x_final(i))
            if i < length(x_final)
                fprintf(', ')
            end            
        end
        fprintf('] \n')
        
        fprintf('\n***** NELDER-MEAD SIMPLEX MINIMIZER COMPLETE *****\n\n')
    end
    
end

