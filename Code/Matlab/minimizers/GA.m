function [ x_final, f_final, x_all, f_all ] = GA( f, bnds, linear_constraints, conv_crit, pop_size, crossover_fraction, crossover_method, mutation_prob, seed, to_parallelize, to_plot, to_print )
%GA Summary of this function goes here
%
%   f - function to minimize. It cannot return NaN or Inf; this will cause
%   errors.
%
%   bnds - nx2 matrix of the lower and upper limits of each variable
%
%   conv_crit - four element vector holding information about when to
%   stop. The first element is the maximum number of generations. The
%   second element is the fraction of the population that must be
%   identical to halt the optimizer. The third element is the radius in the
%   parameter space that defines a ball of "convergence". The fourth
%   element is the objective function value that, once reached, halts the
%   algorithm.
%
%   pop_size - the number of individuals in the population
%
%   crossover_fraction - the number of individuals that are allowed to
%   crossover each iteration
%
%   crossover_method - string 'single', 'double', 'random', or
%   'continuous' describing the method of crossover for mating individuals.
%
%   mutation_prob - probability that an individual mutates
%
%   seed - random number seed to use
%
%   to_parallelize - boolean of whether to parallelize function evaluation
%   or not. 
%
%   to_plot - boolean of whether to plot results or not. If to_plot = 0,
%   nothing will be plotted. If to_plot = 1, the end results will be
%   plotted. If to_plot = 2, plots will be updated after every geenration.
%
%   to_print - boolean of whether to plot messages to the user or not.
    
    %Parse our input arrays to make sure the default values are applied if
    %information is missing.
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
    x0 = mean(bnds,2);
    allow_infeasible_start = true;
    [x_nonbasic_to_full_map,~,nb_var] = elim_var(A,b,x0,allow_infeasible_start);
    
    nb_bnds = bnds(nb_var,:);
    
%     f_con = @(x) f(x_nonbasic_to_full_map(x));
    f_con = @(x) objective_function(f,x_nonbasic_to_full_map,eye(length(nb_var)),x,1e3);

    %Parse our convergence criteria
    switch length(conv_crit)
        case 1
            max_generations = conv_crit(1);
            identical_pop_fraction = 0.95;
            proximity_radius = 1e-6;
            objective_threshold = Inf;
        case 2
            max_generations = conv_crit(1);
            identical_pop_fraction = conv_crit(2);
            proximity_radius = 1e-6;
            objective_threshold = Inf;
        case 3
            max_generations = conv_crit(1);
            identical_pop_fraction = conv_crit(2);
            proximity_radius = conv_crit(3);
            objective_threshold = Inf;
        case 4
            max_generations = conv_crit(1);
            identical_pop_fraction = conv_crit(2);
            proximity_radius = conv_crit(3);
            objective_threshold = conv_crit(4);
        otherwise
            max_generations = 1000;
            identical_pop_fraction = 0.95;
            proximity_radius = 1e-6;
            objective_threshold = Inf;
    end
    
    if isempty(pop_size)
        pop_size = 1000;
    end
    
    if isempty(crossover_fraction)
        crossover_fraction = .5;
    end
    
    if isempty(crossover_method)
        crossover_method = 'continuous';
    end
    
    if isempty(mutation_prob)
        mutation_prob = .001;
    end
    
    if isempty(seed)
    else
        %Set the random number generator seed
        rng(ceil(seed));
    end
    
    %Define n
    n = size(nb_bnds,1);
    
    %Limit the crossover fraction
    crossover_fraction = min(crossover_fraction,1);
    
    %Initialize the population
    current_population = NaN(n,pop_size);
    spaced_params = NaN(n,pop_size);
    for i=1:n
        %We want all of our properties to come from a grid, to aid
        %convergence. First, create as many unique values for each
        %parameter as we need, then save them in the population in a random
        %order.
        dim_range = nb_bnds(i,2) - nb_bnds(i,1);
        buffer = dim_range/(2*pop_size);
        spaced_params(i,:) = linspace(nb_bnds(i,1)+buffer,nb_bnds(i,2)-buffer,pop_size);
        current_population(i,:) = spaced_params(i,randperm(pop_size));
    end
    
    %Initialize the fitness vector
    fitness = NaN(pop_size,max_generations);
    
    %Initialize the distance (from best solution) vector for convergence
    distance = NaN(ceil(identical_pop_fraction*pop_size),max_generations);
    
    %Initialize stopping criteria
    %max_it_met = 0;    
    %x_tol_met = 0;
    i = 1;
    
    %Function evaluations
    f_calls = 0;
    
    %Initialize time vectors
    eval_time = NaN(max_generations,1);
    crossover_time = NaN(max_generations,1);
    mutate_time = NaN(max_generations,1);
    check_conv_time = NaN(max_generations,1);
    
    if to_print
        disp('***** BEGINNING GENETIC ALGORITHM *****')
    end
    %Main loop
    converged = false;
    while ~converged
        if to_print == 2
            fprintf('Beginning iteration %i.\n',i);  
        end
        
        tic;
        %Calculate fitness. Since we assume we're minimizing, we compute
        %the fitness as the negative objective function. This turns the
        %problem into a maximization problem.
        
        %However, we do not need to calculate the fitness of every
        %individual; the top (1-crossover_fraction) will be the same. For
        %example, if 1/4 the pop can reproduce, the bottom 1/4 will be
        %replaced but the top 3/4 will be the same. This is aided by the
        %fact that we do not sort by fitness until after this step.       
        if i > 1
            num_to_evaluate = ceil((1 - crossover_fraction) * pop_size);
            inds_to_evaluate = union(inds(1:num_to_evaluate),individuals_mutated);
            
            num_to_evaluate = length(inds_to_evaluate);
            %Copy fitnesses from the previous generation, but overwrite
            %those we need to evaluate with NaN.
            fitness(:,i) = fitness(:,i-1);
            fitness(inds_to_evaluate,i) = NaN;
            
            %Also keep track of those indices we've already evaluated. If
            %we've already done so, don't evaluate the function. 
            inds_already_evaluated = (1:pop_size)';
            inds_already_evaluated(inds_to_evaluate) = [];
            num_already_evaluated = length(inds_already_evaluated);
        else
            %For the first iteration, we need to evaluate every individual.
            num_to_evaluate = pop_size;
            num_already_evaluated = 0;
            inds_to_evaluate = 1:pop_size;
            inds_already_evaluated = [];
        end
        
        pop_to_evaluate = current_population(:,inds_to_evaluate);
        pop_already_evaluated = current_population(:,inds_already_evaluated);
        
        sliced_fitness = zeros(num_to_evaluate,1);
       
        if to_parallelize
            %First, find the member of the population that we need to
            %evaluate, and store that as a boolean vector need_to_evaluate.
            %Assume we need to evaluate all to start.
            need_to_evaluate = ones(num_to_evaluate,1);
            for j=1:num_to_evaluate
                k = 0;
                %Compare this individual with those that we've already
                %evaluated. If there is a repeat, then simply assign the
                %same fitness value to the fitness vector.
                while (k < num_already_evaluated) && need_to_evaluate(j)
                    k = k + 1;
                    if isequal(pop_to_evaluate(:,j),pop_already_evaluated(:,k))
                        need_to_evaluate(j) = false;
                        sliced_fitness(j) = fitness(inds_already_evaluated(k),i);
                    end
                end
            end
            %Now, scan through the original list of individuals to
            %evaluate, and skip those that we've already evaluated. Perform
            %the evaluations in a parfor loop.
            parfor j=1:num_to_evaluate
                if need_to_evaluate(j)
                    sliced_fitness(j) = -f_con(pop_to_evaluate(:,j));
                    f_calls = f_calls + 1;
                end
            end
        else
            %For the number of individuals to evaluate...
            for j=1:num_to_evaluate
                k = 0;
                need_to_evaluate = true;
                %...compare that individual against all of the individuals
                %that we've already evaluated. If it is a duplicate, set
                %the "need_to_evaluate" boolean to false.
                while k < num_already_evaluated && need_to_evaluate
                    k = k + 1;
                    if isequal(pop_to_evaluate(:,j),pop_already_evaluated(:,k))
                        need_to_evaluate = false;
                    end
                end
                %Now, only evaluate the individual if it is novel.
                if need_to_evaluate
                    sliced_fitness(j) = -f_con(pop_to_evaluate(:,j));
                    f_calls = f_calls + 1;
                else
                    sliced_fitness(j) = fitness(inds_already_evaluated(k),i);
                end
            end
        end

        fitness(inds_to_evaluate,i) = sliced_fitness;
        
        %Sort the population based on fitness for crossover. This will put
        %the highest values (the most fit) at the bottom.
        [vals,inds] = sort(fitness(:,i));
        %++++++++++++++++++++++++++++++++++++++++++
        
        if any(isnan(fitness(:,i)))
            keyboard
            error('One or more fitnesses are NaN (GA.m).')
        end
        
        eval_time(i) = toc();
        conv_time_start = tic;
        
        if any(fitness(:,i) >= objective_threshold)
            converged = true;
            disp('GA stopped because maximum objective was reached.')
        end
        
        %Check convergence based on x values
        j_lim = ceil(identical_pop_fraction*pop_size);
        j_lim = min(pop_size-2,j_lim);
        for j=1:j_lim
            %Use the 2 norm to find the largest variation of an
            %individual from the current best solution
            distance(j,i) = norm(current_population(:,inds(pop_size-j-1)) - current_population(:,inds(end)),2);
        end
        %Compare these values, for each variable, to the convergence
        %tolerance specified by the user. 
        x_tol_met = all((distance(1:j_lim,i) < proximity_radius));
        
        if x_tol_met
            converged = true;
            if to_print
                disp('Solution stalled because x stopped changing.')
            end
        end
        
        if i >= max_generations
            converged = true;
            if to_print
                disp('Maximum generations reached.')
            end
        end
        
        check_conv_time(i) = toc(conv_time_start);
        tic;
        
        if ~converged
            %Create a CDF for mating probabilities. 
            normalized_vals = cumsum(vals/sum(vals));

            %As time goes on, we want to favor the most successful solutions
            %more. Presumably they are closer to a solution, so we weight them
            %with a sigmoid. Initially it is flat and all solutions are favored
            %based on their fitness. As the algorithm progresses, more
            %successful solutions' fitness is doubled while less successful
            %solutions' fitness approaches 0.
            normalized_vals = normalized_vals.* (1 + min(2*i/max_generations,1)*erf(.03*((1:pop_size)-pop_size/2)))';
            normalized_vals = max(0,cumsum(normalized_vals/sum(normalized_vals)));

            %Generate random numbers to use for mating combinations
            mate_probability = rand(ceil(crossover_fraction*pop_size/2),2);

            %First, assign the current population to the next population. We
            %will now replace the least fit (counting upwards along inds[])
            %with the new individuals.
            next_population = current_population;

            %Determine parents and create new individuals
        
            for j=1:floor(crossover_fraction*pop_size/2)
                parent_a = find(normalized_vals > mate_probability(j,1),1,'first');
                parent_b = find(normalized_vals > mate_probability(j,2),1,'first');
                if strcmp(crossover_method,'single')
                    %If we are doing a single crossover, randomly find a boundary
                    %and exchange the genes on either chromosome after that index.
                    if n > 1
                        gene = randi(n-1);
                        try next_population(:,inds(2*j-1)) = [current_population(1:gene,inds(parent_a));current_population(gene+1:end,inds(parent_b))];
                        catch
                            error('Mating not performed properly; ensure the objective function does not return Inf or NaN.');
                        end
                        next_population(:,inds(2*j)  ) = [current_population(1:gene,inds(parent_b));current_population(gene+1:end,inds(parent_a))];
                    else
                        next_population(:,inds(2*j-1)) = current_population(inds(parent_a));
                        next_population(inds(2*j)  ) = current_population(inds(parent_b));
                    end
                elseif strcmp(crossover_method,'double')
                    %If we are doing a double crossover, use the above method with
                    %two boundaries, swapping a middle segment of the chromosomes.
                    gene = sort(ceil(rand(1,2)*(n-1)));
                    next_population(:,inds(2*j-1)) = [current_population(1:gene(1),inds(parent_a));current_population(gene(1)+1:gene(2),inds(parent_b));current_population(gene(2)+1:end,inds(parent_a))];
                    next_population(:,inds(2*j)  ) = [current_population(1:gene(1),inds(parent_b));current_population(gene(1)+1:gene(2),inds(parent_a));current_population(gene(2)+1:end,inds(parent_b))];

                elseif strcmp(crossover_method,'random')
                    %If we do random crossover, the elements to swap are
                    %randomly selected. A random number of random indices are
                    %selected and swapped.
                    gene = ceil(rand(ceil(rand*n),1)*n);
                    next_population(:,inds(2*j-1)) = current_population(:,inds(parent_a));
                    next_population(gene,inds(2*j-1)) = current_population(gene,inds(parent_b));
                    next_population(:,inds(2*j)) = current_population(:,inds(parent_b));
                    next_population(gene,inds(2*j)) = current_population(gene,inds(parent_a));

                elseif strcmp(crossover_method,'continuous')
                    %If we use a continuous crossover, then every element
                    %in the chromosome is a random combination of the two
                    %parents. Each crossover produces two offspring, one
                    %using gene and one using 1-gene.
                    gene = rand(n,1);
                    next_population(:,inds(2*j-1)) = gene.*current_population(:,inds(parent_a)) + (1-gene).*current_population(:,inds(parent_b));
                    next_population(:,inds(2*j  )) = gene.*current_population(:,inds(parent_b)) + (1-gene).*current_population(:,inds(parent_a));

                elseif strcmp(crossover_method,'continuous_AJH')
                    %If we use a continuous crossover, then a random number of
                    %elements
                    %in the chromosome is a random combination of the two
                    %parents. Each crossover produces two offspring, one
                    %using gene and one using 1-gene.
                    gene_nums = [1:n]';
                    %If an index is in genes_to_cross, then that parameter will
                    %be randomized in the children. Otherwise, it is copied
                    %from the parent.
                    genes_to_cross = gene_nums(rand(n,1) > .5+zeros(n,1));
                    gene = 1+zeros(n,1);

                    %If a value in gene equals 1, the child will be a copy of
                    %one parent for that index. If a value does not equal one,
                    %then the two children will have some complementary values
                    %for that index, bounded by the parents' values for that
                    %index.
                    gene(genes_to_cross) = rand(size(genes_to_cross));            

                    next_population(:,inds(2*j-1)) = gene.*current_population(:,inds(parent_a)) + (1-gene).*current_population(:,inds(parent_b));
                    next_population(:,inds(2*j  )) = gene.*current_population(:,inds(parent_b)) + (1-gene).*current_population(:,inds(parent_a));

                else
                    error('A valid crossover method was not entered.')
                end
            end
        end
        
        crossover_time(i) = toc();
        tic;
        
        %Compute mutations. Mutations are not common in GAs, but do help to
        %keep the solution from stagnating. We calculate the number of
        %mutations, as a fraction of the total number of properties in the
        %population. We then generate random properties to mutate, which
        %are changed randomly within their bounds.
        children_inds = inds(1:floor(crossover_fraction*pop_size));
        num_children = length(children_inds);
        
        to_mutate = (rand(n,num_children) < mutation_prob);
        %Find all the indices that we need to mutate.
        inds_to_mutate = find(to_mutate);

        properties_mutated = mod(inds_to_mutate-1,n)+1;
        %Whichever individuals are mutated will need to be reevaluated.
        %Keep track of these (the columns that were mutated), so we can
        %evaluate them in our evaluation stage.
        individuals_mutated = inds(ceil(inds_to_mutate/n));
        
        for j=1:length(inds_to_mutate)
            next_population(properties_mutated(j),individuals_mutated(j)) = spaced_params(properties_mutated(j),randi(pop_size));
%             next_population(inds_to_mutate(j)) = spaced_params(properties_mutated(j),randi(pop_size));
        end
        
%         num_mutations = ceil(mutation_prob*n*pop_size);
%         inds_to_mutate = randi(n*pop_size,num_mutations,1);
%         properties_mutated = mod(inds_to_mutate-1,n)+1;
%         %Whichever individuals are mutated will need to be reevaluated.
%         %Keep track of these (the columns that were mutated), so we can
%         %evaluate them in our evaluation stage.
%         individuals_mutated = ceil(inds_to_mutate/n)  
%         for j=1:length(inds_to_mutate)
%             next_population(inds_to_mutate(j)) = spaced_params(properties_mutated(j),randi(pop_size));
%         end
        
        mutate_time(i) = toc();
       
        
        if to_plot == 2
            if i==1
                fit_handle = figure;
                maxfit_handle = figure;
                time_handle = figure;
            end
            
            figure(fit_handle)
            hold on
            title('Population fitness')
            plot(i,mean(fitness(:,i)),'b.','Linewidth',5)
            plot(i,min(fitness(:,i)),'r.')
            plot(i,max(fitness(:,i)),'g.')
            legend('Mean','Min','Max')
            drawnow
            hold off
            
            figure(maxfit_handle)
            hold on
            plot(i,max(fitness(:,i)),'g.');
            title('Fittest individual') 
            drawnow
            hold off 
            
            figure(time_handle)
            hold on
            plot(i,eval_time(i),'b.','Linewidth',2)
            plot(i,crossover_time(i)+eval_time(i),'g.','Linewidth',2)
            plot(i,mutate_time(i)+crossover_time(i)+eval_time(i),'r.','Linewidth',2)
            plot(i,check_conv_time(i)+mutate_time(i)+crossover_time(i)+eval_time(i),'cyan.','Linewidth',2)
            title_str = ['Stacked line graph of runtime (total = ',num2str(sum(eval_time(1:i))+sum(crossover_time(1:i))+sum(mutate_time(1:i))+sum(check_conv_time(1:i))),'s)'];
            title(title_str);
            legend('Function eval','Crossover','Mutation','Check convergence')
            drawnow
            hold off
        end
        if to_print == 2
            fprintf('%i function calls so far.\n',f_calls);  
        end
        
        i = i + 1;
        
        %We only want to overwrite the current population if we haven't
        %converged. If we have converged, then we need to keep the current
        %population so we output the correct data.
        if ~converged
            current_population = next_population;
        end
    end

    max_ind = find(isnan(fitness(1,:)),1,'first')-1;
    if isempty(max_ind)
        max_ind = max_generations;
    end
    fitness = fitness(:,1:max_ind);

    if to_plot == 1
        figure
        clf
        hold on
        title('Population fitness')
        for j=1:max_ind
            plot(j,mean(fitness(:,j)),'b.','Linewidth',5)
            plot(j,min(fitness(:,j)),'r.')
            plot(j,max(fitness(:,j)),'g.')
        end
        legend('Mean','Min','Max')
        drawnow
        hold off

        figure
        clf
        hold on
        max_fit = NaN(max_ind,1);
        for j=1:max_ind
            max_fit(j) = max(fitness(:,j));
        end
        semilogy(1:max_ind,max_fit,'g.');
        drawnow
        hold off
        title('Fittest individual')    

        figure
        clf
        hold on
        plot(eval_time(1:max_ind),'Linewidth',2)
        plot(crossover_time(1:max_ind)+eval_time(1:max_ind),'g','Linewidth',2)
        plot(mutate_time(1:max_ind)+crossover_time(1:max_ind)+eval_time(1:max_ind),'r','Linewidth',2)
        plot(check_conv_time(1:max_ind)+mutate_time(1:max_ind)+crossover_time(1:max_ind)+eval_time(1:max_ind),'cyan','Linewidth',2)
        title_str = ['Stacked line graph of runtime (total = ',num2str(sum(eval_time(1:max_ind))+sum(crossover_time(1:max_ind))+sum(mutate_time(1:max_ind))+sum(check_conv_time(1:max_ind))),'s)'];
        title(title_str);
        ylim([0,1.1*max(check_conv_time(1:max_ind)+mutate_time(1:max_ind)+crossover_time(1:max_ind)+eval_time(1:max_ind))])
        legend('Function eval','Crossover','Mutation','Check convergence')
        drawnow
        hold off
    end
    
    x_final = x_nonbasic_to_full_map(current_population(:,inds(end)));
    f_final = -f(x_final);
    
    x_all = x_nonbasic_to_full_map(current_population(:,inds));
    f_all = -fitness(inds,end);
    
    if to_print
        disp(['The solution was reached in ',num2str(i-1),' generations and ',num2str(f_calls),' function evaluations,'])
        disp(['a ',num2str(100*(1-f_calls/(pop_size*(i-1)))),'% reduction from expected.'])
        fprintf('***** GENETIC ALGORITHM COMPLETE *****\n\n')
    end
end

