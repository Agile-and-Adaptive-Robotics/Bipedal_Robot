function [ nonbasic_to_full_map, x0_nonbasic, nonbasic_var, basic_var ] = elim_var( A, b, x0, allow_infeasible_start )
%ELIM_VAR Linear equality constraints can be added to nonlinear
%minimization problems by eliminating variables. This function returns a
%function that maps the reduced set of nonbasic variables back into the
%full independent variable vector
%   Detailed explanation goes here
    
    m = size(b,1);
    n = size(x0,1);
    basic_var = [];
    
    %First ensure that A has full row rank
    if m > 0 && size(A,1) < size(A,2) && size(A,1) == m && rank(A) == m
        %Our A matrix is of the correct dimensions and rank
%         if all(A*x0-b == zeros(m,1))
        if allow_infeasible_start || all(abs(A*x0-b) <= 1e-3)
            %Our initial point is feasible.
        
            %Separate A into two matrices. A will be rectangular, m x n. We
            %desired A to be rearranged into a full ran m x m matrix, and leftover
            %columns.
            all_combs_pre = combnk(n:-1:1,m);
            all_combs = all_combs_pre(:,m:-1:1);

            unmet = 1;
            i = 1;
            while unmet && i <= size(all_combs,1)
                %Split up A, test rank of left square sub-matrix
                B_try = A(:,all_combs(i,:));
                rank_B = rank(A(:,all_combs(i,:)));

                if rank_B == m
                    %We have sufficiently divided the A matrix
                    %Find the basic and nonbasic variables, as well as B and N.
                    B = B_try;
                    basic_var = all_combs(i,:);
                    nonbasic_var = setdiff(1:n,basic_var);
                    N = A(:,nonbasic_var);

                    %Our original, un-eliminated x vector will now be in a
                    %different order; store that mapping
                    reorder = [basic_var,nonbasic_var];

                    unmet = 0;
                end

                %Increment the counter
                i = i + 1;
            end

            %Use the reordering of x by basic and nonbasic variables to create
            %a permutation matrix to map the nonbasic variables back into the
            %full x vector, in the proper order.
            perm_mat = zeros(n);
            for i=1:n
                perm_mat(reorder(i),i) = 1;
            end

            %We need to solve for some terms in our transform function. See
            %page 87 in Griva, Nash and Sofer.
            B_inv_b = linsolve(B,b);
            B_inv_N = linsolve(B,N);

            %Output the function that will take in nonbasic variables only, and
            %return the full x vector.
            nonbasic_to_full_map = @(x)perm_mat*([-B_inv_N;eye(n-m)]*x+repmat([B_inv_b;zeros(n-m,1)],1,size(x,2)));
            
            x0_nonbasic = x0(nonbasic_var);
            nonbasic_var = nonbasic_var';
        else
            warning('Initial point is infeasible. No constraints will be applied.')
            nonbasic_to_full_map = @(x)x;
            x0_nonbasic = x0;
            nonbasic_var = (1:n)';
        end
    else
        if m == 0
            %Do nothing
        elseif rank(A) ~= m
            warning('Redundant constraints are present. No constraints will be applied.')
        elseif size(A,1) >= size(A,2)
            warning('Problem is overconstrained. No constraints will be applied.')
        elseif size(A,1) ~= m
            warning('Constraint equation is not the proper size. No constraints will be applied.')
        end
        nonbasic_to_full_map = @(x)x;
        x0_nonbasic = x0;
        nonbasic_var = (1:n)';
    end

end

