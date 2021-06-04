function A = psi_optimize_balanced(m,n,r,c)
    % If m <= n, one may call this function with m,n and r,c interchanged
    % and transpose the output matrix to get the required matrix. This
    % leads to faster processing times in general. For example, if m = 300,
    % n = 100, r = 10, and c = 3, one may use the snippet:
    % B = psi_optimize_balanced(1000, 3000, 3, 10);
    % A = transpose(B);
    % This functions returns a psi optimized matrix with the input params.
    
    tic

    % Seed for the random number generator
    rng(1);

    % Number of non-zero elements in A
    NN = n*c;

    % Minimum possible value of the optimization metric
    MinMetric = n*c*(r+c-1);

    % Initialize the matrix A
    A = spalloc(m,n,NN);
    i = 0;
    for j=1:n
       if i+c <= m
           A(i+1:i+c,j) = 1;
           i = i+c;
       else
          A(i+1:m,j) = 1;
          A(1:i+c-m,j) = 1;
          i = i+c-m;
       end
    end

    % Compute the initial Gram matrix
    G = full(A'*A);

    % Compute the initial value of the optimization metric
    Metric = sum(G.^2,'all');

    % Row-column indices of non-zero elements in A
    [ni,nj] = find(A);

    % Total number of iterations
    NumIter = 10000000;

    % Print the initial value of the optimization metric
    fprintf("Iteration 0: %d\n",Metric);

    % Default mode of the algorithm where flips are made at random
    mode = 0;

    % Interval at which to print
    Interval = 1000;

    % Variable to keep track of the metric at regular intervals
    IntermediateMetricValue = Metric;

    % The main optimization loop
    for i=1:NumIter
        flag = 0;
        while flag == 0
            x = randi(NN);
            y = randi(NN);

            % When mode one is switched on, choose the pair to be flipped
            % so that at least one of the edges being removed is part of a
            % cycle of length four.
            if mode == 1
                Gs = G - diag(diag(G));
                [c1,c2] = find(Gs > 1,1);
                if size(c1) > 0
                    nnr1 = find(A(:,c1));
                    nnr2 = find(A(:,c2));
                    nnr = intersect(nnr1,nnr2);
                    for j=1:size(ni)
                        if ismember(ni(j), nnr) && (nj(j) == c1 || nj(j) == c2)
                            x = j;
                            break;
                        end
                    end
                end
            end

            a1 = ni(x);
            a2 = ni(y);
            b1 = nj(x);
            b2 = nj(y);
            if A(a1,b2) == 0 && A(a2,b1) == 0
                flag = 1;
            end
        end

        A1 = A(a1,:)';
        A2 = A(a2,:)';
        A1(b1) = 0;
        A2(b2) = 0;
        NNR1 = find(A1);
        NNR2 = find(A2);
        NNR = union_sorted(NNR1,NNR2);

        B1 = G(:,b1);
        B2 = G(:,b2);
        B1(NNR1) = B1(NNR1) - 1;
        B1(NNR2) = B1(NNR2) + 1;
        B2(NNR2) = B2(NNR2) - 1;
        B2(NNR1) = B2(NNR1) + 1;

        Old = 2 * (sum(G(NNR,b1).^2) + sum(G(NNR,b2).^2));
        New = 2 * (sum(B1(NNR).^2) + sum(B2(NNR).^2));
        DeltaMetric = New - Old;

        if DeltaMetric < 0
            ni(x) = a1;
            nj(x) = b2;
            ni(y) = a2;
            nj(y) = b1;
            A(a1,b1) = 0;
            A(a2,b2) = 0;
            A(a1,b2) = 1;
            A(a2,b1) = 1;
            G(b1,NNR) = B1(NNR)';
            G(b2,NNR) = B2(NNR)';
            G(NNR,b1) = B1(NNR);
            G(NNR,b2) = B2(NNR);
            Metric = Metric + DeltaMetric;
        end

        if Metric == MinMetric
            break;
        elseif mod(i,Interval) == 0
            fprintf("Iteration %d: %d\n",i,Metric);
            if Metric == IntermediateMetricValue
                mode = 1;
                G = sparse(G);
            end
            IntermediateMetricValue = Metric;
        elseif mode == 1
            fprintf("Iteration %d: %d\n",i,Metric);
        end
    end

    t = toc;

    fprintf("Found an optimal matrix.\n");
    fprintf("Time taken: %.4f seconds.\n",t);
    fprintf("Total number of iterations: %d iterations.\n",i);