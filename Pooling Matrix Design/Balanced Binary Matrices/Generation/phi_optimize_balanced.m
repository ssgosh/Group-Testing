function B = phi_optimize_balanced(B, C)
    % B is an m x n balanced binary matrix and C is an n x n binary side
    % information matrix whose (i,j)th entry is 1 if there is a contact
    % between individuals i and j, and 0 otherwise. This function preserves
    % the psi value of the matrix while it optimizes on phi value.

    tic

    % Number of columns in B
    [~,n] = size(B);
    
    % Column weight
    c = sum(B(:,1));
    
    % Compute the gram matrix
    G = B' * B;
    
    % Total number of iterations
    NumIter = 1000000;

    % Print the initial Metric Value
    GMet = sum(G(C), 'all');
    fprintf("Iteration 0: %d\n", GMet);
    
    % The main optimization loop
    for i=1:NumIter
        x = randi(n);
        y = randi(n);
        while x == y
            x = randi(n);
            y = randi(n);
        end
        Gx = G(:,x);
        Gy = G(:,y);
        Cx = C(:,x);
        Cy = C(:,y);
        DeltaMet = sum(Gx(Cy),'all') + sum(Gy(Cx),'all') - sum(Gx(Cx),'all') - sum(Gy(Cy),'all');
        if C(x,y) == 1
           DeltaMet = DeltaMet - 2*c + 2*G(x,y);
        end
        
        if DeltaMet < 0
            GMet = GMet + 2*DeltaMet;
            B(:,[x y]) = B(:,[y x]);
            G(:,[x y]) = G(:,[y x]);
            G([x y],:) = G([y x],:);
        end
        
        if mod(i, 1000) == 0
            fprintf("Iteration %d: %d\n",i,GMet);
        end

        if GMet == 0
            break;
        end
    end
    
    t = toc;
    
    fprintf("Found an optimal matrix.\n");
    fprintf("Time taken: %.4f seconds.\n",t);
    fprintf("Total number of iterations: %d iterations.\n",i);