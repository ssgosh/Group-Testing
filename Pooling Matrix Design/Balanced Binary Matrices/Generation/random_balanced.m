function A = random_balanced(m,n,r,c)
    % This function returns a random balanced binary matrix with the input
    % params.

    % Number of non-zero elements in A
    NN = n*c;

    % Initialize the matrix A
    A = spalloc(m,n,NN);
    i = 0;
    for k=1:n
       if i+c <= m
           A(i+1:i+c,k) = 1;
           i = i+c;
       else
          A(i+1:m,k) = 1;
          A(1:i+c-m,k) = 1;
          i = i+c-m;
       end
    end

    % Row-column indices of non-zero elements in A
    [ni,nj] = find(A);
    
    % Total number of iterations
    NumIter = m*n;

    % The main optimization loop
    for i=1:NumIter
        flag = 0;
        while flag == 0
            x = randi(NN);
            y = randi(NN);
            a1 = ni(x);
            a2 = ni(y);
            b1 = nj(x);
            b2 = nj(y);
            if A(a1,b2) == 0 && A(a2,b1) == 0
                flag = 1;
            end
        end
        ni(x) = a1;
        nj(x) = b2;
        ni(y) = a2;
        nj(y) = b1;
        A(a1,b1) = 0;
        A(a2,b2) = 0;
        A(a1,b2) = 1;
        A(a2,b1) = 1;
    end
end