function fname = ct_data_gen_binary_fn(num)
    % Seed for the random number generator
    seed = randi(10000);
    rng(seed);
    
    % Choice of hyperparameters
    % N = 1000;  
    N = 100;
    K = 100;
    f = 1;
    k1 = 1;
    k2 = 14;
    k3 = 14;
    alpha_t = 1;
    alpha_s = 1;
    lambda_coeff = 4;
    lambda_t_denom = 5000;
    lambda_s_denom = 50;
    lambda = 1 / (lambda_coeff * lambda_t_denom^alpha_t * lambda_s_denom^alpha_s);
    numOffDiagonal = 3;
    outp = 0; %1/5000;
    
    % Variables
    CT = cell(K, 1);
    X = cell(K, 1);
    Ls = cell(K, 1);
    
    % Clique size distribution
    probabilities = [0.04 0.12 0.12 0.21 0.21 0.3/7.0 0.3/7.0 0.3/7.0 0.3/7.0 0.3/7.0 0.3/7.0 0.3/7.0];
    numPeople = [1 2 3 4 5 6 7 8 9 10 11 12];
    %probabilities = [1.0];
    %numPeople = [20];
    
    % Grouping individuals into cliques
    ngroups = 0;
    prgroups = [];
    count = 0;
    for i=1:N
        count1 = randsample(numPeople, 1, true, probabilities);
        %count1 = 20;
        ngroups = ngroups + 1;
        if count + count1 >= N
            prgroups = [prgroups [count + 1 N]'];
            break;
        else
            prgroups = [prgroups [count + 1 count + count1]'];
            count = count + count1;
        end   
    end
    
    % Block Diagonal Base matrix
    BlockDiagonalBase = zeros(N,N);
    for i=1:ngroups
        groupSize = prgroups(2,i) - prgroups(1,i) + 1;
        BlockDiagonalBase(prgroups(1,i):prgroups(2,i),prgroups(1,i):prgroups(2,i)) = ones(groupSize, groupSize) - eye(groupSize);
    end
    
    % Choose f people to be infected initially
    ListInfected = randperm(N);
    ListInfected = ListInfected(1:f);
    
    % Set infection dates to large values initially (outside the range of consideration)
    InfectionDates = (K+1)*ones(N,1);
    InfectionDates(ListInfected) = zeros(f,1);
    
    % Choose viral loads uniformly at random in (1,32768)
    ViralLoads = zeros(N,1);
    ViralLoads(ListInfected) = 1 + 32767*rand(f,1);
    
    % Simulate for K days
    for k=1:K
        % Contact matrix for day k
        A = BlockDiagonalBase;
        count = 0;
        while count < numOffDiagonal
            i = randsample(N,1);
            j = randsample(N,1);
            if i ~= j && A(i,j) ~= 1
                A(i,j) = 1;
                A(j,i) = 1;
                count = count + 1;
            end
        end
    
        % Duration of contact, signal strength and level of contact matrices (each with basis corresponding to ones in A)
        T = zeros(N,N);
        S = zeros(N,N);
        L = zeros(N,N);
        %for i=1:N
           %for j=1:N
         for i=1:(N-1)
           for j=(i+1):N
               if (A(i,j) == 1)
                   T(i,j) = 1 + 9999*rand;
                   S(i,j) = 1 + 99*rand;
                   L(i,j) = T(i,j)^alpha_t * S(i,j)^alpha_s;
                   % Make the T, S and L matrices symmetric
                   T(j, i) = T(i, j);
                   S(j, i) = S(i, j);
                   L(j, i) = L(i, j);
               end
           end
        end
    
        % Recovery, infection due to contacts
        for i=1:length(ListInfected)
            ind = ListInfected(i);
            % More infections happen only if viral load is non-zero
            % Earlier code was also correct since we also check for
            % infection date
            if ViralLoads(ind) ~= 0
                for j=1:N
                    if (A(ind, j) == 1 && k1 <= k - InfectionDates(ind) && k - InfectionDates(ind) < k2 && ~ismember(j,ListInfected))
                        if (rand < 1 - exp(-lambda*L(ind, j)))
                            ListInfected = [ListInfected j];
                            ViralLoads(j) = 1 + 32767*rand;
                            InfectionDates(j) = k;
                        end
                    end
                end
            end
            % Set viral load to 0 if recovered
            % If k2 is equal to 14, then it is assumed that recovery
            % happens at the end of the day, hence this individual can
            % infect another person at the 14th day
            if (InfectionDates(ind) + k3 <= k)
               ViralLoads(ind) = 0;
            end            
        end
    
        % Infections from outside contact
        for i=1:N
            if ~ismember(i,ListInfected)
                if rand < outp
                    ListInfected = [ListInfected i];
                    ViralLoads(i) = 1 + 32767*rand;
                    InfectionDates(i) = k;
                end
            end
        end
    
        totalpos(k) = nnz(ViralLoads);
        
        X{k} = ViralLoads;
        CT{k} = sparse(A);
        Ls{k} = sparse(L);
    end
    
    [~,maxind] = max(totalpos);
    fprintf('maxind: %d \n', maxind);
    leftind = max(1, maxind - 24);
    fprintf('leftind: %d \n', leftind);
    
    % subtotalpos = totalpos(leftind:maxind+25);
    % fprintf('%d, %d\n',mean(subtotalpos),std(subtotalpos));
    
    fprintf("Seed: %d, Max positives: %d\n", seed, max(totalpos));
    
%     plot(totalpos);
    
    fname = sprintf(['../Data/SIR/%d/%03d_ct_data_binary_family_structured_'...
        'N_%d_K_%d_'...
        'lambda_factor_%d_'...
        'numOffDiagonal_%d_' ...
        'f_%d_k1_%d_k2_%d_outp_%f_' ...
        'seed_%d.mat'], ...
        N, num,  N, K, lambda_coeff, numOffDiagonal, f, k1, k2, outp, seed);
    fprintf("%s\n", fname);
    save(fname);
end