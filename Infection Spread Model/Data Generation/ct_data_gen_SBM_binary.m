clear;
clc;
close all;

% Seed for the random number generator
% seed = 1234;
seed = randi(10000);
% seed = 6551;
rng(seed);

% Choice of hyperparameters
% Stochastic Block Model parameters
% Contact trace graph on each day is a stochastic block model with below
% parameters
N = 100; % Number of nodes of the graph
r = 25; % Number of clusters

% For N = 1000, these settings seem reasonable
% in_p = 0.01; % within-cluster edge probability
% out_p = 0.0001; % inter-cluster edge probability

% For N = 100, these settings seem reasonable
% Below is basically an Erdos-Renyi graph, since in_p and out_p are equal
% in_p = 0.001; % within-cluster edge probability
% out_p = in_p; % inter-cluster edge probability

in_p = 0.5; % within-cluster edge probability
out_p = 0.0006; % inter-cluster edge probability

K = 100;
f = 1;
k1 = 1;
k2 = 14;
alpha_t = 1;
alpha_s = 1;
alpha_v = 1;
lambda_coeff = 1;
lambda_t_denom = 5000;
lambda_s_denom = 50;
lambda = 1 / (lambda_coeff * lambda_t_denom^alpha_t * lambda_s_denom^alpha_s);
%lambda = 1/(1.5 * 5000^alpha_t * 50^alpha_s * 16384.5^alpha_v);
%lambda = 1/(6 * 5000^alpha_t * 50^alpha_s * 16384.5^alpha_v);
numOffDiagonal = 3;
outp = 0; % 1/5000;

% Variables
CT = cell(K, 1);
X = cell(K, 1);
Ls = cell(K, 1);



% Choose f people to be infected initially
ListInfected = randperm(N);
ListInfected = ListInfected(1:f);
% ListInfected = [ 1 ];

% Set infection dates to large values initially (outside the range of consideration)
InfectionDates = (K+1)*ones(N,1);
InfectionDates(ListInfected) = zeros(f,1);

% Choose viral loads uniformly at random in (1,32768)
ViralLoads = zeros(N,1);
ViralLoads(ListInfected) = 1 + 32767*rand(f,1);

% Simulate for K days
for k=1:K
    % Contact matrix for day k
    A = PPM(N, r, in_p, out_p);

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

        % More infections from this person is possible only if their viral
        % load is not zero. This covers for recoveries. In the infection
        % spread model based on viral loads, the level of contact
        % automatically became 0, hence this if was not needed.
        if ViralLoads(ind) ~= 0
            for j=1:N
                if (A(ind, j) == 1 && k1 <= k - InfectionDates(ind) && k - InfectionDates(ind) < k2 && ~ismember(j,ListInfected))
                    if (rand < 1 - exp(-lambda*L(ind, j)))
                    % Commented out above if statement to make infections
                    % deterministic
                    % if (rand < 1)
                        ListInfected = [ListInfected j];
                        ViralLoads(j) = 1 + 32767*rand;
                        InfectionDates(j) = k;
                    end
                end
            end
        end
        % Recovery only happens at the end of a day. Hence any contacts on
        % the day of recovery may still lead to new infections. This is
        % consistent with previous code, where the level of contact was
        % computed before updating viral loads
        if (InfectionDates(ind) + 14 <= k)
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

%subtotalpos = totalpos(leftind:maxind+25);
%fprintf('%d, %d\n',mean(subtotalpos),std(subtotalpos));

fprintf("Seed: %d, Max positives: %d\n", seed, max(totalpos));
plot(totalpos)
fname = sprintf( '../Data/ct_data_SBM_binary_SIR_N_%d_r_%d_in_p_%f_out_p_%f_K_%d_seed_%d.mat', N, r, in_p, out_p, K, seed);
fprintf("%s\n", fname);
save(fname);