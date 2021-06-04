clear;
clc;
close all;

% Seed for the random number generator
rng(7);

% Choice of hyperparameters
N = 1000;
K = 500;
f = 3;
k1 = 3;
k2 = 8;
alpha_t = 1;
alpha_s = 1;
alpha_v = 1;
lambda = 1/(1.5 * 5000^alpha_t * 50^alpha_s * 16384.5^alpha_v);
numOffDiagonal = 50;
outp = 1/5000;

% Variables
CT = cell(K,1);
X = cell(K,1);

% Clique size distribution
probabilities = [0.04 0.12 0.12 0.21 0.21 0.3/7.0 0.3/7.0 0.3/7.0 0.3/7.0 0.3/7.0 0.3/7.0 0.3/7.0];
numPeople = [1 2 3 4 5 6 7 8 9 10 11 12];

% Grouping individuals into cliques
ngroups = 0;
prgroups = [];
count = 0;
for i=1:N
    count1 = randsample(numPeople, 1, true, probabilities);
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
    for i=1:N
       for j=1:N
           if (A(i,j) == 1)
               T(i,j) = 1 + 9999*rand;
               S(i,j) = 1 + 99*rand;
               L(i,j) = T(i,j)^alpha_t * S(i,j)^alpha_s * max(ViralLoads(i),ViralLoads(j))^alpha_v;
           end
       end
    end

    % Recovery, infection due to contacts
    for i=1:length(ListInfected)
        ind = ListInfected(i);
        if (InfectionDates(ind) + 14 <= k)
           ViralLoads(ind) = 0;
        end
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
end

[~,maxind] = max(totalpos);
subtotalpos = totalpos(maxind-24:maxind+25);
fprintf('%d, %d\n',mean(subtotalpos),std(subtotalpos));

save('../Data/ct_data.mat');