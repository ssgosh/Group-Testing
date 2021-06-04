clear;
clc;
close all;

addpath('Helper Functions');

rng(1); % seed for the random number generator

nsig    = 100;  % number of signals
s       = 40;   % sparsity
sigval  = 0.1;  % standard deviation of multiplicative noise
tau     = 0.2;  % prediction threshold
q       = 0.95; % RT-PCR cycle amplification factor

for ss=1:nsig
    % Read pooling matrix from file
    load('../Pooling Matrix Design/Balanced Binary Matrices/Examples/psi_opt_bal_300x1000.mat', 'B');
    A = full(B);
    [m,n] = size(A);
    
    % Signal generation
    x = zeros(n,1);
    indices = randperm(n);
    x(indices(1:s)) = 1+32767*rand(s,1)'; 

    % Generate noisy measurement vector
    eps = sigval*randn(m,1);
    y = (A*x).*((1+q).^eps);

    % COMP algorithm
    final_x = ones(n,1);
    nz_indices_y = find(y>0);       % nonzero measurements
    z_indices_y = find(y==0);       % zero-valued measurements
    for j=1:length(z_indices_y)
       zi = z_indices_y(j);
       final_x(find(A(zi,:)>0)) = 0; 
    end
    A = A(:,final_x>0);                 % remove entries of x which are obvious negatives
    nz_indices_x = find(final_x ~=  0); % store indices of putative positives
    y = y(nz_indices_y);                % remove entries of y corresponding to zero-valued measurements
    A = A(nz_indices_y,:);              % remove rows of A corresponding to zero-valued measurements
    [m,n] = size(A);                    % updated size of A
    fprintf('\nSignal %d: %d %d;', ss, m, n);

    % COMP result metrics
    rmseval_comp(ss) = norm(x-final_x,2)/norm(x,2);
    fn_comp(ss) = length(find(x > 0 & final_x == 0));
    fp_comp(ss) = length(find(x == 0 & final_x > 0));
    sen_comp(ss) = length(find(x > 0 & final_x > 0))/length(find(x > 0)); 
    spec_comp(ss) = length(find(x == 0 & final_x == 0))/length(find(x == 0)); 
    fprintf (' %d %d;',fn_comp(ss),fp_comp(ss));

    % Stuff for cross-validation
    mr = ceil(0.9*m); mcv = m-mr;
    indices = randperm(m);
    yr = y(indices(1:mr));
    ycv = y(indices(mr+1:end));
    Ar = A(indices(1:mr),:);
    Acv = A(indices(mr+1:end),:);
    lambda_min = max([sigval*sqrt(log(n))-8,0.001]);
    lambda_max = sigval*sqrt(log(n))+8;
    lambdas = lambda_min:0.05:lambda_max;
    l_lambdas = length(lambdas);

    % Cross Validation
    cv_error = length(l_lambdas);
    for j=1:l_lambdas
        x_est = l1_ls_nonneg(Ar,Ar',mr,n,yr,lambdas(j),0.001,1);
        cv_error(j) = norm(ycv-Acv*x_est,2)/mcv;
    end

    % Find best lambda from CV
    [minval,minind] = min(cv_error);

    % Solve the non-negative lasso optimization problem
    x_est = l1_ls_nonneg(A,A',m,n,y,lambdas(minind),0.001,1);
    
    % Uncomment the next two lines if you wish to switch on clustering
%     [~,centroids] = kmeans(x_est,2);
%     tau = sum(centroids)/2;

    % Classification by thresholding the output vector
    x_est(x_est < tau) = 0;
    final_x(:) = 0; final_x(nz_indices_x) = x_est;
    
    % Union with sure positives declared by Definite Defectives
    for i=1:length(nz_indices_y)
        if sum(B(i,:), 'all') == 1 && final_x(nz_indices_x(find(B(i,:)>0))) == 0
            final_x(nz_indices_x(find(B(i,:)>0))) = 10^(-8);
        end
    end
    
    % COMP NNLASSO result metrics
    rmseval(ss) = norm(x-final_x,2)/norm(x,2);
    fn(ss) = length(find(x > 0 & final_x == 0));
    fp(ss) = length(find(x == 0 & final_x > 0));
    sen(ss) = length(find(x > 0 & final_x > 0))/length(find(x > 0)); 
    spec(ss) = length(find(x == 0 & final_x == 0))/length(find(x == 0));
    fprintf (' %d %d',fn(ss),fp(ss));
end

% Print final result metrics
fprintf('\nCOMP: s = %d, [%.3f,%.3f] & [%.3f,%.3f] & [%.3f,%.3f] & [%.3f,%.3f] & [%.3f,%.3f]',...
s,mean(rmseval_comp),std(rmseval_comp), mean(fn_comp),std(fn_comp), mean(fp_comp),std(fp_comp), mean(sen_comp),std(sen_comp), mean(spec_comp),std(spec_comp));
fprintf('\nCOMP NNLASSO: s = %d, [%.3f,%.3f] & [%.3f,%.3f] & [%.3f,%.3f] & [%.3f,%.3f] & [%.3f,%.3f]',...
s,mean(rmseval),std(rmseval), mean(fn),std(fn), mean(fp),std(fp), mean(sen),std(sen), mean(spec),std(spec));
