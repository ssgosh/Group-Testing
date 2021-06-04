current_path = pwd;

% Replace the path in the next line with a path to your cvx folder
cd ~/cvx;
cvx_setup;
cd(current_path);

addpath('../../Testing Algorithms/Helper Functions');

mat_prefix = '../../Pooling Matrix Design/Balanced Binary Matrices/Examples/psi_opt_bal_';
matrix_names = [mat_prefix + "375x1000"; mat_prefix + "300x1000"; mat_prefix + "250x1000";
    mat_prefix + "200x1000"; mat_prefix + "150x1000"; mat_prefix + "100x1000"];
matrix_name = string(matrix_names(mat_index,:));

load('../../Infection Spread Model/Data/ct_data_' + string(dat_index) + '.mat');

fid = fopen('../Results/CTR-' + string(dat_index) + string(mat_index) + '.txt','a');

x_est_lasso = cell(50,1);
x_est_glasso = cell(50,1);
x_est_oglasso = cell(50,1); 

for ss=1:50
    k = maxind-25+ss;
    
    tau = 0.2;
    
    load(matrix_name + '.mat','B'); % load matrix
    B = full(B);
    CCT = CT{k};
    for km=max(k-6,1):(k-1)
        CCT = CCT|CT{km};
    end
    [m,n] = size(B); % number of measurements, number of signal elements
    x = X{k};
    
    totalpositive(ss) = nnz(x);
    
    sigval = 0.1;
    q = 0.95;
    eps = sigval*randn(m,1);
    y = (B*x).*((1+q).^eps); % noisy measurement generation
    
    % COMP algorithm
    final_x = ones(n,1);
    nz_indices_y = find(y > 0); % measurements which were nonzero
    z_indices_y = find(y==0); % measurements which were zero-valued
    for i=1:length(z_indices_y)
       zi = z_indices_y(i);
       final_x(find(B(zi,:)>0)) = 0;
    end
    B(:,final_x == 0) = []; % remove entries of x which are obvious negatives
    nz_indices_x = find(final_x ~=  0); % store indices of putative positives
    y = y(nz_indices_y); % remove entries of y corresponding to zero-valued measurements
    B = B(nz_indices_y,:); % remove rows of A corresponding to zero-valued measurements
    [m,n] = size(B); % updated size of A
    fprintf(fid,'\nSignal %d: %d %d;', ss, m, n);
    
    % Updating groups
    groups = [];
    count = 0;
    numgroups = 0;
    for i=1:ngroups
        count1 = 0;
        for j=prgroups(1,i):prgroups(2,i)
            if ismember(j, nz_indices_x)
                count1 = count1 + 1;
            end
        end
        if count1 > 0
            groups = [groups [count + 1 count + count1]'];
            numgroups = numgroups + 1;
            count = count + count1;
        end
    end
    
    CCT = CCT(nz_indices_x, nz_indices_x);
    Cliques = maximalCliques(CCT,'v2');
    [~,numCliques] = size(Cliques);
    
    CPeople = Cliques(:,1);
    for i=2:numCliques
        CPeople = CPeople|Cliques(:,i);
    end
    for i=1:n
        if CPeople(i) == 0
            Cliques = [Cliques double(1:n == i)'];
        end
    end
    [~,numCliques] = size(Cliques);
    
    FinCliques = cell(0);
    for i=1:numCliques
        clique = find(Cliques(:,i) ==  1);
        if length(clique) ~= 0
            FinCliques = [FinCliques; clique];
        end
    end
    [numCliques,~] = size(FinCliques);

    % COMP result metrics
    rmseval_comp(ss) = norm(x-final_x,2)/norm(x,2);
    fn_comp(ss) = length(find(x > 0 & final_x == 0));
    fp_comp(ss) = length(find(x == 0 & final_x > 0));
    sen_comp(ss) = length(find(x > 0 & final_x > 0))/length(find(x > 0)); 
    spec_comp(ss) = length(find(x == 0 & final_x == 0))/length(find(x == 0)); 
    fprintf(fid,' %d %d;',fn_comp(ss),fp_comp(ss));
    
    % Stuff for cross-validation
    mr = ceil(0.9*m); mcv = m-mr;
    indices = randperm(m);
    yr = y(indices(1:mr));
    ycv = y(indices(mr+1:end));
    Br = B(indices(1:mr),:);
    Bcv = B(indices(mr+1:end),:);
    lambda_min = max([sigval*sqrt(log(n))-8,0.001]);
    lambda_max = sigval*sqrt(log(n))+8;
    lambdas = lambda_min:0.05:lambda_max;
    l_lambdas = length(lambdas);
    
    Bp = B(:,FinCliques{1});
    for i=2:numCliques
        Bp = [Bp B(:,FinCliques{i})];
    end
    Brp = Br(:,FinCliques{1});
    for i=2:numCliques
        Brp = [Brp Br(:,FinCliques{i})];
    end
    Bcvp = Bcv(:,FinCliques{1});
    for i=2:numCliques
        Bcvp = [Bcvp Bcv(:,FinCliques{i})];
    end
    [~,np] = size(Bp);
    
    tic
    
    % Cross Validation
    cv_error = length(l_lambdas);
    for i=1:l_lambdas
        cvx_begin quiet
            variable x_est(np)
            expression xs(numCliques)
            count = 0;
            for j=1:numCliques
                xs(j) = norm(x_est(count+1:count+length(FinCliques{j})), 2);
                count = count+length(FinCliques{j});
            end
            minimize( norm(yr-Brp*x_est, 2) + lambdas(i)*sum(xs) )
            subject to
                x_est >= 0;
        cvx_end
        cv_error(i) = norm(ycv-Bcvp*x_est,2)/mcv;
    end

    % Find best lambda from CV
    [~, minind] = min(cv_error);
    
    % Solve the sparse group lasso optimization problem
    cvx_begin quiet
        variable x_est(np)
        expression xs(numCliques)
        count = 0;
        for j=1:numCliques
            xs(j) = norm(x_est(count+1:count+length(FinCliques{j})), 2);
            count = count+length(FinCliques{j});
        end
        minimize( norm(y-Bp*x_est, 2) + lambdas(minind)*sum(xs) )
        subject to
            x_est >= 0;
    cvx_end
    
    count = 0;
    w_est = zeros(n,1);
    for j=1:numCliques
        x_est_group = zeros(n,1);
        x_est_group(FinCliques{j}) = x_est(count+1:count+length(FinCliques{j}));
        w_est = w_est + x_est_group;
        count = count+length(FinCliques{j});
    end
    x_est = w_est;
    
    % COMP OGLASSO result metrics
    x_est(x_est < tau) = 0;
    final_x(:) = 0; final_x(nz_indices_x) = x_est;
    for i=1:length(nz_indices_y)
        if sum(B(i,:), 'all') == 1 && final_x(nz_indices_x(find(B(i,:)>0))) == 0
            final_x(nz_indices_x(find(B(i,:)>0))) = y(i);
        end
    end
    x_est_oglasso{ss} = final_x;
    rmseval_oglasso(ss) = norm(x-final_x,2)/norm(x,2);
    fn_oglasso(ss) = length(find(x > 0 & final_x == 0));
    fp_oglasso(ss) = length(find(x == 0 & final_x > 0));
    sen_oglasso(ss) = length(find(x > 0 & final_x > 0))/length(find(x > 0)); 
    spec_oglasso(ss) = length(find(x == 0 & final_x == 0))/length(find(x == 0)); 
    fprintf(fid,' %d %d;', fn_oglasso(ss), fp_oglasso(ss));
    
    toc
    
    tic
    
    % Cross Validation
    cv_error = length(l_lambdas);
    for i=1:l_lambdas
        cvx_begin quiet
            variable x_est(n)
            expression xs(numgroups)
            for j=1:numgroups
                xs(j) = norm(x_est(groups(1,j):groups(2,j)), 2);
            end
            minimize( norm(yr-Br*x_est, 2) + lambdas(i)*sum(xs) )
            subject to
                x_est >= 0;
        cvx_end
        cv_error(i) = norm(ycv-Bcv*x_est,2)/mcv;
    end
    
    % Find best lambda from CV
    [~, minind] = min(cv_error);

    % Solve the sparse group lasso optimization problem
    cvx_begin quiet
        variable x_est(n)
        expression xs(numgroups)
        for j=1:numgroups
            xs(j) = norm(x_est(groups(1,j):groups(2,j)), 2);
        end
        minimize( norm(y-B*x_est, 2) + lambdas(minind)*sum(xs) )
        subject to
            x_est >= 0;
    cvx_end
    
    % COMP GLASSO result metrics
    x_est(x_est < tau) = 0;
    final_x(:) = 0; final_x(nz_indices_x) = x_est;
    for i=1:length(nz_indices_y)
        if sum(B(i,:), 'all') == 1 && final_x(nz_indices_x(find(B(i,:)>0))) == 0
            final_x(nz_indices_x(find(B(i,:)>0))) = y(i);
        end
    end
    x_est_glasso{ss} = final_x;
    rmseval_glasso(ss) = norm(x-final_x,2)/norm(x,2);
    fn_glasso(ss) = length(find(x > 0 & final_x == 0));
    fp_glasso(ss) = length(find(x == 0 & final_x > 0));
    sen_glasso(ss) = length(find(x > 0 & final_x > 0))/length(find(x > 0)); 
    spec_glasso(ss) = length(find(x == 0 & final_x == 0))/length(find(x == 0)); 
    fprintf(fid,' %d %d;', fn_glasso(ss), fp_glasso(ss));
    
    toc
    
    tic

    % Cross Validation
    cv_error = length(l_lambdas);
    for i=1:l_lambdas
        x_est = l1_ls_nonneg(Br,Br',mr,n,yr,lambdas(i),0.001,1);
        cv_error(i) = norm(ycv-Bcv*x_est,2)/mcv;
    end

    % Find best lambda from CV
    [~,minind] = min(cv_error);

    % Solve the non-negative lasso optimization problem
    x_est = l1_ls_nonneg(B,B',m,n,y,lambdas(minind),0.001,1);

    % COMP NNLASSO result metrics
    x_est(x_est < tau) = 0;
    final_x(:) = 0; final_x(nz_indices_x) = x_est;
    for i=1:length(nz_indices_y)
        if sum(B(i,:), 'all') == 1 && final_x(nz_indices_x(find(B(i,:)>0))) == 0
            final_x(nz_indices_x(find(B(i,:)>0))) = y(i);
        end
    end
    x_est_lasso{ss} = final_x;
    rmseval_lasso(ss) = norm(x-final_x,2)/norm(x,2);
    fn_lasso(ss) = length(find(x > 0 & final_x == 0));
    fp_lasso(ss) = length(find(x == 0 & final_x > 0));
    sen_lasso(ss) = length(find(x > 0 & final_x > 0))/length(find(x > 0));
    spec_lasso(ss) = length(find(x == 0 & final_x == 0))/length(find(x == 0));
    fprintf(fid,' %d %d',fn_lasso(ss),fp_lasso(ss));
    
    toc
end

% Print final result metrics
fprintf(fid,'\nNumber of infected people: %.5f, %.5f',mean(totalpositive),std(totalpositive));
fprintf(fid,'\nCOMP: [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f]',...
mean(rmseval_comp),std(rmseval_comp), mean(fn_comp),std(fn_comp), mean(fp_comp),std(fp_comp), mean(sen_comp),std(sen_comp), mean(spec_comp),std(spec_comp));
fprintf(fid,'\nCOMP SQRT OVERLAPPING GROUP LASSO: [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f]',...
mean(rmseval_oglasso),std(rmseval_oglasso), mean(fn_oglasso),std(fn_oglasso), mean(fp_oglasso),std(fp_oglasso), mean(sen_oglasso),std(sen_oglasso), mean(spec_oglasso),std(spec_oglasso));
fprintf(fid,'\nCOMP SQRT GROUP LASSO: [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f]',...
mean(rmseval_glasso),std(rmseval_glasso), mean(fn_glasso),std(fn_glasso), mean(fp_glasso),std(fp_glasso), mean(sen_glasso),std(sen_glasso), mean(spec_glasso),std(spec_glasso));
fprintf(fid,'\nCOMP LASSO: [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f] & [%.5f,%.5f]',...
mean(rmseval_lasso),std(rmseval_lasso), mean(fn_lasso),std(fn_lasso), mean(fp_lasso),std(fp_lasso), mean(sen_lasso),std(sen_lasso), mean(spec_lasso),std(spec_lasso));

fclose(fid);

save('../Results/ctr-' + string(dat_index) + string(mat_index) + '.mat');

exit;