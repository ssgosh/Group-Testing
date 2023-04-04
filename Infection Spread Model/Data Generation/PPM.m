% Stochastic block model
% n nodes
% k communities, of size n / k each.
% p = intra-community link prob
% q = inter-community link prob

% Following generates SBM with same expected number of edges as the family
% graph with CT, accumulating 20 days before day 38.
% imshow(SBM(1000,20,0.1228,0.0016))
function A = PPM(n, k, p, q)
    s = n / k;
    % Set parameters for Bernoulli trials
    P = q * ones([n,n]);
    for i = 1:k
        rows = s*(i-1)+1:s*i;
        P(rows, rows) = p;
    end

    A = binornd(ones([n,n]), P);
    A = tril(A);
    A = A - diag(diag(A));
    A = A + A';
end