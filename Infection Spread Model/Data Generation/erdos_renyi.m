function A = erdos_renyi(n, p)
    W = binornd(1,p,n,n);
    W = tril(W);
    W = W - diag(diag(W));
    A = W + W';
end
