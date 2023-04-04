
function A = SBM(BlockDiagonalBase, in_p, out_p)
    % Probability matrix for creating SBM for family-structured graph
    sz = size(BlockDiagonalBase);
    N = sz(2);
    P = in_p * BlockDiagonalBase + out_p * ones(N, N) - out_p * BlockDiagonalBase;
    A = binornd(ones([N,N]), P);
    A = tril(A);
    A = A - diag(diag(A));
    A = A + A';
end