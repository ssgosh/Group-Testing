function c = union_sorted(a,b)
    % This function has been reproduced from https://github.com/tminka/lightspeed
    c = sort([a(~ismember_sorted(a,b)); b]);