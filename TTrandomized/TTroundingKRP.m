function X = TTroundingKRP(Y, l)
%TTroundingKRP Performs TT-rounding using Khatri-Rao Product (KRP) contractions.
%
%   X = TTroundingKRP(Y, l) compresses a tensor in Tensor Train (TT) format
%   by rounding its ranks using a Khatri-Rao Product-based sketching approach.
%
%   INPUT:
%       Y : Cell array of size N x 1 containing the TT-cores of the input tensor.
%           Each Y{n} is a 2D array representing the vertical unfolding of the n-th TT-core.
%       l : Target rank(s). Can be either:
%               - A scalar (applied uniformly to all TT ranks),
%               - A vector of length N specifying rank at each position.
%
%   OUTPUT:
%       X : Cell array of size N x 1 containing the rounded TT-cores.


% Get the number of TT-cores and mode sizes
[N, I, ~] = TTsizes(Y);
    
% Compute right-to-left partial contractions using KRP skecthing up to rank max(l)
W = KRPPartialContractionsRL(Y,max(l));

% Initialize output cell array
X = cell(N,1);

% Left-to-right orthogonalization and TT-rounding
for n = 1 : N - 1
    % Determine target rank for current step
    if length(l) == 1
        rank = l;
    else
        rank = l(n + 1);
    end
    
    % QR decomposition to obtain orthonormal core X{n}
    [X{n}, ~] = qr(Y{n}*W{n}(:,1:rank), 0);

    % Project Y{n} onto X{n} and update next core
    M = X{n}' * Y{n};
    Y{n + 1} = h2v(M * v2h(Y{n + 1}, I(n+1)), I(n+1));
end
% Update Final TT-core
X{N} = Y{N};
end
