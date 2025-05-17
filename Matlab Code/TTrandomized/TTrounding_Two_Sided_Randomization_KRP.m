function X = TTrounding_Two_Sided_Randomization_KRP(Y, l)
%TTrounding_Two_Sided_Randomization_KRP Performs randomized TT-rounding using Khatri-Rao based sketching from both sides.
%
%   X = TTrounding_Two_Sided_Randomization_KRP(Y, l) compresses a tensor 
%   in Tensor Train (TT) format by applying a two-sided randomized 
%   rounding scheme based on Khatri-Rao product sketching.
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

%% Get the number of TT-cores and mode sizes
[N,I,~] = TTsizes(Y);
    
% Heuristic choice of oversampling ranks
rho = [1; floor(2*l(2:N)'); 1];
    
% Compute Khatri-Rao partial contractions from both sides
WL = KRPPartialContractionsLR(Y,max(l)); 
WR = KRPPartialContractionsRL(Y, max(rho));
    
% Initialize output tensor cores
X = cell(N,1);

% First core
[U,S,V] = svd(WL{1}(1:l(2),:) * WR{1}(:,1:rho(2)), 'econ');
S = diag(diag(S).^(-0.5));
L = WR{1}(:,1:rho(2))*V*S;
R = S*U'*WL{1}(1:l(2),:);
X{1} = Y{1} * L;

% Middle cores
for n=2:N-1
   [U,S,V] = svd(WL{n}(1:l(n+1),:) * WR{n}(:,1:rho(n+1)), 'econ');
   S = diag(diag(S).^(-0.5));
   L = WR{n}(:,1:rho(n+1))*V*S;
   X{n} = h2v(R * v2h(Y{n} * L, I(n)), I(n));
   R = S*U'*WL{n}(1:l(n+1),:);
end
% Last core
X{N} = h2v( R * v2h(Y{N}, I(N)), I(N));
end