function W = KRPPartialContractionsLR(Y,r,n)
%KRPPartialContractionsLR Computes partial contractions of cores from left to right.
%
%   W = KRPPartialContractionsLR(Y, r, n) computes partial tensor contractions 
%   from mode 1 to mode n-1 using vertically unfolded tensor cores.
%
%   Inputs:
%     - Y: A cell array of size N-by-1, where each cell contains the vertical 
%          unfolding of a tensor core. These cores are typically from a 
%          tensor-train or related decomposition.
%     - r: A positive integer specifying the target rank.
%     - n(optional): A positive integer in the range [2, N], indicating up to which 
%          core (from the left) the partial contractions are computed.
%
%   Output:
%     - W: A cell array of size (N-1)-by-1 containing the partial contractions.
%          W{i} stores the result of contracting cores from the left up to core n-1.

%% Get the number of TT-cores and mode sizes
[N,I,~] = TTsizes(Y);

%Default Parameter
if nargin < 3 
    n = N;
end

X = KRPrand(I,r);
W = cell(N-1,1);
for i = 1:1:n-1
   if i == 1
       W{i} = X{i}'*Y{i};
   else
       W{i} = khatrirao(X{i},W{i-1}')'*Y{i};
   end
end
end