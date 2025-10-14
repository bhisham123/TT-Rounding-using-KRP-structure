function W = KRPPartialContractionsRL(Y,r,n,sphere)
%KRPPartialContractionsRL Computes partial contractions of cores from right to left.
%
%   W = KRPPartialContractionsRL(Y, r, n) computes partial tensor contractions 
%   from mode N to mode n+1 using vertically unfolded tensor cores.
%
%   Inputs:
%     - Y: A cell array of size N-by-1, where each cell contains the vertical 
%          unfolding of a tensor core. These cores are typically from a 
%          tensor-train or related decomposition.
%     - r: A positive integer specifying the target rank.
%     - n(optional): A positive integer in the range [1, N-1], indicating up to which 
%          core (from the right) the partial contractions are computed.
%     - sphere (optional): Binary variable. If True use spherical distribution oteherwise use Gaussian.
%                             Default value is false. 
%   Output:
%     - W: A cell array of size (N-1)-by-1 containing the partial contractions.
%          W{i} stores the result of contracting cores from the right up to core n+1.

% ----------- Default Parameters --------------

if nargin < 4
    sphere = false;
end
if nargin < 3 
    n = 1;
end


[N,I,~] = TTsizes(Y);

if isscalar(r)
    l = ones(N+1, 1);
    l(2:N) = r;
else
    l = r;
end

X = KRPrand(I,l,n,sphere);
W = cell(N-1,1);
for i = N:-1:n+1
   temp_r = max(l(2:i));
   if i == N
       W{i-1} = v2h(Y{i},I(i))*(X{i-1}(:,1:temp_r));
   else
       W{i-1} = v2h(Y{i},I(i))*khatrirao(W{i}(:,1:temp_r),X{i-1}(:,1:temp_r));
   end
end
end