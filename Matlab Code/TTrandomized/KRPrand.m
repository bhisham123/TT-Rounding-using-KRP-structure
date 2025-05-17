function X = KRPrand(I,r,n)
%KRPRrand Generates random matrices with normal distribution for partial Khatri-Rao product contractions.
%
%   X = KRPRrand(I, R, n) returns a cell array of random matrices used for
%   partial contractions from mode N down to mode n+1.
%
%   Inputs:
%     - I: A vector of length N specifying the mode sizes of a tensor.
%     - r: A scalar specifying the ranks. 
%     - n (optional): A positive integer in the range [1, N-1], indicating the point up
%          to which the matrices should be generated (i.e., for modes N to n+1).
%
%   Output:
%     - X: A cell array of size (N-1)-by-1 where X{k} contains a random matrix
%          of size I(k+1) × r, for k = N-1 down to n.


if nargin < 3 
    n = 1;
end

N = length(I);
X = cell(N-1,1);
for k = N:-1:n+1
    X{k-1} = randn(I(k), r);
end
end