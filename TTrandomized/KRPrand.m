function X = KRPrand(I,r,n,sphere)
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
%     - sphere (optional): Binary variable. True means use spherical distribution and false means Gaussian.
%                             Default value is false.
%
%   Output:
%     - X: A cell array of size (N-1)-by-1 where X{k} contains a random matrix
%          of size I(k+1) × r, for k = N-1 down to n.


if nargin < 3 
    n = 1;
end
if nargin < 4 
    sphere = false;
end


N = length(I);
if isscalar(r)
    l = ones(N+1, 1);
    l(2:N) = r;
else
    l = r;
end


X = cell(N-1,1);
for k = N:-1:n+1
    temp_r = max(l(2:k));
    X{k-1} = randn(I(k), temp_r);
    if sphere
            X{k-1} = X{k-1} ./ vecnorm(X{k-1}, 2, 1);
            X{k-1} = X{k-1} * sqrt(I(k));
    end
end
end