function [G] = nestedContractionsRL(X, Y)
% This subroutine computes the sequence of contractions between the TT 
% tensor $X$ and the subparts of the TT tensor $Y$ from Right to Left

% Function handle Horizontal to Vertical unfolding and vice versa
H2V =@(X, I) reshape(X, [size(X, 1) * I, size(X, 2)/I]);
V2H =@(X, I) reshape(X, [size(X, 1)/I, I * size(X, 2)]);

N = size(X, 1); % number of modes
rx = zeros(N + 1, 1); % TT ranks of X
rx(1) = 1; % Left boundary condition
rx(N + 1) = 1; % Right boundary condition

ry = zeros(N + 1, 1); % TT ranks of Y
ry(1) = 1; % Left boundary condition
ry(N + 1) = 1; % Right boundary condition

I = zeros(N, 1); % modes sizes
I(1) = size(X{1}, 1);

for i = 2 : N
    rx(i) = size(X{i - 1}, 2);
    ry(i) = size(Y{i - 1}, 2);
    I(i) = size(X{i}, 1)/rx(i);
end

G = cell(N, 1);
% Note that G{1} will be the inner product between X and Y

% forming G_n = contractions between cores (n+1, \ldots, N) of X and Y 
G{N} = V2H(X{N}, I(N)) * V2H(Y{N}, I(N))';
for i = N - 1 : -1 : 1
    temp = X{i} * G{i + 1};
    G{i} = V2H(temp, I(i))* V2H(Y{i}, I(i))';
end

