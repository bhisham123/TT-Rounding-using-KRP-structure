function X = KTRrand(I, R)

N = length(I);
if( length(R) == 1)
    R = R * ones(N -1, 1);
end

assert(length(R) == N-1, "Please provide conforming sizes of the modes sizes and the TT ranks");
% assert(R(1) == 1 && R(N+1) == 1, "Boundary TT ranks must be 1");

X = cell(N-1, 1);
for k = N:-1:2
    % X{k-1} = 1/sqrt(I(k)*R(k-1))*randn(I(k), R(k-1));
    X{k-1} = randn(I(k), R(k-1));
end