function M = Explicit_mat(Y,I)
    N = length(Y);
    M = Y{1}*v2h(Y{2},I(2));
    for n = 3:N
        M = M*kron(v2h(Y{n},I(n)), eye(prod(I(2:n-1))));
    end
end