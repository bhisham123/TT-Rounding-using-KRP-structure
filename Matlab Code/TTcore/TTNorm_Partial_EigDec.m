function T = TTNorm_Partial_EigDec(X)
    [N,I,~] = TTsizes(X);
    T = cell(N-1,1);
    % T2 = cell(N-1,1);
    Z = v2h(X{N},I(N));
    T{N} = Z*Z';
    % T2{N} = Z*Z';
    for n = N-1:-1:1
        [V,E] = eig(T{n+1});
        BigE = abs(diag(E))>(1e-16*abs(E(end)));
        R = V(:,BigE)*sqrt(E(BigE,BigE));
        % size(R)
        Z = v2h(X{n}*R,I(n));
        T{n} = Z*Z';
        % T{n} = (T{n} + T{n}')/2;
        norm(T{n+1}-(V*E*V'),'fro')/norm(T{n+1},'fro')
        % AA = v2h(X{n}*(V*E*V'),I(n))*v2h(X{n},I(n))';
        % T{n} = v2h(X{n}*T{n+1},I(n))*v2h(X{n},I(n))';
        % norm(AA - T{n},'fro')
        % T{n} = (T{n} + T{n}')/2;
        % figure, semilogy(eig(T{n}))
        % T2{n} = v2h(X{n}*T2{n+1},I(n))*v2h(X{n},I(n))';
        % T2{n} = (T2{n} + T2{n}')/2;
        % norm(T{n}-T2{n},'fro')
    end
    T{1} = sqrt(abs(T{1}));
    % T2{n} = sqrt(abs(T2{1}));
end