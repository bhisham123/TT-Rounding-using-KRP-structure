function T = TTNorm_Partial_Chol(X)
    [N,I,~] = TTsizes(X);
    T = cell(N-1,1);
    % T1 = cell(N-1,1);
    Z = v2h(X{N},I(N));
    T{N} = Z*Z';
    % T1{N} = Z*Z';
    for n = N-1:-1:1
        R = cholcov(T{n+1});
        Z = v2h(X{n}*R',I(n));
        T{n} = Z*Z';

        % T1{n} = v2h(X{n}*T1{n+1},I(n))*v2h(X{n},I(n))';
        % T1{n} = (T1{n} + T1{n}')/2;
        % norm(T{n}-T1{n},'fro')  
    end
    T{1} = sqrt(abs(T{1}));  
    % T1{1} = sqrt(abs(T1{1}));  
end