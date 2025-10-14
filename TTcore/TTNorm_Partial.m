function T = TTNorm_Partial(X)
    [N,I,~] = TTsizes(X);
    T = cell(N-1,1);
    T{N} = v2h(X{N},I(N)) * v2h(X{N},I(N))';
    for n = N-1:-1:1
        T{n} = v2h(X{n}*T{n+1},I(n))*v2h(X{n},I(n))';
        T{n} = (T{n} + T{n}')/2;
    end
    T{1} = sqrt(abs(T{1}));  
end