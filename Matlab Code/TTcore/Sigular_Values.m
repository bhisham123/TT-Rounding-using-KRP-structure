function Sing_values = Sigular_Values(Y,flag)

    if nargin < 2
        flag = true;  
    end
    
    
    X = TTorthogonalizeRL(Y);
    [N,I,~] = TTsizes(X);
    
     
    Sing_values = cell(N-1,1);
    for n = 1:N-1
        % [X{n},R] = qr(X{n}, 0);
        [U, S, V] = svd(X{n}, 'econ');
        Sing_values{n} = diag(S);
        X{n} = U;
        X{n+1} = S * V' * v2h(X{n+1}, I(n+1));
        X{n+1} = h2v(X{n+1}, I(n+1));
    end
    


end