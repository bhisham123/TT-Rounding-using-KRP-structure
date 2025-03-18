function X = TTrounding_KTR(Y, l)

    [N,I,~] = TTsizes(Y);
    W = KRPpartialContractionsRL(Y,max(l));
    X = cell(N,1);
    for n = 1 : N - 1
        if ( length(l) == 1)
            rank = l;
        else
            rank = l(n+1);
        end
        
        
        [X{n}, ~] = qr(Y{n}*W{n}(:,1:rank), 0);    
        M = X{n}' * Y{n};
        Y{n + 1} = h2v(M * v2h(Y{n + 1}, I(n+1)), I(n+1));
    end
    X{N} = Y{N};
end
