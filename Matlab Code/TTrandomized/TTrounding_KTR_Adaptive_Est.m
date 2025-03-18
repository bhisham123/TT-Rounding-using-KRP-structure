function X = TTrounding_KTR_Adaptive_Est(Y,Y_norm, tol, l)
    
    if (nargin == 3)
        l = 0;
    end
    
    [N,I,~] = TTsizes(Y);
    
    tol_local = (tol/ sqrt(N - 1)) * Y_norm ; % local tolerance
    
    
    b = 10; % block size
    W = KRPpartialContractionsRL(Y,b); 
    
    X = Y;
    % Err = [];
    for n = 1 : N - 1
        b = 10;
        reduced_cols = b;

        Current_core = X{n};
        [rows,cols] = size(Current_core);
        Next_core    = X{n+1};

        if b > cols
            b = cols;
        end
    
        S = Current_core*W{n}(:,1:b); 
        [X{n}, ~] = qr(S, 0);
    
        %compute the low rank approximation 
        M = X{n}' * Current_core;
        X{n+1} =  h2v(M * v2h(Next_core, I(n+1)), I(n+1));
        
        if size(W{n},2) < reduced_cols+b
            W_new = KRPpartialContractionsRL_New(Y,b,n);
            for j = n:N-1
                W{j} = [W{j},W_new{j}];
            end
        end

        S = Current_core*W{n}(:,reduced_cols+1:reduced_cols+b);
        S = S - X{n}*(X{n}'*S);

        while (norm(S,'fro')/sqrt(b) > tol_local) 
            if (reduced_cols + b) > cols 
                b = cols - reduced_cols;
                if b == 0
                    break;
                end
            end
            S = S(:,1:b);
            %Re-orthogonalization step
            S = S - X{n}*(X{n}'*S);

            [Q,~] = qr(S, 0);
            X{n} = [X{n},Q];

            M = Q' * Current_core;
            V = M * v2h(Next_core, I(n+1));
            X{n+1} =  h2v([v2h(X{n+1}, I(n+1));V], I(n+1));
            
            reduced_cols = reduced_cols +b;
            if size(W{n},2) < reduced_cols+b
                W_new = KRPpartialContractionsRL_New(Y,b,n);
                for j = n:N-1
                    W{j} = [W{j},W_new{j}];
                end
            end
            S = Current_core*W{n}(:,reduced_cols+1:reduced_cols+b);
            S = S - X{n}*(X{n}'*S);
        end
        % Err(n) = norm(S,'fro')/sqrt(b);
    end
end
