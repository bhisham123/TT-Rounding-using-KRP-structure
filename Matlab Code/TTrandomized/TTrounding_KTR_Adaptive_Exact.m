function X = TTrounding_KTR_Adaptive_Exact(Y, tol, l)

    if (nargin == 2)
        l = 0;
    end
    
    [N,I,~] = TTsizes(Y);
    % 
    % T = TTNorm_Partial(Y);
    T = TTNorm_Partial_Chol(Y);
    % T = TTNorm_Partial_EigDec(Y);
    
    Y_norm = T{1};
    
    tol_local = (tol/ sqrt(N - 1)) * Y_norm ; % local tolerance
    
    b = 10; % block size
    
    W = KRPpartialContractionsRL(Y,b); 
    
    X = Y;
    % Err = [];
    for n = 1 : N - 1
        b = 10;
        reduced_cols = b;

        if n > 1
            actual_norm_sq = est_norm_sq;  
        else
            actual_norm_sq = Y_norm^2;
        end
    
    
        current_core = X{n};
        [rows,cols] = size(current_core);
        next_core = X{n+1};
    
        S = current_core*W{n}(:,1:b); 
        [X{n}, ~] = qr(S, 0);
    
    
        %compute the low rank approximation 
        M = X{n}' * current_core;
        X{n+1} =  h2v(M * v2h(next_core, I(n+1)), I(n+1));
        

        P = M*T{n+1};
        est_norm_sq = v2h(P, size(X{n},2))*v2h(M,size(X{n},2))';
        diff = sqrt(abs(actual_norm_sq - est_norm_sq));
      
        while (diff > tol_local)  
            if (reduced_cols + b) > cols 
                b = cols - reduced_cols;
                if b == 0
                    break;
                end
            end
            if size(W{n},2) < reduced_cols+b
                W_new = KRPpartialContractionsRL_New(Y,b,n);
                for j = n:N-1
                    W{j} = [W{j},W_new{j}];
                end
            end
           
            S = current_core*W{n}(:,reduced_cols+1:reduced_cols+b);
            S = S - X{n}*(X{n}'*S);
            
            [Q,~] = qr(S, 0);
            X{n} = [X{n},Q];
    
            M = Q' * current_core;
            V = M * v2h(next_core, I(n+1));
            X{n+1} =  h2v([v2h(X{n+1}, I(n+1));V], I(n+1));
            
            est_norm_sq = est_norm_sq + v2h(M*T{n+1}, b)*v2h(M,b)';
            diff = sqrt(abs(actual_norm_sq - est_norm_sq));
            reduced_cols = reduced_cols + b;
        end
        % Err(n) = diff;
    end
end
