function X = TTrounding_Orthogonalize_then_Randomize_Adaptive(Y, tol, l)
    if (nargin == 2)
        l = 0;
    end
    
    X = TTorthogonalizeRL(Y);
    [N,I,rx] = TTsizes(X);
    
    
    normX = norm(X{1}, 'fro');
    local_tol = tol * normX /sqrt(N - 1);    % Truncation threshold
    
    % Err = [];
    
    % Rounding the tensor cores form left to right
    for n = 1 : N-1
        b = 10;
        reduced_cols = b;
        if n > 1
            actual_norm_sq = est_norm_sq;  %TTnorm(X(n:N));
        else
            actual_norm_sq = normX^2;
        end
       
        current_core = X{n};
        [rows,cols] = size(current_core);
        next_core    = X{n+1};
        
        if b  > cols 
            b = cols;
        end
    
        [X{n},~] = qr(current_core * randn(rx(n+1), b), 0);
        
        M = X{n}'*current_core;
        X{n+1} = h2v(M * v2h(next_core, I(n+1)), I(n+1));
        est_norm_sq = norm(M,'fro')^2;
        % est_norm_sq = norm(X{n+1},'fro')^2;
        diff = sqrt(abs(actual_norm_sq - est_norm_sq));
    
        while (diff > local_tol) 
            if (reduced_cols + b) > cols 
                b = cols - reduced_cols;
                if b == 0
                    break;
                end
            end
    
            S = current_core * randn(rx(n+1), b);
            S = S - X{n}*(X{n}'*S);
            % Re-orthognalization Step
            S = S - X{n}*(X{n}'*S);
            [Q,~] = qr(S, 0);
            X{n} = [X{n},Q];
            
            M = Q'*current_core;
            T = M * v2h(next_core, I(n+1));
            X{n+1} = h2v([v2h(X{n+1}, I(n+1)); T],I(n+1));
    
            est_norm_sq = est_norm_sq + norm(M,'fro')^2;
            % est_norm_sq = est_norm_sq + norm(T,'fro')^2;
    
            % err = actual_norm_sq - est_norm_sq;
            % if err < 0
            %     err
            %     size(X{n},2)
            %     size(Y{n},2)
            %     normality_error = norm(X{n}'*X{n}-eye(size(X{n},2)),'fro')
            % end
            
            diff = sqrt(abs(actual_norm_sq - est_norm_sq));
            reduced_cols = reduced_cols +b;
        end 
        % Err(n) = diff;
    end
end