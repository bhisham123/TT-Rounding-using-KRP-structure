function X = TTrounding_Orthogonalize_then_Randomize_Adaptive2(Y, tol, init_f, incr_f, tol_scale, do_rounding)
%TTrounding_Orthogonalize_then_Randomize_Adaptive: Adaptive TT-rounding which first perform orthogonalization followed by random sketching.
%
%   X = TTrounding_Orthogonalize_then_Randomize_Adaptive(Y, tol, init_f, incr_f, use_svd_filtering)
%   
%   Performs adaptive Tensor-Train (TT) rounding by first perfroming orthogonalizing
%   of the input tensor from right to left followed by randomized projections from left to right
%   This method avoids explicit full SVDs by using random sketching and uses incremental QR-based 
%   basis construction.
%
%   Inputs:
%     - Y                : A TT-tensor represented as a cell array of vertically 
%                          unfolded cores (size N-by-1).
%     - tol              : Desired relative approximation tolerance (Frobenius norm).
%     - init_f           : Initial factor to determine starting block size for projection, 
%                          computed as floor(rank(n+1)/init_f). Default value is 5.
%     - incr_f           : Incremental factor used to compute the block size of
%                          incremental QR, i.e., b = floor(R(n+1)/incr_f).
%                          The default value of incr_f is 20.
% 
%   Output:
%     - X                : A TT-tensor (cell array) with reduced TT-ranks 
%                          that approximates the original tensor Y within the specified tolerance.
%



% ----------- Default Parameters --------------
if nargin < 6 || isempty(do_rounding)
    do_rounding = false; 
end
if nargin < 5 || isempty(tol_scale)
    tol_scale = 1; 
end
if nargin < 4 || isempty(incr_f)
    incr_f = 0.05;
end
if nargin < 3 || isempty(init_f)
    init_f = 0.1;
end


X = TTorthogonalizeRL(Y); %Right to left orthogonalization of cores
[N,I,rank] = TTsizes(X);


normX = norm(X{1}, 'fro'); % compute norm of tensor X

tau = tol * normX /sqrt(N - 1);  %error distributed to each core
    
% Rounding the tensor cores form left to right
for n = 1 : N-1
    current_core = X{n};
    [rows,cols] = size(current_core);
    max_rank = min(rows,cols);
    next_core    = X{n+1};

    Init_b = max(floor(max_rank*init_f),2); %initial block size
    orthogonal_cols = 0; 

    [X{n},~] = qr(current_core * randn(rank(n+1), Init_b), 0); %tall and skinny QR
    
    %compute the low rank approximation 
    M = X{n}'*current_core;
    X{n+1} = h2v(M * v2h(next_core, I(n+1)), I(n+1)); 

    orthogonal_cols = orthogonal_cols + Init_b; %update the orthogonal columns count
    
    % est_norm_sq = norm(M,'fro')^2;
    % diff = sqrt(abs(actual_norm_sq - est_norm_sq));

    b = max(floor(max_rank*incr_f),2); %block size for increment

    S = current_core * randn(rank(n+1), b); %Sketch compuation
    S = S - X{n}*(X{n}'*S); % Remove the directions already computed (orthogonalization)

    while (norm(S,'fro')/sqrt(b) > tau/tol_scale) 
        b = min(b,max_rank - orthogonal_cols); 
        if b <= 0
            break
        end

        S = S(:,1:b);
        [Q,~] = qr(S, 0); %tall and skinny QR

        %Internal Reorthogonalization
        [Q,~] = qr(Q - X{n}*(X{n}'*Q),0);

        X{n} = [X{n},Q]; %Add new directions
        
        M = Q'*current_core;
        T = M * v2h(next_core, I(n+1));
        X{n+1} = h2v([v2h(X{n+1}, I(n+1)); T],I(n+1));

        orthogonal_cols = orthogonal_cols + b;    %update the orthogonal columns count

        % % Error compuation
        % est_norm_sq = est_norm_sq + norm(M,'fro')^2;
        % diff = sqrt(abs(actual_norm_sq - est_norm_sq));

        S = current_core * randn(rank(n+1), b); %Sketch compuation
        S = S - X{n}*(X{n}'*S); % Remove the directions already computed (orthogonalization)
    end     
end
% TT-rounding of X, which is already orthogonalized.
if do_rounding
    normX = norm(X{N}, 'fro');
    tau = tol * normX / sqrt(N - 1);                    % Truncation threshold
    for n = N:-1:2
        [Q,R] = qr(v2h(X{n}, I(n))', 0);
        [U, S, V] = svd(R, 'econ');
        rank = trunc(diag(S), tau);
        X{n} = h2v( (Q*U(:, 1:rank))', I(n));
        try
        X{n-1} = X{n-1} * (V(:, 1:rank) * S(1:rank, 1:rank));
        catch
            warning('Some problem');
            break;
        end
    end
end
end