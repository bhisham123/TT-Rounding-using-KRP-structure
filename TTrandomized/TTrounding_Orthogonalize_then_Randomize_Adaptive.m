function X = TTrounding_Orthogonalize_then_Randomize_Adaptive(Y, tol, init_f, incr_f)
%TTrounding_Orthogonalize_then_Randomize_Adaptive: Adaptive TT-rounding which first perform orthogonalization followed by random sketching.
%
%   X = TTrounding_Orthogonalize_then_Randomize_Adaptive(Y, tol, init_f, incr_f)
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
if nargin < 4
    incr_f = 0.05;
end
if nargin < 3
    init_f = 0.1;
end


X = TTorthogonalizeRL(Y); %Right to left orthogonalization of cores
[N,I,rank] = TTsizes(X);


normX = norm(X{1}, 'fro'); % compute norm of tensor X

tau = tol * normX /sqrt(N - 1);  %error distributed to each core
    
% Rounding the tensor cores form left to right
for n = 1 : N-1
    if n > 1
        FullNormSq = ProjNormSq;  
    else
        FullNormSq = normX^2;
    end
    current_core = X{n};
    [rows,cols] = size(current_core);
    max_rank = min(rows,cols);
    next_core    = X{n+1};

    % Init_b = max(floor(max_rank*init_f),1); %initial block size
    Init_b = 1;
    orthogonal_cols = 0; 

    [X{n},~] = qr(current_core * randn(rank(n+1), Init_b), 0); %tall and skinny QR
    
    %compute the low rank approximation 
    M = X{n}'*current_core;
    X{n+1} = h2v(M * v2h(next_core, I(n+1)), I(n+1)); 
    
    orthogonal_cols = orthogonal_cols + Init_b; %update the orthogonal columns count
    
    ProjNormSq = norm(M,'fro')^2;
   

    % b = max(floor(max_rank*incr_f),1); %block size for increment
    b = 1;

    while (sqrt(FullNormSq - ProjNormSq) > tau) 
        b = min(b,max_rank - orthogonal_cols); 
        if b <= 0
            break
        end
        S = current_core * randn(rank(n+1), b);
        S = S - X{n}*(X{n}'*S);
        [Q,~] = qr(S, 0); %tall and skinny QR

        %Internal Reorthogonalization
        [Q,~] = qr(Q - X{n}*(X{n}'*Q),0);

        X{n} = [X{n},Q]; %Add new directions
        
        M = Q'*current_core;
        T = M * v2h(next_core, I(n+1));
        X{n+1} = h2v([v2h(X{n+1}, I(n+1)); T],I(n+1));

        orthogonal_cols = orthogonal_cols + b;    %update the orthogonal columns count
        
        ProjNormSq = ProjNormSq + norm(M,'fro')^2;
    end     
end
end