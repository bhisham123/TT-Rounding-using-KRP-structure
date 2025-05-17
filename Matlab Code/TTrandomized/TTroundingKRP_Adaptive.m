function X = TTroundingKRP_Adaptive(Y,tol, init_f,incr_f)
%TTroundingKRP_Adaptive: Adaptive TT-rounding using Khatri-Rao-based projections.
%
%   X = TTroundingKRP_Adaptive(Y, tol, init_f, incr_f, use_svd_filtering)
%   performs adaptive Tensor-Train (TT) rounding on a TT-format tensor Y 
%   using randomized projections based on the Khatri-Rao product and an 
%   incremental block expansion strategy.
%
%   Inputs:
%     - Y                 : A TT-tensor represented as a cell array of cores 
%                          (each vertically unfolded). Size N-by-1.
%     - tol               : Desired approximation accuracy in Frobenius norm.
%     - init_f (optional) : Initial factor to compute the starting block size 
%                          for projections, i.e., b_init = floor(R(n+1)/init_f).
%                          The default value of init_f is 5.
%     - incr_f (optional) : Incremental factor used to compute the block size of
%                          incremental QR, i.e., b = floor(R(n+1)/incr_f).
%                          The default value of incr_f is 20.
%
%
%   Output:
%     - X : A TT-tensor (cell array of cores) with reduced TT-ranks, 
%           approximating the original tensor Y up to the specified tolerance.
%
%   Description:
%     The function iteratively projects each TT core onto an orthonormal basis
%     constructed using randomized range finding. It begins with a small number
%     of random projections determined by `init_f`, and adaptively increases the 
%     projection dimension using `incr_f` until the approximation error falls 
%     below a locally computed tolerance.
%
%     The contraction matrices used in the Khatri-Rao product computations are 
%     updated in blocks using 'KRPpartialContractionsRL_New' when needed.

% ----------- Default Parameters --------------
if nargin < 4 || isempty(incr_f)
    incr_f = 20;
end
if nargin < 3 || isempty(init_f)
    init_f = 5;
end

[N,I,r] = TTsizes(Y);
W = cell(N-1,1);

X = Y;
%Rouding tensor cores from left to right
for n = 1 : N - 1
    Init_b = max(floor(r(n+1)/init_f),2); %initial block size
    orthogonal_cols = 0; 
    
    Current_core = X{n};
    [~,cols] = size(Current_core);
    Next_core    = X{n+1};
    
    if size(W{n},2) < Init_b
         W_new = KRPPartialContractionsRL(Y,Init_b-size(W{n},2),n);
        for j = n:N-1
            W{j} = [W{j},W_new{j}];
        end
    end
    
    S = Current_core*W{n}(:,1:Init_b); %Compute sketch
    
    if n ==1
        normX_est = norm(S,'fro')/sqrt(Init_b); % tensor norm estimation
        tau = tol*normX_est/sqrt(N-1);  %error distributed to each core
    end
    
    [X{n}, ~] = qr(S, 0); %tall and skinny QR 
    
    %compute the low rank approximation 
    M = X{n}' * Current_core;
    X{n+1} =  h2v(M * v2h(Next_core, I(n+1)), I(n+1));
    
    orthogonal_cols = orthogonal_cols + Init_b; %update the orthogonal columns count
    
    b = max(floor(r(n+1)/incr_f),2); %block size for increment
    
    % compute partial sketches
    if size(W{n},2) < orthogonal_cols + b
        W_new = KRPPartialContractionsRL(Y,b,n);
        for j = n:N-1
            W{j} = [W{j},W_new{j}];
        end
    end
    
    S = Current_core*W{n}(:,orthogonal_cols+1:orthogonal_cols+b); %compute sketch
    S = S - X{n}*(X{n}'*S); % Remove the directions already computed (orthogonalization)
    
    while (norm(S,'fro')/sqrt(b) > tau) 
        b = min(b,cols - orthogonal_cols); 
        if b <= 0
            break
        end
        
        S = S(:,1:b);
        [Q,~] = qr(S, 0); %tall and skinny QR
        
        %Internal Reorthogonalization
        [Q,~] = qr(Q - X{n}*(X{n}'*Q),0);
    
        X{n} = [X{n},Q]; %Add new directions
    
        M = Q' * Current_core;
        V = M * v2h(Next_core, I(n+1));
        X{n+1} =  h2v([v2h(X{n+1}, I(n+1));V], I(n+1));
    
        orthogonal_cols = orthogonal_cols + b; %update the orthogonal columns count
    
        % Compute partial sketches
        if size(W{n},2) < orthogonal_cols+b
            W_new = KRPPartialContractionsRL(Y,b,n);
            for j = n:N-1
                W{j} = [W{j},W_new{j}];
            end
        end
        S = Current_core*W{n}(:,orthogonal_cols+1:orthogonal_cols+b); %compute sketch
        S = S - X{n}*(X{n}'*S); % Remove the directions already computed (orthogonalization)
    end
     % if norm(X{n}'*X{n} - eye(size(X{n},2)),'fro') > 1e-7
     %    disp([' n = ', num2str(n), ', Othogonality Error ', num2str(norm(X{n}'*X{n} - eye(size(X{n},2)),'fro')), ', Condition number = ', num2str(cond(Current_core),'%.2e')])
     %    disp('')
     % end
end
end



