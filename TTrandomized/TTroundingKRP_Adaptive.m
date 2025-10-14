function X = TTroundingKRP_Adaptive(Y,tol,init_f,incr_f,min_samples,tol_scale)
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
%                          for projections, i.e., b_init = floor(R(n+1)*init_f).
%                          The default value of init_f is 0.1.
%     - incr_f (optional) : Incremental factor used to compute the block size of
%                          incremental QR, i.e., b = floor(R(n+1)*incr_f).
%                          The default value of incr_f is 0.05.
%     - min_samples (optional): minimum number of smaples to use for norm estimation, Default: 20.

%
%   Output:
%     - X : A TT-tensor (cell array of cores) with reduced TT-ranks, 
%           approximating the original tensor Y up to the specified tolerance.
%
% ----------- Default Parameters --------------
if nargin < 6 
    tol_scale = 1; 
end
if nargin < 5
    min_samples = 20; 
end

if nargin < 4 
    incr_f = 0.05;
end
if nargin < 3 
    init_f = 0.1;
end

[N,I,r] = TTsizes(Y);
W = cell(N-1,1);

X = Y;

sample_size = max(max(floor(r(2:N)*init_f)),min_samples);
W = KRPPartialContractionsRL(Y, sample_size,1);
normX_est = norm(X{1}*W{1},'fro')/sqrt(sample_size); % tensor norm estimation
tau = tol*normX_est/sqrt(N-1);  %error distributed to each core

%Rouding tensor cores from left to right
for n = 1 : N - 1
    Current_core = X{n};
    max_cols = min(size(Current_core));
    Next_core    = X{n+1};
    
    Init_b = max(floor(max_cols*init_f),1); %initial block size
    orthogonal_cols = 0; %track number of basis computed so far
    
    S = Current_core*W{n}(:,1:Init_b); %Compute sketch

    [X{n}, ~] = qr(S, 0); %reduced QR decomposition 
    
    %compute the low rank approximation 
    M = X{n}' * Current_core;
    X{n+1} =  h2v(M * v2h(Next_core, I(n+1)), I(n+1));
    
    orthogonal_cols = orthogonal_cols + Init_b; %update the basis count    
    
    b = max(floor(max_cols*incr_f),1); %block size for increment
    sample_size = max(b,min_samples);  %samples for norm estimation
    
    % compute partial sketches
    if size(W{n},2) < orthogonal_cols + sample_size
        W_new = KRPPartialContractionsRL(Y,(orthogonal_cols + sample_size)-size(W{n},2),n);
        for j = n:N-1
            W{j} = [W{j},W_new{j}];
        end
    end
    
    S = Current_core*W{n}(:,orthogonal_cols+1:orthogonal_cols+sample_size); %compute sketch
    S = S - X{n}*(X{n}'*S); % Remove the directions already computed (orthogonalization)
    
    while (norm(S,'fro')/sqrt(sample_size) > tau/tol_scale) 
        b = min(b,max_cols - orthogonal_cols); 
        if b <= 0
            break
        end
        
        S = S(:,1:b);
        [Q,~] = qr(S, 0); %reduced QR decomposition
        
        %Internal Reorthogonalization
        [Q,~] = qr(Q - X{n}*(X{n}'*Q),0);
    
        X{n} = [X{n},Q]; %Add new directions
    
        M = Q' * Current_core;
        V = M * v2h(Next_core, I(n+1));
        X{n+1} =  h2v([v2h(X{n+1}, I(n+1));V], I(n+1));
    
        orthogonal_cols = orthogonal_cols + b; %update the orthogonal columns count
    
        sample_size = max(b,min_samples); %samples for norm estimation

        % Compute partial sketches
        if size(W{n},2) < orthogonal_cols+sample_size
            W_new = KRPPartialContractionsRL(Y,(orthogonal_cols + sample_size)-size(W{n},2),n);
            for j = n:N-1
                W{j} = [W{j},W_new{j}];
            end
        end
        S = Current_core*W{n}(:,orthogonal_cols+1:orthogonal_cols+sample_size); %compute sketch
        S = S - X{n}*(X{n}'*S); % Remove the directions already computed (orthogonalization) 
    end
     if norm(X{n}'*X{n} - eye(size(X{n},2)),'fro') > 1e-7
        disp([' n = ', num2str(n), ', Othogonality Error ', num2str(norm(X{n}'*X{n} - eye(size(X{n},2)),'fro')), ', Condition number = ', num2str(cond(Current_core),'%.2e')])
        disp('')
     end
end
end



